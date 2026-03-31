# ==============================================================================
# run_cal.R — Main calibration script (generic for any phase)
#
# USAGE:
#   Rscript run_cal.R config/phase1.yaml       # from terminal
#
#   PHASE_CONFIG <- "config/phase2.yaml"        # from RStudio/VSCode
#   source("run_cal.R")
# ==============================================================================

.bak_config <- if (exists("PHASE_CONFIG")) PHASE_CONFIG else NULL
.bak_master <- if (exists(".MASTER_STATE")) .MASTER_STATE else NULL
rm(list = ls())
gc()
if (!is.null(.bak_config)) PHASE_CONFIG <- .bak_config
if (!is.null(.bak_master)) .MASTER_STATE <- .bak_master
rm(.bak_config, .bak_master)

library(SWATrunR)
library(yaml)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(hydroGOF)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. READ CONFIG
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (exists("PHASE_CONFIG")) {
  config_path <- PHASE_CONFIG
} else if (length(args) >= 1) {
  config_path <- args[1]
} else {
  stop("Usage: Rscript run_cal.R config/phaseX.yaml\n  Or set PHASE_CONFIG before source().")
}

cfg <- yaml::read_yaml(config_path)
cat("=== CAL Phase", cfg$phase, "- Type:", cfg$cal$type, "===\n")

# ------------------------------------------------------------------------------
# 2. LOAD MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/param_inherit.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. VERSIONING
# ------------------------------------------------------------------------------
ver <- create_versioned_folder(cfg$project_path, cfg$cal$prefix)
cat("TAG:", ver$tag, "\nFOLDER:", ver$path, "\n")

# ------------------------------------------------------------------------------
# 4. LOAD GSA RESULTS (Morris ranking)
# ------------------------------------------------------------------------------
gsa_results <- find_gsa_results(
  cfg$project_path,
  cfg$gsa$prefix,
  cfg$phase
)

cat("GSA results from:", gsa_results$folder, "\n")

# ------------------------------------------------------------------------------
# 5. PARAMETERS
# ------------------------------------------------------------------------------
if (cfg$phase == 1) {
  # Phase 1: parameters from YAML + extras for calibration
  gsa_params <- tibble(
    short    = sapply(cfg$parameters, `[[`, "short"),
    par_name = sapply(cfg$parameters, `[[`, "par_name"),
    min      = sapply(cfg$parameters, `[[`, "min"),
    max      = sapply(cfg$parameters, `[[`, "max")
  )

  if (!is.null(cfg$cal_extra_params)) {
    extras <- tibble(
      short    = sapply(cfg$cal_extra_params, `[[`, "short"),
      par_name = sapply(cfg$cal_extra_params, `[[`, "par_name"),
      min      = sapply(cfg$cal_extra_params, `[[`, "min"),
      max      = sapply(cfg$cal_extra_params, `[[`, "max")
    )
    param_info <- bind_rows(gsa_params, extras)
  } else {
    param_info <- gsa_params
  }
} else {
  # Phases 2-3: inherit param_info from GSA
  param_info <- gsa_results$param_info %>%
    as_tibble() %>%
    select(short, par_name, min, max)
}

cat("Calibration parameters:", nrow(param_info), "\n")
print(param_info[, c("short", "min", "max")])

# ------------------------------------------------------------------------------
# 6. GENERATE PARAMETER SAMPLES (runif)
# ------------------------------------------------------------------------------
n_sim <- cfg$cal$n_simulations

cal_params <- map_dfc(seq_len(nrow(param_info)), function(i) {
  runif(n_sim, min = param_info$min[i], max = param_info$max[i])
})
names(cal_params) <- param_info$par_name
cal_params <- as_tibble(cal_params)

cat("Simulations:", n_sim, "\n")

# ------------------------------------------------------------------------------
# 7. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- lapply(cfg$cal_outputs, function(o) {
  define_output(file = o$file, variable = o$variable, unit = o$unit)
})

# ------------------------------------------------------------------------------
# 8. CLEAN CACHE
# ------------------------------------------------------------------------------
cache_folders <- list.files(cfg$project_path, pattern = "run_|lhoat_|soft_|hard_",
                            full.names = TRUE, include.dirs = TRUE)
cache_folders <- cache_folders[file.info(cache_folders)$isdir]
if (length(cache_folders) > 0) unlink(cache_folders, recursive = TRUE, force = TRUE)

# ------------------------------------------------------------------------------
# 9. RUN SWAT+
# ------------------------------------------------------------------------------
cal_start  <- cfg$cal$start_date
cal_end    <- cfg$cal$end_date
cal_warmup <- cfg$cal$warmup

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = cal_params,
  start_date   = cal_start,
  end_date     = cal_end,
  years_skip   = cal_warmup,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 10. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs:", length(valid_runs), "\n")

# ==============================================================================
# 11. METRICS + FILTERING (depends on calibration type)
# ==============================================================================

if (cfg$cal$type == "soft") {
  # ---------- SOFT CALIBRATION (Phase 1) ----------
  ranking_balance <- calc_soft_metrics(sim, valid_runs, cfg$cal$targets, cfg$cal$threshold_pct)
  satisfactory <- ranking_balance %>% filter(satisfactory_run)

  cat("Satisfactory runs:", nrow(satisfactory), "\n")

  # Monthly streamflow diagnostics for satisfactory runs
  monthly_metrics <- tibble(run_id = character(), NSE_mon = numeric(), KGE_mon = numeric())
  if (nrow(satisfactory) > 0 && "vazao_exutorio" %in% names(sim$simulation)) {
    obs_data <- read.csv(file.path(cfg$project_path, cfg$obs_files$monthly), stringsAsFactors = FALSE) %>%
      mutate(Date = as.Date(Date)) %>%
      filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
      arrange(Date)

    comp_mon <- build_comparativo(
      sim$simulation$vazao_exutorio, obs_data, satisfactory$run_id, floor_month = TRUE
    )

    monthly_metrics <- map_dfr(satisfactory$run_id, function(run) {
      tibble(
        run_id  = run,
        NSE_mon = round(NSE(sim = comp_mon[[run]], obs = comp_mon$Flow), 3),
        KGE_mon = round(KGE(sim = comp_mon[[run]], obs = comp_mon$Flow), 3)
      )
    })
  }

  # Consolidated results
  params_out <- cal_params %>% mutate(run_index = row_number())
  results <- ranking_balance %>%
    left_join(params_out, by = "run_index") %>%
    left_join(monthly_metrics, by = "run_id")

  satisfactory_results <- results %>% filter(satisfactory_run)

  # Plots
  g_scatter <- plot_balance_scatter(ranking_balance, cfg$cal$targets,
                                     paste("Phase", cfg$phase, "- behavioural domain"))
  print(g_scatter)
  save_tiff(g_scatter, paste0("balance_scatter_phase", cfg$phase, ".tif"), ver$path)

  if (nrow(satisfactory) > 0 && "vazao_exutorio" %in% names(sim$simulation)) {
    g_hydro <- plot_hydrograph(comp_mon, satisfactory$run_id,
                                paste("Phase", cfg$phase, "soft calibration - monthly streamflow"))
    print(g_hydro)
    save_tiff(g_hydro, paste0("hydrograph_phase", cfg$phase, ".tif"), ver$path, width = 20, height = 10)
  }

  # Morris ranking for new ranges
  ranking_morris <- gsa_results$ranking %>%
    filter(if ("metric" %in% names(.)) metric == "obj_total" else TRUE)

  # New ranges for next phase
  if (nrow(satisfactory_results) > 0) {
    new_ranges <- calc_new_ranges(
      satisfactory_results,
      param_info$par_name,
      param_info,
      ranking_morris
    )
    print(new_ranges)
    write.csv(new_ranges,
              file.path(ver$path, "novos_ranges_phase2_from_morris_behavioral.csv"),
              row.names = FALSE)
  }

} else {
  # ---------- HARD CALIBRATION (Phases 2-3) ----------
  suffix <- cfg$metric_suffix  # "_mon" or "_day"

  # Primary scale (calibration target)
  primary_key    <- if (cfg$temporal_scale == "daily") "daily" else "monthly"
  primary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_day" else "vazao_exutorio_mon"

  obs_primary <- read.csv(file.path(cfg$project_path, cfg$obs_files[[primary_key]]),
                           stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date))

  if (cfg$temporal_scale == "monthly") {
    obs_primary <- obs_primary %>% mutate(Date = floor_date(Date, "month"))
  }

  obs_primary <- obs_primary %>%
    filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
    arrange(Date)

  floor_mon <- cfg$temporal_scale == "monthly"
  comp_primary <- build_comparativo(sim$simulation[[primary_output]], obs_primary, valid_runs, floor_mon)

  filt <- filter_invalid_runs(comp_primary, valid_runs)
  valid_runs <- filt$runs_ok
  comp_primary <- comp_primary %>% select(Date, Flow, all_of(valid_runs))

  # Primary metrics
  metrics <- calc_hard_metrics(comp_primary, valid_runs, suffix = suffix)
  metrics <- classify_hard_runs(metrics, cfg$cal$threshold_nse, cfg$cal$threshold_kge,
                                 cfg$cal$threshold_abs_pbias, suffix = suffix)

  # Consolidated results
  params_out <- cal_params %>% mutate(run_index = seq_len(n())) %>% relocate(run_index)
  results <- metrics %>% left_join(params_out, by = "run_index")
  satisfactory_results <- results %>% filter(satisfactory_run)

  cat("Satisfactory runs:", nrow(satisfactory_results), "\n")

  # Secondary scale diagnostics
  secondary_key    <- if (cfg$temporal_scale == "daily") "monthly" else "daily"
  secondary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_mon" else "vazao_exutorio_day"

  if (secondary_output %in% names(sim$simulation)) {
    obs_secondary <- read.csv(file.path(cfg$project_path, cfg$obs_files[[secondary_key]]),
                               stringsAsFactors = FALSE) %>%
      mutate(Date = as.Date(Date))

    if (secondary_key == "monthly") {
      obs_secondary <- obs_secondary %>% mutate(Date = floor_date(Date, "month"))
    }

    obs_secondary <- obs_secondary %>%
      filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
      arrange(Date)

    sec_suffix <- if (secondary_key == "monthly") "_mon" else "_day"
    comp_sec <- build_comparativo(sim$simulation[[secondary_output]], obs_secondary,
                                   valid_runs, floor_month = (secondary_key == "monthly"))
    metrics_sec <- calc_hard_metrics(comp_sec, valid_runs, suffix = sec_suffix)
    sec_metric_cols <- grep(paste0("(NSE|KGE|PBIAS)", sec_suffix), names(metrics_sec), value = TRUE)
    results <- results %>% left_join(metrics_sec %>% select(run_id, all_of(sec_metric_cols)), by = "run_id")
  }

  # Water balance diagnostic (Phase 3)
  if (!is.null(cfg$diagnostic_balance) && "precip" %in% names(sim$simulation)) {
    diag_balance <- calc_balance_diagnostic(sim, valid_runs, cfg$diagnostic_balance)
    results <- results %>% left_join(diag_balance, by = "run_id")
  }

  satisfactory_results <- results %>% filter(satisfactory_run)

  # Plots
  nse_col <- paste0("NSE", suffix)
  kge_col <- paste0("KGE", suffix)

  g_scatter <- plot_performance_scatter(metrics, nse_col, kge_col,
                                         cfg$cal$threshold_nse, cfg$cal$threshold_kge,
                                         paste(cfg$temporal_scale, "performance space - Phase", cfg$phase))
  print(g_scatter)
  save_tiff(g_scatter, paste0("scatter_performance_phase", cfg$phase, ".tif"), ver$path)

  g_hydro <- plot_hydrograph(
    comp_primary,
    if (nrow(satisfactory_results) > 0) satisfactory_results$run_id else NULL,
    paste("Phase", cfg$phase, "-", cfg$temporal_scale, "streamflow")
  )
  print(g_hydro)
  save_tiff(g_hydro, paste0("hydrograph_phase", cfg$phase, ".tif"), ver$path, width = 20, height = 10)

  # Parameter boxplot
  if (nrow(satisfactory_results) > 0) {
    g_box <- plot_param_boxplot(satisfactory_results, param_info,
                                 paste("Behavioural parameters - Phase", cfg$phase))
    print(g_box)
    save_tiff(g_box, paste0("boxplot_params_phase", cfg$phase, ".tif"), ver$path, width = 25, height = 10)
  }

  # New ranges for next phase
  ranking_morris <- gsa_results$ranking %>%
    filter(if ("metric" %in% names(.)) metric == "obj_total" else TRUE)

  if (nrow(satisfactory_results) > 0) {
    new_ranges <- calc_new_ranges(
      satisfactory_results,
      param_info$par_name,
      param_info,
      ranking_morris
    )
    print(new_ranges)

    next_phase <- cfg$phase + 1
    ranges_filename <- paste0("novos_ranges_phase", next_phase, "_from_",
                               tolower(cfg$cal$type), "_phase", cfg$phase, ".csv")
    write.csv(new_ranges, file.path(ver$path, ranges_filename), row.names = FALSE)
  }
}

# ==============================================================================
# 12. SAVE RESULTS
# ==============================================================================
phase_label <- cfg$phase

saveRDS(
  list(
    script_version       = cfg$script_version,
    config               = cfg,
    gsa_source           = gsa_results$folder,
    sim                  = sim,
    param_info           = param_info,
    results              = results,
    satisfactory_results = satisfactory_results
  ),
  file.path(ver$path, paste0("resultado_cal_phase", phase_label, ".rds"))
)

write.csv(results,
          file.path(ver$path, paste0("results_cal_phase", phase_label, ".csv")),
          row.names = FALSE)

if (nrow(satisfactory_results) > 0) {
  write.csv(satisfactory_results,
            file.path(ver$path, paste0("satisfactory_results_phase", phase_label, ".csv")),
            row.names = FALSE)
}

cat("\nFiles saved to:", ver$path, "\n")
cat("CAL Phase", phase_label, "COMPLETED.\n")
