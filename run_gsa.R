# ==============================================================================
# run_gsa.R — Main GSA script (generic for any phase)
#
# USAGE:
#   Rscript run_gsa.R config/phase1.yaml       # from terminal
#
#   PHASE_CONFIG <- "config/phase1.yaml"        # from RStudio/VSCode
#   source("run_gsa.R")
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
library(lhs)
library(remotes)
library(purrr)


# ------------------------------------------------------------------------------
# 1. READ CONFIG
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (exists("PHASE_CONFIG")) {
  config_path <- PHASE_CONFIG
} else if (length(args) >= 1) {
  config_path <- args[1]
} else {
  stop("Usage: Rscript run_gsa.R config/phaseX.yaml\n  Or set PHASE_CONFIG before source().")
}

cfg <- yaml::read_yaml(config_path)
cat("=== GSA Phase", cfg$phase, "===\n")
cat("Config:", config_path, "\n")

# ------------------------------------------------------------------------------
# 2. LOAD MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/lhoat_engine.R")
source("R/param_inherit.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")
if (isTRUE(cfg$check_range$enabled)) source("R/check_range.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. VERSIONING
# ------------------------------------------------------------------------------
ver <- create_versioned_folder(cfg$project_path, cfg$gsa$prefix)
cat("TAG:", ver$tag, "\nFOLDER:", ver$path, "\n")

# ------------------------------------------------------------------------------
# 4. PARAMETERS
# ------------------------------------------------------------------------------
if (cfg$phase == 1) {
  # Phase 1: parameters defined in YAML
  param_info <- tibble(
    short    = sapply(cfg$parameters, `[[`, "short"),
    par_name = sapply(cfg$parameters, `[[`, "par_name"),
    min      = sapply(cfg$parameters, `[[`, "min"),
    max      = sapply(cfg$parameters, `[[`, "max")
  )
} else {
  # Phases 2-3: inherit ranges from previous phase
  inherited <- find_inherited_ranges(
    cfg$project_path,
    cfg$inherit$from_prefix,
    cfg$inherit$ranges_file
  )

  param_map <- unlist(cfg$param_map)

  extra_params <- NULL
  if (!is.null(cfg$extra_params)) {
    extra_params <- tibble(
      short    = sapply(cfg$extra_params, `[[`, "short"),
      par_name = sapply(cfg$extra_params, `[[`, "par_name"),
      min      = sapply(cfg$extra_params, `[[`, "min"),
      max      = sapply(cfg$extra_params, `[[`, "max"),
      grupo    = NA_character_,
      q_low    = NA_real_,
      q_high   = NA_real_,
      post_min = NA_real_,
      post_max = NA_real_,
      mediana  = NA_real_
    )
  }

  param_info <- build_param_info(inherited$data, param_map, extra_params)
}

p <- nrow(param_info)
cat("Parameters:", p, "\n")
print(param_info)

# ------------------------------------------------------------------------------
# 5. CHECK PARAMETER RANGE (Phase 1 only)
# ------------------------------------------------------------------------------
if (isTRUE(cfg$check_range$enabled)) {
  check_result <- check_parameter_range(param_info, cfg$project_path, cfg$check_range$map)
}

# ------------------------------------------------------------------------------
# 6. GENERATE LH-OAT TRAJECTORIES
# ------------------------------------------------------------------------------
gsa_start  <- if (!is.null(cfg$gsa$start_date)) cfg$gsa$start_date else cfg$start_date
gsa_end    <- if (!is.null(cfg$gsa$end_date))   cfg$gsa$end_date   else cfg$end_date
gsa_warmup <- if (!is.null(cfg$gsa$warmup))     cfg$gsa$warmup     else cfg$warmup

design_norm <- generate_lhoat_trajectories(
  m_traj      = cfg$gsa$m_traj,
  p           = p,
  delta       = cfg$gsa$delta,
  param_names = param_info$short
)

scaled <- scale_parameters(design_norm, param_info)

cat("Expected runs:", nrow(scaled$swat_params), "\n")

# ------------------------------------------------------------------------------
# 7. CLEAN CACHE
# ------------------------------------------------------------------------------
cache_folders <- list.files(cfg$project_path, pattern = "run_|lhoat_",
                            full.names = TRUE, include.dirs = TRUE)
cache_folders <- cache_folders[file.info(cache_folders)$isdir]
if (length(cache_folders) > 0) unlink(cache_folders, recursive = TRUE, force = TRUE)

# ------------------------------------------------------------------------------
# 8. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- lapply(cfg$gsa_outputs, function(o) {
  define_output(file = o$file, variable = o$variable, unit = o$unit)
})

# ------------------------------------------------------------------------------
# 9. RUN SWAT+
# ------------------------------------------------------------------------------
sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = scaled$swat_params,
  start_date   = gsa_start,
  end_date     = gsa_end,
  years_skip   = gsa_warmup,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 10. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs:", length(valid_runs), "\n")

# ------------------------------------------------------------------------------
# 11. METRICS
# ------------------------------------------------------------------------------
if (cfg$phase == 1) {
  # Phase 1: annual water balance ratios
  metrics <- calc_gsa_phase1_metrics(sim, valid_runs, cfg$gsa_targets)
} else {
  # Phases 2-3: streamflow metrics (NSE/KGE/PBIAS)
  obs_file_key <- if (cfg$temporal_scale == "daily") "daily" else "monthly"
  obs_path <- file.path(cfg$project_path, cfg$obs_files[[obs_file_key]])
  obs_data <- read.csv(obs_path, stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date))

  if (cfg$temporal_scale == "monthly") {
    obs_data <- obs_data %>%
      mutate(Date = lubridate::floor_date(Date, "month"))
  }

  obs_data <- obs_data %>%
    filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(gsa_end)) %>%
    arrange(Date)

  floor_mon <- cfg$temporal_scale == "monthly"
  comparativo <- build_comparativo(sim$simulation$vazao_exutorio, obs_data, valid_runs, floor_mon)

  filt <- filter_invalid_runs(comparativo, valid_runs)
  valid_runs <- filt$runs_ok
  comparativo <- comparativo %>% select(Date, Flow, all_of(valid_runs))

  metrics <- calc_hard_metrics(comparativo, valid_runs, suffix = cfg$metric_suffix)
  metrics <- metrics %>% mutate(run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_id)))
}

# ------------------------------------------------------------------------------
# 12. ATTACH TRAJECTORY + COMPUTE EE
# ------------------------------------------------------------------------------
metrics <- design_norm %>%
  select(run_index, traj_id, step_id, all_of(param_info$short)) %>%
  inner_join(metrics, by = "run_index") %>%
  arrange(traj_id, step_id)

ee_all <- calc_all_ee(metrics, cfg$gsa_metrics, param_info$short)
sens   <- compute_sensitivity(ee_all)

cat("\n=== GLOBAL RANKING ===\n")
print(sens$ranking)

# ------------------------------------------------------------------------------
# 13. MORRIS SCREENING PLOT
# ------------------------------------------------------------------------------
ranking_classified <- classify_morris(sens$ranking)
g_morris <- plot_morris(ranking_classified, paste("Phase", cfg$phase))
print(g_morris)
save_tiff(g_morris, paste0("Morris_sensitivity_screening_phase", cfg$phase, ".tif"), ver$path)

# ------------------------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------------------------
phase_label <- cfg$phase

saveRDS(
  list(
    script_version    = cfg$script_version,
    config            = cfg,
    sim               = sim,
    ranking           = sens$ranking,
    sensitivity       = sens$sensitivity,
    metrics           = metrics,
    ee_all            = ee_all,
    param_info        = param_info,
    lhoat_design_norm = design_norm,
    swat_params       = scaled$swat_params
  ),
  file.path(ver$path, paste0("resultado_lhoat_phase", phase_label, ".rds"))
)

write.csv(sens$ranking, file.path(ver$path, paste0("ranking_global_lhoat_phase", phase_label, ".csv")), row.names = FALSE)
write.csv(metrics,      file.path(ver$path, paste0("metricas_lhoat_phase", phase_label, ".csv")), row.names = FALSE)
write.csv(ee_all,       file.path(ver$path, paste0("ee_all_lhoat_phase", phase_label, ".csv")), row.names = FALSE)
write.csv(param_info,   file.path(ver$path, paste0("param_info_lhoat_phase", phase_label, ".csv")), row.names = FALSE)

cat("\nFiles saved to:", ver$path, "\n")
cat("GSA Phase", phase_label, "COMPLETED.\n")
