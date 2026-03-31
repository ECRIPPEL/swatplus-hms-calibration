# ==============================================================================
# plots.R — Publication-ready scientific plots
# Morris screening, hydrographs, scatter, boxplot
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# Register font on Windows
if (.Platform$OS.type == "windows") {
  windowsFonts(TimesNR = windowsFont("Times New Roman"))
}

BASE_FAMILY <- "TimesNR"

# ==============================================================================
# MORRIS SCREENING PLOT
# ==============================================================================

#' Morris screening plot (mu* vs sigma)
#'
#' @param ranking Classified tibble (with grupo column, from classify_morris)
#' @param phase_label Phase label (e.g., "Phase 1")
#' @return ggplot object
plot_morris <- function(ranking, phase_label = "Phase") {
  mu_med    <- mean(ranking$mu_star, na.rm = TRUE)
  sigma_med <- mean(ranking$sigma,   na.rm = TRUE)

  ranking <- ranking %>%
    mutate(
      plot_label = if_else(mu_star > 0 | sigma > 0, parameter, "")
    )

  ggplot(ranking, aes(x = mu_star, y = sigma)) +
    geom_hline(yintercept = sigma_med, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = mu_med, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(
      aes(fill = grupo, shape = grupo),
      color = "black", size = 3.8, stroke = 0.5, alpha = 0.8
    ) +
    geom_text_repel(
      aes(label = plot_label),
      family = BASE_FAMILY, size = 4, color = "black",
      box.padding = 0.5, point.padding = 0.3,
      segment.color = "grey70", segment.size = 0.4,
      min.segment.length = 0.2, max.overlaps = Inf
    ) +
    scale_fill_manual(values = c(
      "Low importance" = "grey85",
      "High importance" = "grey55",
      "High importance + interaction" = "black"
    )) +
    scale_shape_manual(values = c(
      "Low importance" = 21,
      "High importance" = 22,
      "High importance + interaction" = 24
    )) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = paste("Morris sensitivity screening -", phase_label),
      x = expression(mu^"*"),
      y = expression(sigma),
      fill = NULL, shape = NULL
    ) +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      plot.title  = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      axis.title  = element_text(size = 14, color = "black"),
      axis.text   = element_text(size = 12, color = "black"),
      legend.position = "top",
      legend.text = element_text(size = 11),
      axis.line   = element_line(color = "black", linewidth = 0.5),
      axis.ticks  = element_line(color = "black")
    )
}

# ==============================================================================
# HYDROGRAPH — observed vs simulated (all runs + optional behavioural highlight)
# ==============================================================================

#' Hydrograph plot: observed vs simulated streamflow
#'
#' @param comparativo data.frame with Date, Flow, and run_* columns
#' @param behavioral_runs Character vector of satisfactory run IDs (NULL = no highlight)
#' @param title_label Plot title
#' @return ggplot object
plot_hydrograph <- function(comparativo, behavioral_runs = NULL,
                            title_label = "Observed vs simulated streamflow") {
  all_runs <- setdiff(names(comparativo), c("Date", "Flow"))

  sim_long <- comparativo %>%
    select(-Flow) %>%
    pivot_longer(cols = -Date, names_to = "run_id", values_to = "sim")

  if (!is.null(behavioral_runs)) {
    sim_long <- sim_long %>%
      mutate(
        group = ifelse(run_id %in% behavioral_runs, "Behavioral", "Non-behavioral")
      )

    g <- ggplot() +
      geom_line(
        data = sim_long %>% filter(group == "Non-behavioral"),
        aes(x = Date, y = sim, group = run_id),
        colour = "grey85", alpha = 0.60, linewidth = 0.22
      ) +
      geom_line(
        data = sim_long %>% filter(group == "Behavioral"),
        aes(x = Date, y = sim, group = run_id),
        colour = "#5DADE2", alpha = 0.75, linewidth = 0.30
      )
  } else {
    g <- ggplot() +
      geom_line(
        data = sim_long,
        aes(x = Date, y = sim, group = run_id),
        color = "lightblue", linewidth = 0.4, alpha = 0.7
      )
  }

  g +
    geom_line(
      data = comparativo,
      aes(x = Date, y = Flow),
      colour = "black", linewidth = 0.45
    ) +
    labs(title = title_label, x = NULL, y = expression(Streamflow ~ (m^3/s))) +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title   = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 11, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      plot.margin  = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# SCATTER — NSE vs KGE performance space
# ==============================================================================

#' Performance space scatter plot (NSE vs KGE)
#'
#' @param metrics Tibble with NSE, KGE, satisfactory_run
#' @param nse_col Name of the NSE column (e.g., "NSE_mon")
#' @param kge_col Name of the KGE column (e.g., "KGE_mon")
#' @param threshold_nse NSE reference line
#' @param threshold_kge KGE reference line
#' @param title_label Plot title
#' @return ggplot object
plot_performance_scatter <- function(metrics, nse_col, kge_col,
                                     threshold_nse, threshold_kge,
                                     title_label = "Performance space") {
  ggplot() +
    geom_point(
      data = metrics %>% filter(!satisfactory_run),
      aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape = 21, size = 2.4, stroke = 0.20, color = "black", fill = "grey85", alpha = 0.75
    ) +
    geom_point(
      data = metrics %>% filter(satisfactory_run),
      aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape = 21, size = 2.8, stroke = 0.22, color = "black", fill = "#6BAED6", alpha = 0.90
    ) +
    geom_vline(xintercept = threshold_nse, linetype = "dashed", color = "grey35") +
    geom_hline(yintercept = threshold_kge, linetype = "dashed", color = "grey35") +
    labs(title = title_label, x = "NSE", y = "KGE") +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text  = element_text(size = 11, colour = "black"),
      legend.position = "none"
    )
}

# ==============================================================================
# SCATTER — ET/P vs WYLD/P (water balance)
# ==============================================================================

#' Water balance scatter plot (ET/P vs WYLD/P)
#'
#' @param ranking Tibble with rto_et, rto_wyld, satisfactory_run
#' @param targets List with et_rto, wyld_rto
#' @param title_label Plot title
#' @return ggplot object
plot_balance_scatter <- function(ranking, targets, title_label = "Water balance domain") {
  target_df <- tibble(rto_et = targets$et_rto, rto_wyld = targets$wyld_rto)

  ggplot() +
    geom_point(
      data = ranking %>% filter(!satisfactory_run),
      aes(x = rto_et, y = rto_wyld),
      shape = 21, size = 2.0, stroke = 0.18, colour = "grey35", fill = "grey88", alpha = 0.60
    ) +
    geom_point(
      data = ranking %>% filter(satisfactory_run),
      aes(x = rto_et, y = rto_wyld),
      shape = 21, size = 2.5, stroke = 0.22, colour = "black", fill = "grey35", alpha = 0.92
    ) +
    stat_ellipse(
      data = ranking %>% filter(satisfactory_run),
      aes(x = rto_et, y = rto_wyld),
      type = "norm", level = 0.95, linewidth = 0.60, colour = "grey10"
    ) +
    geom_vline(xintercept = targets$et_rto, linetype = "dashed", linewidth = 0.40, colour = "grey35") +
    geom_hline(yintercept = targets$wyld_rto, linetype = "dashed", linewidth = 0.40, colour = "grey35") +
    geom_point(data = target_df, aes(x = rto_et, y = rto_wyld),
               shape = 4, size = 3.4, stroke = 0.95, colour = "red") +
    labs(title = title_label, x = "ET / P", y = "WYLD / P") +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      plot.margin = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# BOXPLOT — Normalized behavioural parameters
# ==============================================================================

#' Normalized boxplot of behavioural parameter distributions
#'
#' @param satisfactory_results Tibble of satisfactory runs
#' @param param_info Tibble with short, par_name, min, max
#' @param title_label Plot title
#' @return ggplot object
plot_param_boxplot <- function(satisfactory_results, param_info,
                               title_label = "Behavioural parameter distributions") {
  prior_df <- param_info %>%
    transmute(parameter = short, par_name, prior_min = min, prior_max = max)

  param_long <- satisfactory_results %>%
    select(all_of(param_info$par_name)) %>%
    pivot_longer(cols = everything(), names_to = "par_name", values_to = "value") %>%
    left_join(param_info %>% select(short, par_name), by = "par_name") %>%
    rename(parameter = short) %>%
    filter(!is.na(parameter), !is.na(value))

  param_norm <- param_long %>%
    left_join(prior_df, by = c("parameter", "par_name")) %>%
    mutate(
      value_norm = case_when(
        is.na(prior_min) | is.na(prior_max) ~ NA_real_,
        prior_max <= prior_min ~ NA_real_,
        TRUE ~ (value - prior_min) / (prior_max - prior_min)
      )
    ) %>%
    filter(!is.na(value_norm))

  param_norm$parameter <- factor(param_norm$parameter, levels = param_info$short)

  ggplot(param_norm, aes(x = parameter, y = value_norm)) +
    geom_jitter(width = 0.22, alpha = 0.28, color = "black", size = 1.1) +
    geom_boxplot(
      width = 0.28, fill = "white", color = "black",
      linewidth = 0.40, outlier.shape = NA, alpha = 0.85
    ) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "grey45", linewidth = 0.45) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey65", linewidth = 0.40) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title_label, x = NULL, y = "Normalized parameter value") +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, face = "italic", size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.line    = element_line(color = "black", linewidth = 0.5),
      axis.ticks   = element_line(color = "black")
    )
}

# ==============================================================================
# HELPER — Save plot as TIFF (publication-ready)
# ==============================================================================

#' Save plot as TIFF for publication
#'
#' @param plot ggplot object
#' @param filename File name (without path)
#' @param save_path Destination folder
#' @param width Width in cm (default 18)
#' @param height Height in cm (default 15)
save_tiff <- function(plot, filename, save_path, width = 18, height = 15) {
  ggsave(
    filename    = file.path(save_path, filename),
    plot        = plot,
    width       = width,
    height      = height,
    units       = "cm",
    dpi         = 300,
    compression = "lzw",
    device      = "tiff"
  )
}
