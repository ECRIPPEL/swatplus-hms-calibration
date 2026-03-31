# ==============================================================================
# param_inherit.R — Automatic discovery of inherited parameter ranges
# ==============================================================================

library(dplyr)

#' Find the inherited ranges CSV from a previous phase
#'
#' @param base_path SWAT+ project base directory (e.g., "C:/uru")
#' @param from_prefix Folder prefix of the source phase (e.g., "SOFT_phase1_v")
#' @param ranges_file Name of the ranges CSV file
#' @return List with $path (CSV path), $folder (source folder), $data (tibble)
find_inherited_ranges <- function(base_path, from_prefix, ranges_file) {
  folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
  folders <- folders[grepl(paste0("[/\\\\]", from_prefix, "[0-9]+$"), folders)]

  if (length(folders) == 0) {
    stop(paste0("No folder ", from_prefix, "* found in: ", base_path))
  }

  versions <- suppressWarnings(
    as.numeric(sub(paste0(".*", from_prefix), "", basename(folders)))
  )

  order_idx <- order(versions, decreasing = TRUE, na.last = NA)
  folders_sorted <- folders[order_idx]

  for (folder_i in folders_sorted) {
    candidate <- file.path(folder_i, ranges_file)
    if (file.exists(candidate)) {
      data <- read.csv(candidate, stringsAsFactors = FALSE)

      stopifnot(!any(is.na(data$new_min)))
      stopifnot(!any(is.na(data$new_max)))
      stopifnot(all(data$new_min < data$new_max))

      cat("Inherited ranges from: ", folder_i, "\n", sep = "")
      cat("File: ", candidate, "\n", sep = "")
      cat("Parameters: ", nrow(data), "\n", sep = "")

      return(list(path = candidate, folder = folder_i, data = data))
    }
  }

  stop(paste0("File ", ranges_file, " not found in any ", from_prefix, "* folder"))
}

#' Find param_info and ranking from the latest GSA of a given phase
#'
#' @param base_path SWAT+ project base directory
#' @param gsa_prefix GSA folder prefix (e.g., "GSA_phase2_v")
#' @param phase_num Phase number (2 or 3)
#' @return List with $folder, $param_info, $ranking
find_gsa_results <- function(base_path, gsa_prefix, phase_num) {
  folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
  folders <- folders[grepl(paste0("[/\\\\]", gsa_prefix, "[0-9]+$"), folders)]

  if (length(folders) == 0) {
    stop(paste0("No folder ", gsa_prefix, "* found in: ", base_path))
  }

  versions <- suppressWarnings(
    as.numeric(sub(paste0(".*", gsa_prefix), "", basename(folders)))
  )

  order_idx <- order(versions, decreasing = TRUE, na.last = NA)
  folders_sorted <- folders[order_idx]

  pat_param   <- paste0("^param_info_lhoat_phase", phase_num, "(_.*)*\\.csv$")
  pat_ranking <- paste0("^ranking_global_lhoat_phase", phase_num, "(_.*)*\\.csv$")
  pat_metric  <- paste0("^metricas_lhoat_phase", phase_num, "(_.*)*\\.csv$")

  for (folder_i in folders_sorted) {
    cand_param   <- list.files(folder_i, pattern = pat_param, full.names = TRUE)
    cand_ranking <- list.files(folder_i, pattern = pat_ranking, full.names = TRUE)
    cand_metric  <- list.files(folder_i, pattern = pat_metric, full.names = TRUE)

    if (length(cand_param) > 0 && length(cand_ranking) > 0) {
      param_info <- read.csv(cand_param[1], stringsAsFactors = FALSE)
      ranking    <- read.csv(cand_ranking[1], stringsAsFactors = FALSE)

      stopifnot(all(c("short", "par_name", "min", "max") %in% names(param_info)))
      stopifnot(all(param_info$min < param_info$max))

      cat("GSA results from: ", folder_i, "\n", sep = "")

      return(list(
        folder     = folder_i,
        param_info = param_info,
        ranking    = ranking,
        metrics_file = if (length(cand_metric) > 0) cand_metric[1] else NA_character_
      ))
    }
  }

  stop(paste0("No param_info + ranking found in any ", gsa_prefix, "* folder"))
}

#' Build param_info from inherited ranges + parameter map + optional extras
#'
#' @param ranges_data Tibble of inherited ranges (columns: parameter, new_min, new_max)
#' @param param_map Named vector (short -> SWATrunR par_name)
#' @param extra_params Optional tibble of extra parameters
#' @return Tibble param_info ready for use
build_param_info <- function(ranges_data, param_map, extra_params = NULL) {
  param_info <- ranges_data %>%
    mutate(
      par_name = unname(param_map[parameter]),
      short    = parameter,
      min      = new_min,
      max      = new_max
    ) %>%
    filter(!is.na(par_name)) %>%
    select(short, par_name, min, max,
           any_of(c("grupo", "q_low", "q_high", "post_min", "post_max", "mediana")))

  if (!is.null(extra_params) && nrow(extra_params) > 0) {
    extra_params <- extra_params %>%
      filter(!short %in% param_info$short)
    param_info <- bind_rows(param_info, extra_params)
  }

  stopifnot(length(unique(param_info$short)) == nrow(param_info))
  stopifnot(all(param_info$min < param_info$max))

  param_info
}
