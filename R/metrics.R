# ==============================================================================
# metrics.R — Performance metrics and water balance computation
# ==============================================================================

library(dplyr)
library(tibble)
library(purrr)
library(hydroGOF)

#' Compute NSE, KGE, PBIAS, and obj_total for each run
#'
#' @param comparativo data.frame with Date, Flow, and run_* columns
#' @param valid_runs Character vector of run names
#' @param suffix Metric column suffix (e.g., "_mon", "_day")
#' @return Tibble with metrics per run
calc_hard_metrics <- function(comparativo, valid_runs, suffix = "_mon") {
  map_dfr(valid_runs, function(run_nm) {
    sim_vals <- comparativo[[run_nm]]
    obs_vals <- comparativo$Flow

    if (all(is.na(sim_vals)) || length(na.omit(sim_vals)) < 2 ||
        sd(sim_vals, na.rm = TRUE) == 0) {
      out <- tibble(
        run_id    = run_nm,
        NSE       = NA_real_,
        KGE       = NA_real_,
        PBIAS     = NA_real_,
        absPBIAS  = NA_real_,
        err_nse   = NA_real_,
        err_kge   = NA_real_,
        err_pbias = NA_real_,
        obj_total = NA_real_,
        run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
      )
      names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
      return(out)
    }

    nse_val   <- NSE(sim = sim_vals, obs = obs_vals)
    kge_val   <- KGE(sim = sim_vals, obs = obs_vals)
    pbias_val <- pbias(sim = sim_vals, obs = obs_vals)

    nse_n   <- as.numeric(nse_val)
    kge_n   <- as.numeric(kge_val)
    pbias_n <- as.numeric(pbias_val)

    out <- tibble(
      run_id    = run_nm,
      NSE       = nse_n,
      KGE       = kge_n,
      PBIAS     = pbias_n,
      absPBIAS  = abs(pbias_n),
      err_nse   = pmax(0, 1 - pmin(nse_n, 1)),
      err_kge   = pmax(0, 1 - pmin(kge_n, 1)),
      err_pbias = abs(pbias_n) / 100,
      obj_total = pmax(0, 1 - pmin(nse_n, 1)) + pmax(0, 1 - pmin(kge_n, 1)) + abs(pbias_n) / 100,
      run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
    )

    names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
    out
  })
}

#' Classify runs as satisfactory (hard calibration thresholds)
#'
#' @param metrics Tibble of metrics (output of calc_hard_metrics)
#' @param threshold_nse Minimum NSE threshold
#' @param threshold_kge Minimum KGE threshold
#' @param threshold_abs_pbias Maximum |PBIAS| threshold
#' @param suffix Column suffix (e.g., "_mon", "_day")
#' @return Tibble with added column satisfactory_run
classify_hard_runs <- function(metrics, threshold_nse, threshold_kge,
                               threshold_abs_pbias, suffix = "_mon") {
  nse_col <- paste0("NSE", suffix)
  kge_col <- paste0("KGE", suffix)

  metrics %>%
    mutate(
      satisfactory_run =
        .data[[nse_col]] >= threshold_nse &
        .data[[kge_col]] >= threshold_kge &
        absPBIAS <= threshold_abs_pbias
    ) %>%
    arrange(desc(satisfactory_run), desc(.data[[nse_col]]), desc(.data[[kge_col]]), absPBIAS)
}

#' Compute annual water balance metrics (soft calibration — Phase 1)
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param targets List with et_rto, wyld_rto
#' @param threshold_pct Maximum relative error per component (%)
#' @return Tibble with ranking and satisfactory_run flag
calc_soft_metrics <- function(sim, valid_runs, targets, threshold_pct) {
  ranking <- tibble(
    run_id   = valid_runs,
    P_sim    = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip), valid_runs]),
    ET_sim   = as.numeric(sim$simulation$et[nrow(sim$simulation$et), valid_runs]),
    WYLD_sim = as.numeric(sim$simulation$wateryld[nrow(sim$simulation$wateryld), valid_runs])
  ) %>%
    mutate(
      rto_et   = ET_sim / P_sim,
      rto_wyld = WYLD_sim / P_sim,
      target_et   = targets$et_rto,
      target_wyld = targets$wyld_rto,
      err_et   = abs(rto_et - target_et),
      err_wyld = abs(rto_wyld - target_wyld),
      err_rel_et_pct   = 100 * err_et / target_et,
      err_rel_wyld_pct = 100 * err_wyld / target_wyld,
      agg_rel_error_pct = err_rel_et_pct + err_rel_wyld_pct,
      satisfactory_run =
        err_rel_et_pct   <= threshold_pct &
        err_rel_wyld_pct <= threshold_pct,
      run_index = as.numeric(sub("^run_0*([0-9]+).*", "\\1", run_id))
    ) %>%
    arrange(desc(satisfactory_run), agg_rel_error_pct)

  ranking
}

#' Compute GSA Phase 1 metrics (annual water balance ratios)
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param targets List with q_p, et_p, bfi
#' @return Tibble with metrics per run
calc_gsa_phase1_metrics <- function(sim, valid_runs, targets) {
  tibble(
    run_id = valid_runs,
    precip = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip), valid_runs]),
    et     = as.numeric(sim$simulation$et[nrow(sim$simulation$et), valid_runs]),
    surq   = as.numeric(sim$simulation$surq[nrow(sim$simulation$surq), valid_runs]),
    latq   = as.numeric(sim$simulation$latq[nrow(sim$simulation$latq), valid_runs]),
    gw_q   = as.numeric(sim$simulation$gw_q[nrow(sim$simulation$gw_q), valid_runs])
  ) %>%
    mutate(
      q_total   = surq + latq + gw_q,
      q_p       = q_total / precip,
      et_p      = et / precip,
      bfi       = gw_q / q_total,
      err_q_p   = abs(q_p - targets$q_p),
      err_et_p  = abs(et_p - targets$et_p),
      err_bfi   = abs(bfi - targets$bfi),
      obj_total = err_q_p + err_et_p + err_bfi,
      run_index = as.numeric(sub("^run_0*([0-9]+).*", "\\1", run_id))
    )
}

#' Compute water balance diagnostic (Phase 3)
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param targets List with et_rto, wyld_rto
#' @return Tibble with ET/P, WYLD/P, and error metrics
calc_balance_diagnostic <- function(sim, valid_runs, targets) {
  tibble(
    run_id   = valid_runs,
    P_sim    = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip), valid_runs]),
    ET_sim   = as.numeric(sim$simulation$et[nrow(sim$simulation$et), valid_runs]),
    WYLD_sim = as.numeric(sim$simulation$wateryld[nrow(sim$simulation$wateryld), valid_runs])
  ) %>%
    mutate(
      ET_P   = ET_sim / P_sim,
      WYLD_P = WYLD_sim / P_sim,
      target_ET_P   = targets$et_rto,
      target_WYLD_P = targets$wyld_rto,
      err_ET_P_abs   = abs(ET_P - target_ET_P),
      err_WYLD_P_abs = abs(WYLD_P - target_WYLD_P),
      err_ET_P_pct   = 100 * err_ET_P_abs / target_ET_P,
      err_WYLD_P_pct = 100 * err_WYLD_P_abs / target_WYLD_P,
      agg_rel_error_pct = err_ET_P_pct + err_WYLD_P_pct
    )
}
