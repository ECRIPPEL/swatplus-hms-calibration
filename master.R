# ==============================================================================
# master.R — Orchestrator: GSA1 -> CAL1 -> GSA2 -> CAL2 -> GSA3 -> CAL3
#
# Runs all 6 steps sequentially.
# Can be interrupted and resumed from any step.
#
# USAGE:
#   source("master.R")                    # run all steps
#   START_FROM <- 3; source("master.R")   # resume from step 3 (GSA2)
# ==============================================================================

# Starting step (default = 1)
if (!exists("START_FROM")) START_FROM <- 1

steps <- list(
  list(step = 1, label = "GSA Phase 1",  script = "run_gsa.R", config = "config/phase1.yaml"),
  list(step = 2, label = "CAL Phase 1",  script = "run_cal.R", config = "config/phase1.yaml"),
  list(step = 3, label = "GSA Phase 2",  script = "run_gsa.R", config = "config/phase2.yaml"),
  list(step = 4, label = "CAL Phase 2",  script = "run_cal.R", config = "config/phase2.yaml"),
  list(step = 5, label = "GSA Phase 3",  script = "run_gsa.R", config = "config/phase3.yaml"),
  list(step = 6, label = "CAL Phase 3",  script = "run_cal.R", config = "config/phase3.yaml")
)

cat("==============================================================\n")
cat(" SWAT+ HIERARCHICAL MULTI-SCALE CALIBRATION\n")
cat(" GSA1 -> CAL1 -> GSA2 -> CAL2 -> GSA3 -> CAL3\n")
cat("==============================================================\n")
cat(" Starting from step:", START_FROM, "\n\n")

for (s in steps) {
  if (s$step < START_FROM) {
    cat("[SKIP] Step", s$step, "-", s$label, "\n")
    next
  }

  cat("\n==============================================================\n")
  cat(" STEP", s$step, "/6 -", s$label, "\n")
  cat("==============================================================\n")

  t_start <- Sys.time()

  PHASE_CONFIG <- s$config
  .MASTER_STATE <- list(steps = steps, s = s, START_FROM = START_FROM, t_start = t_start)
  source(s$script, local = FALSE)
  list2env(.MASTER_STATE, envir = globalenv())
  rm(.MASTER_STATE)

  t_end <- Sys.time()
  elapsed <- round(difftime(t_end, t_start, units = "mins"), 1)

  cat("\n[DONE] Step", s$step, "-", s$label, "in", elapsed, "minutes\n")
}

cat("\n==============================================================\n")
cat(" CALIBRATION COMPLETE!\n")
cat("==============================================================\n")
