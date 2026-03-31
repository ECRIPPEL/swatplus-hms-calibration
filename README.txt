satplus-hms-calibration

HIERARCHICAL MULTI-SCALE CALIBRATION OF SWAT+ 
GUIDED BY ADAPTIVE SENSITIVITY ANALYSIS


How to use:

# Run a specific phase:
Rscript run_gsa.R config/phase1.yaml
Rscript run_cal.R config/phase1.yaml

# Or run everything sequentially:
source("master.R")

# Resume from a specific step (e.g., GSA2 = step 3):
START_FROM <- 3
source("master.R")