# swatplus-hms-calibration

A hierarchical multi-scale SWAT+ calibration framework using LH-OAT Morris 
screening and progressive parameter narrowing across three phases.


## Workflow

GSA Phase 1 → CAL Phase 1 → GSA Phase 2 → CAL Phase 2 → GSA Phase 3 → CAL Phase 3


- **Phase 1** — Annual water balance (soft calibration via behavioural filtering)
- **Phase 2** — Monthly streamflow (hard calibration using NSE, KGE, and PBIAS)
- **Phase 3** — Daily streamflow (hard calibration using NSE, KGE, and PBIAS)

## Quick Start

```r
setwd("C:/path/to/your/swatplus-hms-calibration")
source("master.R")

Requirements
R = 4.3.3, SWAT+ >= swatplus-61.0.2.61-ifx-win_amd64-Rel, RTools = 4.3
R packages: SWATrunR, yaml, dplyr, hydroGOF, ggplot2, lhs, and others

This framework relies on [SWATrunR](https://github.com/chrisschuerz/SWATrunR) 
by Christoph Schürz for running SWAT+ simulations in R.

Documentation
See HOW_TO_USE.md for full setup, configuration, and troubleshooting.
