# What If Without the Proportional Hazards Method?

### MATH 5472 Course Project · Survival Analysis under PH and Non-PH Scenarios

Author: **LIN Baoyi (21206217)**\
Reproducibility Website: <https://github.com/BaoyiLINN/MATH-5472-Course-Project>

------------------------------------------------------------------------

## Overview

This repository contains all reproducible code for the course project **"What If Without the Proportional Hazards Method?"**, including:

1.  **Simulation experiments** comparing Cox PH, AFT, and time-varying coefficient models under proportional and non-proportional hazards.
2.  **Real-world analysis** on the Mayo Clinic **Primary Biliary Cirrhosis (PBC)** dataset:
    -   PH diagnostics (Schoenfeld residuals)
    -   Kaplan–Meier stratification
    -   Prediction error curves (PEC) and Integrated Brier Score (IBS)

All figures used in the final report are generated from the scripts in this repository.

------------------------------------------------------------------------

## Repository Structure

```         
.
├── 01_simulation_ph_vs_nonph.R      # Simulation study: PH vs Non-PH scenarios
├── 02_realdata_pbc_analysis.R       # Real-data analysis on PBC dataset
├── figures/                         # Auto-generated figures (after running scripts)
│   ├── simulation/
│   └── pbc/
└── results/                         # RDS / CSV intermediate outputs
```

------------------------------------------------------------------------

# 1. Simulation Study

### File: `01_simulation_ph_vs_nonph.R`

This script reproduces **Section 4.1** of the paper:

-   Scenario A: True PH\
-   Scenario B: Strongly violated PH (effect reversal)\
-   Models compared:
    -   Cox Proportional Hazards
    -   Log-normal AFT
    -   Time-varying coefficient Cox model\
-   Metrics computed:
    -   Bias and RMSE of effect estimates
    -   95% CI coverage
    -   C-index
    -   Integrated Brier Score (IBS)

### To run:

``` bash
Rscript 01_simulation_ph_vs_nonph.R
```

Outputs will be saved into:

```         
figures/simulation/
results/simulation/
```

------------------------------------------------------------------------

# 2. Real-World Data Experiment (PBC)

### File: `02_realdata_pbc_analysis.R`

This script performs:

-   Cox PH model fitting
-   Schoenfeld residual diagnostics\
-   KM curves (high vs low bilirubin)
-   Time-varying coefficient Cox modeling
-   PEC curves and IBS (Cox vs TVC)

### To run:

``` bash
Rscript 02_realdata_pbc_analysis.R
```

Figures will be generated in:

```         
figures/pbc/
```

------------------------------------------------------------------------

## Installation & Dependencies

The analysis uses standard R survival-analysis packages.

### Install all required packages:

``` r
install.packages(c(
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "splines",
  "pec"
))
```

------------------------------------------------------------------------

## How to Fully Reproduce the Paper

1.  Clone this repository:

``` bash
git clone https://github.com/BaoyiLINN/MATH-5472-Course-Project.git
cd MATH-5472-Course-Project
```

2.  Run simulation:

``` bash
Rscript 01_simulation_ph_vs_nonph.R
```

3.  Run PBC analysis:

``` bash
Rscript 02_realdata_pbc_analysis.R
```

4.  All figures will appear in the `figures/` directory.

5.  Insert figures into the LaTeX paper (neurips format).

------------------------------------------------------------------------

## Session Information

To ensure reproducibility, your R session info is automatically saved into:

```         
results/sessionInfo.txt
```

------------------------------------------------------------------------
