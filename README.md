# MATH-5472-Course-Project
Reproducible code for the paper
“What If Without the Proportional Hazards Method?”
Author: Baoyi LIN (21206217)
Course: MATH 5472 — Computer Age Statistical Inference

Overview

This repository contains an end-to-end reproducible pipeline for:

Simulation study comparing:

Cox Proportional Hazards model

Accelerated Failure Time (AFT) model

Time-varying coefficient Cox model (non-PH)

Real-data analysis using the PBC (Primary Biliary Cirrhosis) dataset from the survival R package.

Both experiments support the conclusions in the paper regarding the strengths and limitations of the Proportional Hazards (PH) assumption.

Repository Structure
cox-nonph-analysis/
│
├── R/
│   ├── 01_simulation_ph_vs_nonph.R        # Simulation study
│   ├── 02_realdata_pbc_analysis.R         # Real PBC dataset analysis
│
├── results/
│   ├── simulation/                        # Output from simulations
│   ├── pbc/                               # Output from real data analysis
│
├── figures/                               # Plots used in the paper
│
├── README.md                              # Project documentation
└── .gitignore

Installation

The code requires R ≥ 4.2.0.

Install required packages:

install.packages(c("survival", "dplyr", "purrr", "ggplot2", "survminer", "pec"))

How to Reproduce the Experiments
▶ 1. Simulation Study

Run the following script in R:

source("R/01_simulation_ph_vs_nonph.R")


This script will:

Generate datasets under PH and non-PH scenarios

Fit Cox / AFT / time-varying Cox models

Compute bias, RMSE, coverage, C-index, IBS

Save results to: results/simulation/

Produce figures under: figures/simulation/

▶ 2. Real Data Analysis (PBC dataset)

Run:

source("R/02_realdata_pbc_analysis.R")


This script will:

Load the PBC dataset

Fit Cox model and test PH assumption

Plot Schoenfeld residuals and survival curves

Compare Cox vs AFT vs time-varying Cox

Save outputs to: results/pbc/

Generate figures used in the paper

Reproducibility

This repository follows the NeurIPS reproducibility guidelines:

All analyses are contained in R scripts

All results can be reproduced by running the corresponding script once dependencies are installed

Figures used in the paper are stored in figures/ and regenerated automatically

No external private datasets are required

Contact

For questions, please contact:
Baoyi LIN — blinan@connect.ust.hk