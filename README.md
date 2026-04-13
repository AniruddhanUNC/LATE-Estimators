# Benchmarking LATE Estimators

This repository contains the code and report for a graduate capstone project comparing five estimators for the Local Average Treatment Effect (LATE).

## Methods compared

- Two-Stage Least Squares (2SLS)
- Post-Double-Selection Lasso IV (Post-Lasso IV)
- Double Machine Learning IV (DML-IV)
- Augmented Inverse Probability Weighted IV (AIPW-IV)
- IV Targeted Maximum Likelihood Estimation (IV-TMLE)

All five methods target the LATE under the standard IV assumptions of Imbens and Angrist (1994). The comparison covers bias, RMSE, confidence interval coverage, and CI width across three data-generating processes (linear, nonlinear, and sparse nuisance) and a grid of sample sizes and instrument strengths.

## Repository structure

```
.
├── late_simulation.R              # Simulation code (all five estimators)
├── late_comparison_capstone.pdf   # Full written report
├── sim_results.csv                # simulated values of estimators, bias, MSE, CI lengths, coverage, etc. 
└── README.md
```

## Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c("AER", "hdm", "glmnet", "ranger",
                   "ggplot2", "dplyr", "tidyr", "parallel"))
```

## Running the simulation

Open `late_simulation.R` in R or RStudio and run the entire script. The runtime is about 15 minutes. The script uses `parLapply` for parallelism; if you are on macOS and see `MallocStackLogging` messages in the console, these are harmless macOS system warnings unrelated to the results.

The script prints four plots and a summary table directly to the console. 

## Author

Aniruddhan Ganesaraman — Department of Statistics and Operations Research, UNC Chapel Hill.
