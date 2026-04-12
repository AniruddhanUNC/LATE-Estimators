# Benchmarking LATE Estimators

This repository contains the code and report for a graduate capstone project comparing five instrumental variables estimators for the Local Average Treatment Effect (LATE).

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
├── late_final_data.csv            # simulated values of estimators, bias, MSE, CI lengths, coverage, etc. 
└── README.md
```

## Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c("AER", "hdm", "glmnet", "ranger",
                   "ggplot2", "dplyr", "tidyr", "parallel"))
```

## Running the simulation

Open `late_simulation.R` in R or RStudio and run the entire script. The run time can be over 2 hours. The script uses `parLapply` for parallelism; if you are on macOS and see `MallocStackLogging` messages in the console, these are harmless macOS system warnings unrelated to the results.

The script prints two plots and a summary table directly to the console. 

## Reference

This project is a capstone for a graduate course in Causal Machine Learning at the University of North Carolina at Chapel Hill. The primary references are:

- Imbens & Angrist (1994). *Econometrica*.
- Belloni, Chernozhukov & Hansen (2014). *Review of Economic Studies*.
- Chernozhukov et al. (2018). *Econometrics Journal*.
- Tan (2006). *Journal of the American Statistical Association*.
- van der Laan & Rose (2011). *Targeted Learning*. Springer.

## Author

Aniruddhan Ganesaraman — Department of Statistics and Operations Research, UNC Chapel Hill.
