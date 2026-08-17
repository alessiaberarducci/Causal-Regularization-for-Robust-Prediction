# Run the complete Chamber all-pairs analysis in four readable stages.

# 1. Load settings and data, then center, scale, and spline-expand each pair.
source("script/chamber/chamber_preprocessing.R")

# 2. Fit OLS and causal-regularized estimators and estimate covariance matrices.
source("script/chamber/chamber_estimation.R")

# 3. Select gamma, run grouped Hotelling tests, and summarize each pair.
source("script/chamber/chamber_validation_tests.R")

# 4. Screen all pairs and write the CSV, PDF, and console outputs.
source("script/chamber/chamber_screening.R")
