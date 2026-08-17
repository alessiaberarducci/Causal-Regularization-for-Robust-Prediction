# Paper Example B

# Load estimation, validation and plotting functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/visualizations.R")
source("functions/confidence_intervals.R")

# Simulate the two environments.
set.seed(134)
n <- 100

x1_e1 <- rexp(n)
y_e1 <- rpois(n, lambda = x1_e1)
x2_e1 <- -y_e1 + runif(n, -0.5, 0.5)

x1_e2 <- rpois(n, lambda = 1)
y_e2 <- log1p(x1_e2) + runif(n, -1, 1)
x2_e2 <- -0.5 * y_e2 + runif(n, -0.5, 0.5)

data <- list(
  Xe = cbind(X1 = x1_e2, X2 = x2_e2),
  ye = y_e2,
  Xo = cbind(X1 = x1_e1, X2 = x2_e1),
  yo = y_e1
)

# Select gamma by cross-validation.
gamma_grid <- seq(0, 300, by = 0.1)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
cv_result <- cross_validation(5, data, estimators)
best_gamma <- as.numeric(cv_result$best_gamma)

# Refit on the full sample and compute the regularization path.
m <- build_matrices(data)
beta_cp <- compute_cd(m)
beta_ols <- compute_ols(m)
beta_gamma <- vapply(
  gamma_grid,
  function(gamma) drop(make_estimator(gamma)(data)),
  numeric(m$p)
)

# CR is the estimate on the regularization path at the selected gamma.
best_gamma_index <- which.min(abs(gamma_grid - best_gamma))
beta_cr <- beta_gamma[, best_gamma_index]

# Save both plots; the risk plot focuses on gamma between 0 and 15.
dir.create("visualization", showWarnings = FALSE)
show_plots <- interactive()
plot_cv_risk(
  cv_result, gamma_grid,
  file = "visualization/exampleB_risk.pdf", x_max = 15
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cp, best_gamma,
  file = "visualization/exampleB.pdf", x_max = 50,
  reference_label = "RL"
)
if (show_plots) {
  plot_cv_risk(cv_result, gamma_grid, x_max = 15)
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cp, best_gamma, x_max = 50,
    reference_label = "RL"
  )
}

# Add pointwise Monte Carlo intervals to a separate coefficient-path figure.
ci_repetitions <- as.integer(Sys.getenv("PAPER_EXAMPLE_CI_REPS", "500"))
beta_ci <- beta_path_confidence_intervals(
  simulate_example_b_data, n, gamma_grid,
  repetitions = ci_repetitions, seed = 2027
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cp, best_gamma,
  file = "visualization/exampleB_ci.pdf", x_max = 50,
  beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
  reference_label = "RL"
)
if (show_plots) {
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cp, best_gamma, x_max = 50,
    beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
    reference_label = "RL"
  )
}

# Report the selected gamma and compare the three estimates.
cat("Best gamma:", best_gamma, "\n")
print(cbind(OLS = drop(beta_ols), CR = drop(beta_cr), CD = drop(beta_cp)))
