# Paper Example D

# Load estimation, validation and plotting functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/visualizations.R")
source("functions/confidence_intervals.R")

# Structural equations:
# X1 <- epsilon_1 + shift_1
# Y  <- X1 + epsilon_Y
# X2 <- -Y + shift_2 + epsilon_2
set.seed(123)
n <- 1000

B <- matrix(
  c(
    0, 0, 0,
    0, 0, -1,
    1, 0, 0
  ),
  nrow = 3, byrow = TRUE
)
M <- solve(diag(3) - B)

# Environment e1 has an additive shift only on X1.
epsilon_e1 <- matrix(rnorm(n * 3), ncol = 3)
A_e1 <- cbind(
  A_1 = rnorm(n, mean = 0, sd = 1),
  A_2 = rep(0, n),
  A_Y = rep(0, n)
)
Z_e1 <- t(M %*% t(epsilon_e1 + A_e1))
colnames(Z_e1) <- c("X1", "X2", "Y")

# Environment e2 has an additive shift only on X2.
epsilon_e2 <- matrix(rnorm(n * 3), ncol = 3)
A_e2 <- cbind(
  A_1 = rep(0, n),
  A_2 = rnorm(n, mean = 0, sd = 1),
  A_Y = rep(0, n)
)
Z_e2 <- t(M %*% t(epsilon_e2 + A_e2))
colnames(Z_e2) <- c("X1", "X2", "Y")

data <- list(
  Xe = Z_e1[, c("X1", "X2")],
  ye = Z_e1[, "Y"],
  Xo = Z_e2[, c("X1", "X2")],
  yo = Z_e2[, "Y"]
)

# Select gamma by cross-validation.
gamma_grid <- seq(0, 30, by = 0.1)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
cv_result <- cross_validation(5, data, estimators)
best_gamma <- as.numeric(cv_result$best_gamma)

# Refit on the full sample and compute the regularization path.
m <- build_matrices(data)
beta_cd <- compute_cd(m)
beta_ols <- compute_ols(m)
beta_gamma <- vapply(
  gamma_grid,
  function(gamma) drop(make_estimator(gamma)(data)),
  numeric(m$p)
)

# CR is the estimate on the path at the selected gamma.
best_gamma_index <- which.min(abs(gamma_grid - best_gamma))
beta_cr <- beta_gamma[, best_gamma_index]

# Save the cross-validation and coefficient-path plots.
dir.create("visualization", showWarnings = FALSE)
show_plots <- interactive()
plot_cv_risk(
  cv_result, gamma_grid,
  file = "visualization/exampleD_risk.pdf", x_max = 15
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cd, best_gamma,
  beta_true = c(1, 0), file = "visualization/exampleD.pdf", x_max = 30,
  reference_label = "CP", truth_label = "CP"
)
if (show_plots) {
  plot_cv_risk(cv_result, gamma_grid, x_max = 15)
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cd, best_gamma,
    beta_true = c(1, 0), x_max = 30,
    reference_label = "CP", truth_label = "CP"
  )
}

# Add pointwise Monte Carlo intervals to a separate coefficient-path figure.
ci_repetitions <- as.integer(Sys.getenv("PAPER_EXAMPLE_CI_REPS", "500"))
beta_ci <- beta_path_confidence_intervals(
  simulate_example_d_data, n, gamma_grid,
  repetitions = ci_repetitions, seed = 2028
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cd, best_gamma,
  beta_true = c(1, 0), file = "visualization/exampleD_ci.pdf", x_max = 30,
  beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
  reference_label = "CP", truth_label = "CP"
)
if (show_plots) {
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cd, best_gamma,
    beta_true = c(1, 0), x_max = 30,
    beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
    reference_label = "CP", truth_label = "CP"
  )
}

# Compare OLS, CR and CD with the true coefficient.
cat("Best gamma:", best_gamma, "\n")
print(cbind(
  OLS = drop(beta_ols),
  CR = drop(beta_cr),
  CD = drop(beta_cd),
  true = c(1, 0)
))
