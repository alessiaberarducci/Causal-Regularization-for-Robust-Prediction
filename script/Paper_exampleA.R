# Paper Example A

# Load estimation, validation and plotting functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/visualizations.R")
source("functions/confidence_intervals.R")

# Define the SEM and simulate the two environments.
set.seed(123)
n <- 100

B <- matrix(
  c(
    0, 0, 0, 1,
    0, 0, -1, 0,
    1, 0, 0, 1,
    0, 0, 0, 0
  ),
  nrow = 4, byrow = TRUE
)
M <- solve(diag(4) - B)
intercept <- c(0, 0, 0, 2)

# Environment 1 has no additive intervention.
epsilon_e1 <- matrix(rnorm(n * 4), ncol = 4)
Z_e1 <- t(M %*% t(sweep(epsilon_e1, 2, intercept, "+")))
colnames(Z_e1) <- c("X1", "X2", "Y", "H")

# Environment 2 shifts X1, X2 and Y, but not H.
epsilon_e2 <- matrix(rnorm(n * 4), ncol = 4)
A_e2 <- cbind(
  A_1 = rnorm(n, mean = 2, sd = 1),
  A_2 = rnorm(n, mean = 1, sd = 1),
  A_Y = rnorm(n, mean = 3, sd = 2),
  A_H = rep(0, n)
)
Z_e2 <- t(M %*% t(sweep(epsilon_e2 + A_e2, 2, intercept, "+")))
colnames(Z_e2) <- c("X1", "X2", "Y", "H")

data <- list(
  Xe = Z_e2[, c("X1", "X2")],
  ye = Z_e2[, "Y"],
  Xo = Z_e1[, c("X1", "X2")],
  yo = Z_e1[, "Y"]
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

# CR is the estimate on the regularization path at the selected gamma.
best_gamma_index <- which.min(abs(gamma_grid - best_gamma))
beta_cr <- beta_gamma[, best_gamma_index]

# Evaluate the empirical worst risk over C_{gamma_CV} at OLS, CV and RL.
point_estimates <- list(
  OLS = drop(beta_ols),
  CV = drop(beta_cr),
  RL = drop(beta_cd)
)
worst_risk_points <- do.call(rbind, lapply(names(point_estimates), function(method) {
  beta <- point_estimates[[method]]
  r_delta_plus <- calc_Rdelta(data, beta)
  r_plus <- in_sample_risk(data, beta)
  data.frame(
    method = method,
    Rplus = r_plus,
    Rdelta_plus = r_delta_plus,
    worst_risk = 0.5 * r_plus + 0.5 * best_gamma * r_delta_plus
  )
}))
row.names(worst_risk_points) <- NULL

# Save both plots; the risk plot focuses on gamma between 0 and 15.
dir.create("visualization", showWarnings = FALSE)
write.csv(
  worst_risk_points,
  "visualization/exampleA_worst_risk_points.csv",
  row.names = FALSE
)
show_plots <- interactive()
plot_cv_risk(
  cv_result, gamma_grid,
  file = "visualization/exampleA_risk.pdf", x_max = 15
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cd, best_gamma,
  beta_true = c(1, 0), file = "visualization/exampleA.pdf", x_max = 30,
  reference_label = "RL", truth_label = "CP"
)
if (show_plots) {
  plot_cv_risk(cv_result, gamma_grid, x_max = 15)
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cd, best_gamma,
    beta_true = c(1, 0), x_max = 30,
    reference_label = "RL", truth_label = "CP"
  )
}

# Add pointwise Monte Carlo intervals to a separate coefficient-path figure.
ci_repetitions <- as.integer(Sys.getenv("PAPER_EXAMPLE_CI_REPS", "500"))
beta_ci <- beta_path_confidence_intervals(
  simulate_example_a_data, n, gamma_grid,
  repetitions = ci_repetitions, seed = 2026
)
plot_beta_path(
  beta_gamma, gamma_grid, beta_cd, best_gamma,
  beta_true = c(1, 0), file = "visualization/exampleA_ci.pdf", x_max = 30,
  beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
  reference_label = "RL", truth_label = "CP"
)
if (show_plots) {
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cd, best_gamma,
    beta_true = c(1, 0), x_max = 30,
    beta_lower = beta_ci$beta_lower, beta_upper = beta_ci$beta_upper,
    reference_label = "RL", truth_label = "CP"
  )
}

# Compare OLS, CR and CD.
cat("Best gamma:", best_gamma, "\n")
print(cbind(OLS = drop(beta_ols), CR = drop(beta_cr), CD = drop(beta_cd)))
cat("Worst risk over C_{gamma_CV} at OLS, CV and RL:\n")
print(worst_risk_points, row.names = FALSE)

# Compare the training intervention with a weaker unseen intervention.
F_A <- crossprod(A_e2) / n
A_tilde <- cbind(
  A_1 = rnorm(n, mean = 1, sd = 0.5),
  A_2 = rnorm(n, mean = 0.5, sd = 0.5),
  A_Y = rnorm(n, mean = 1.5, sd = 1),
  A_H = rep(0, n)
)
F_A_tilde <- crossprod(A_tilde) / n

# Check the eigenvalues of F_A_tilde - gamma_star F_A.
intervention_eigen <- eigen(
  F_A_tilde - best_gamma * F_A,
  symmetric = TRUE
)
cat("Eigenvalues of F_A_tilde - best_gamma * F_A:\n")
print(intervention_eigen$values)

# Generate observations from the unseen environment with shifts A_tilde.
set.seed(456)
n_tilde <- nrow(A_tilde)
epsilon_tilde <- cbind(
  eps_1 = rnorm(n_tilde),
  eps_2 = rnorm(n_tilde),
  eps_Y = rnorm(n_tilde),
  eps_H = rnorm(n_tilde)
)
Z_tilde <- t(
  M %*% t(sweep(epsilon_tilde + A_tilde, 2, intercept, "+"))
)
colnames(Z_tilde) <- c("X1", "X2", "Y", "H")
X_tilde <- Z_tilde[, c("X1", "X2"), drop = FALSE]
y_tilde <- Z_tilde[, "Y"]

# Compare prediction risk in the new environment.
risk_tilde <- data.frame(
  estimator = c("OLS", "CR"),
  beta_1 = c(beta_ols[1], beta_cr[1]),
  beta_2 = c(beta_ols[2], beta_cr[2]),
  risk = c(
    risk(X_tilde, y_tilde, beta_ols),
    risk(X_tilde, y_tilde, beta_cr)
  )
)
cat("Risk in the environment with shifts A_tilde:\n")
print(risk_tilde, row.names = FALSE)
