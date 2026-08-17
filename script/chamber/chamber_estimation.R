# Chamber estimation: compute the regularization path, coefficients, and covariance.

gamma_solution_matrix <- function(matrices, beta_cp) {
  # Compute the full regularization path efficiently by eigendecomposition.
  Gplus <- matrices$Gplus
  Gdelta_plus <- matrices$Gdelta_plus
  H <- 0.5 * (
    Gplus + Gdelta_plus + t(Gplus + Gdelta_plus)
  )
  H_decomposition <- eigen(H, symmetric = TRUE)
  H_scale <- max(abs(H_decomposition$values))
  retained <- H_decomposition$values > pseudoinverse_tolerance * H_scale
  if (!any(retained)) stop("Gplus + Gdelta_plus has numerical rank zero.")
  H_inverse_square_root <- sweep(
    H_decomposition$vectors[, retained, drop = FALSE],
    2L,
    sqrt(H_decomposition$values[retained]),
    "/"
  )
  transformed_delta <-
    t(H_inverse_square_root) %*% Gdelta_plus %*% H_inverse_square_root
  transformed_delta <- 0.5 * (
    transformed_delta + t(transformed_delta)
  )
  delta_decomposition <- eigen(transformed_delta, symmetric = TRUE)
  delta_values <- pmin(pmax(delta_decomposition$values, 0), 1)
  projection <- H_inverse_square_root %*% delta_decomposition$vectors

  right_hand_side_0 <- drop(matrices$Zplus)
  right_hand_side_1 <- drop(Gdelta_plus %*% beta_cp)
  solutions <- matrix(
    NA_real_, nrow = nrow(Gplus), ncol = length(gamma_grid)
  )
  zero_index <- which(gamma_grid == 0)
  solutions[, zero_index] <- drop(symmetric_pseudoinverse_solve(
    Gplus, right_hand_side_0
  ))

  positive_indices <- which(gamma_grid > 0)
  positive_gamma <- gamma_grid[positive_indices]
  right_hand_sides <- outer(
    right_hand_side_0, rep(1, length(positive_gamma))
  ) + outer(right_hand_side_1, positive_gamma)
  projected_right_hand_sides <- t(projection) %*% right_hand_sides
  denominators <- 1 + outer(delta_values, positive_gamma - 1)
  solutions[, positive_indices] <- projection %*% (
    projected_right_hand_sides / denominators
  )
  solutions
}

fit_estimators_and_tests <- function(model_data, gamma) {
  # Fit OLS and robust estimators, estimate covariance, and run group tests.

  # Build pooled and difference moment matrices.
  matrices <- build_matrices(model_data)

  # Estimate the robust reference, pooled OLS, and selected CR coefficients.
  beta_cp <- symmetric_pseudoinverse_solve(
    matrices$Gdelta, matrices$Zdelta
  )
  beta_ols <- drop(symmetric_pseudoinverse_solve(
    matrices$Gplus, matrices$Zplus
  ))
  beta_cv <- drop(symmetric_pseudoinverse_solve(
    matrices$Gplus + gamma * matrices$Gdelta_plus,
    matrices$Zplus + gamma * matrices$Gdelta_plus %*% beta_cp
  ))

  identity <- diag(ncol(model_data$Xe))
  Gplus_inverse <- symmetric_pseudoinverse_solve(
    matrices$Gplus, identity
  )
  Gdelta_inverse <- symmetric_pseudoinverse_solve(
    matrices$Gdelta, identity
  )
  Ggamma <- matrices$Gplus + gamma * matrices$Gdelta_plus
  Ggamma_inverse <- symmetric_pseudoinverse_solve(Ggamma, identity)
  A_gamma <- gamma * matrices$Gdelta_plus %*% Gdelta_inverse

  residual_e_ols <- drop(model_data$ye - model_data$Xe %*% beta_ols)
  residual_o_ols <- drop(model_data$yo - model_data$Xo %*% beta_ols)

  # Estimate the OLS coefficient covariance matrix from residual variances.
  variance_Ze_ols <- stats::var(residual_e_ols) *
    crossprod(model_data$Xe) / nrow(model_data$Xe)^2
  variance_Zo_ols <- stats::var(residual_o_ols) *
    crossprod(model_data$Xo) / nrow(model_data$Xo)^2
  variance_beta_ols <- Gplus_inverse %*%
    (variance_Ze_ols + variance_Zo_ols) %*% t(Gplus_inverse)
  variance_beta_ols <- 0.5 * (
    variance_beta_ols + t(variance_beta_ols)
  )

  residual_e_cv <- drop(model_data$ye - model_data$Xe %*% beta_cv)
  residual_o_cv <- drop(model_data$yo - model_data$Xo %*% beta_cv)

  # Estimate the robust coefficient covariance matrix.
  variance_Ze_cv <- stats::var(residual_e_cv) *
    crossprod(model_data$Xe) / nrow(model_data$Xe)^2
  variance_Zo_cv <- stats::var(residual_o_cv) *
    crossprod(model_data$Xo) / nrow(model_data$Xo)^2
  variance_Zgamma <-
    (identity + A_gamma) %*% variance_Ze_cv %*%
      t(identity + A_gamma) +
    (identity - A_gamma) %*% variance_Zo_cv %*%
      t(identity - A_gamma)
  variance_beta_cv <- Ggamma_inverse %*%
    variance_Zgamma %*% t(Ggamma_inverse)
  variance_beta_cv <- 0.5 * (
    variance_beta_cv + t(variance_beta_cv)
  )

  coefficient_groups <- attr(model_data, "coefficient_groups")
  if (is.null(coefficient_groups)) {
    coefficient_groups <- setNames(
      lapply(seq_along(predictors), function(index) {
        ((index - 1L) * terms_per_predictor + 1L):
          (index * terms_per_predictor)
      }),
      predictors
    )
  }
  design_rank <- qr(rbind(model_data$Xe, model_data$Xo))$rank

  # Compute residual degrees of freedom for the F approximation.
  residual_degrees_of_freedom <-
    nrow(model_data$Xe) + nrow(model_data$Xo) - design_rank
  test_table <- compute_chamber_hotelling_t_squared_tests(
    beta_ols = beta_ols,
    covariance_ols = variance_beta_ols,
    beta_cv = beta_cv,
    covariance_cv = variance_beta_cv,
    coefficient_groups = coefficient_groups,
    direct_causes = direct_causes,
    residual_degrees_of_freedom = residual_degrees_of_freedom,
    significance_level = significance_level
  )

  list(
    table = test_table,
    residual_degrees_of_freedom = residual_degrees_of_freedom
  )
}
