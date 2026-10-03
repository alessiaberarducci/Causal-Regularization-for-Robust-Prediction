# Simulation models and pointwise Monte Carlo uncertainty bands.

estimate_beta_path_ci <- function(data, gamma_grid) {
  m <- build_matrices(data)
  beta_cd <- drop(compute_cd(m))

  vapply(gamma_grid, function(gamma) {
    drop(solve(
      m$Gplus + gamma * m$Gdelta_plus,
      m$Zplus + gamma * m$Gdelta_plus %*% beta_cd
    ))
  }, numeric(m$p))
}

beta_path_confidence_intervals <- function(
    simulate_data, n, gamma_grid, repetitions = 500,
    level = 0.95, seed = 2026, return_draws = FALSE) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  first_path <- estimate_beta_path_ci(simulate_data(n), gamma_grid)
  p <- nrow(first_path)
  draws <- array(
    NA_real_,
    dim = c(repetitions, p, length(gamma_grid))
  )
  draws[1, , ] <- first_path

  if (repetitions > 1) {
    for (b in 2:repetitions) {
      draws[b, , ] <- estimate_beta_path_ci(
        simulate_data(n), gamma_grid
      )
    }
  }

  alpha <- (1 - level) / 2
  beta_lower <- beta_upper <- matrix(
    NA_real_, nrow = p, ncol = length(gamma_grid)
  )
  for (j in seq_len(p)) {
    bounds <- apply(
      draws[, j, , drop = FALSE],
      3, quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE
    )
    beta_lower[j, ] <- bounds[1, ]
    beta_upper[j, ] <- bounds[2, ]
  }

  result <- list(
    repetitions = repetitions,
    level = level,
    beta_lower = beta_lower,
    beta_upper = beta_upper
  )
  if (return_draws) result$draws <- draws
  result
}

simulate_example_a_data <- function(n) {
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

  epsilon_e1 <- matrix(rnorm(n * 4), ncol = 4)
  Z_e1 <- t(M %*% t(sweep(epsilon_e1, 2, intercept, "+")))

  epsilon_e2 <- matrix(rnorm(n * 4), ncol = 4)
  A_e2 <- cbind(
    rnorm(n, mean = 2, sd = 1),
    rnorm(n, mean = 1, sd = 1),
    rnorm(n, mean = 3, sd = 2),
    rep(0, n)
  )
  Z_e2 <- t(M %*% t(sweep(epsilon_e2 + A_e2, 2, intercept, "+")))

  list(
    Xe = Z_e2[, 1:2, drop = FALSE],
    ye = Z_e2[, 3],
    Xo = Z_e1[, 1:2, drop = FALSE],
    yo = Z_e1[, 3]
  )
}

simulate_example_b_data <- function(n) {
  x1_e1 <- rexp(n)
  y_e1 <- rpois(n, lambda = x1_e1)
  x2_e1 <- -y_e1 + runif(n, -0.5, 0.5)

  x1_e2 <- rpois(n, lambda = 1)
  x2_e2 <- rexp(n)
  y_e2 <- log1p(x1_e2 + x2_e2) + rnorm(n)

  list(
    Xe = cbind(X1 = x1_e2, X2 = x2_e2),
    ye = y_e2,
    Xo = cbind(X1 = x1_e1, X2 = x2_e1),
    yo = y_e1
  )
}

simulate_example_d_data <- function(n, mean_A1 = 3, mean_A2 = 1) {
  B <- matrix(
    c(
      0, 0, 0,
      0, 0, -1,
      1, 0, 0
    ),
    nrow = 3, byrow = TRUE
  )
  M <- solve(diag(3) - B)

  epsilon_e1 <- matrix(rnorm(n * 3), ncol = 3)
  A_e1 <- cbind(
    rnorm(n, mean = mean_A1, sd = 1), rep(0, n), rep(0, n)
  )
  Z_e1 <- t(M %*% t(epsilon_e1 + A_e1))

  epsilon_e2 <- matrix(rnorm(n * 3), ncol = 3)
  A_e2 <- cbind(
    rep(0, n), rnorm(n, mean = mean_A2, sd = 1), rep(0, n)
  )
  Z_e2 <- t(M %*% t(epsilon_e2 + A_e2))

  list(
    Xe = Z_e1[, 1:2, drop = FALSE],
    ye = Z_e1[, 3],
    Xo = Z_e2[, 1:2, drop = FALSE],
    yo = Z_e2[, 3]
  )
}

simulate_example_c_data <- function(n) {
  B <- matrix(
    c(
      0, 0, 0,
      0, 0, -1,
      1, 0, 0
    ),
    nrow = 3, byrow = TRUE
  )
  M <- solve(diag(3) - B)

  epsilon_e1 <- matrix(rnorm(n * 3), ncol = 3)
  A_e1 <- cbind(
    rnorm(n, 0, 1),
    rnorm(n, 0, sqrt(3)),
    rep(0, n)
  )
  Z_e1 <- t(M %*% t(epsilon_e1 + A_e1))

  epsilon_e2 <- matrix(rnorm(n * 3), ncol = 3)
  A_e2 <- cbind(
    rep(0, n),
    rnorm(n, 0, 1),
    rep(0, n)
  )
  Z_e2 <- t(M %*% t(epsilon_e2 + A_e2))

  list(
    Xe = Z_e1[, 1:2, drop = FALSE],
    ye = Z_e1[, 3],
    Xo = Z_e2[, 1:2, drop = FALSE],
    yo = Z_e2[, 3]
  )
}
