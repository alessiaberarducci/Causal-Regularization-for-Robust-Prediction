# Cross-validation for generalized causal regularization

# Assign balanced fold labels in random order.
split_idx <- function(k, n) {
  if (k < 2L) stop("The number of folds must be at least 2.")
  if (n < k) stop("The number of folds cannot exceed the sample size.")
  sample(rep(seq_len(k), length.out = n))
}

subset_data <- function(data, rows_e, rows_o) {
  # The environments are split independently because their sizes can differ.
  list(
    Xe = data$Xe[rows_e, , drop = FALSE],
    ye = data$ye[rows_e],
    Xo = data$Xo[rows_o, , drop = FALSE],
    yo = data$yo[rows_o]
  )
}

build_matrices <- function(data) {
  # Replace the eigenvalues of Gdelta by their absolute values.
  m <- moments(data)
  eig <- eigen(m$Gdelta, symmetric = TRUE)
  m$Gdelta_plus <- eig$vectors %*%
    diag(abs(eig$values), nrow = length(eig$values)) %*%
    t(eig$vectors)
  m
}

estimate_beta_cp <- function(data) {
  compute_cd(moments(data))
}

calc_Rdelta <- function(data, beta_hat) {
  # Validation risk is measured relative to the validation CP estimate.
  m <- build_matrices(data)
  beta_diff <- drop(beta_hat) - drop(estimate_beta_cp(data))
  drop(t(beta_diff) %*% m$Gdelta_plus %*% beta_diff)
}

make_estimator <- function(gamma) {
  # Return an estimator fixed at the current value of gamma.
  force(gamma)
  function(data) {
    m <- build_matrices(data)
    beta_cp <- estimate_beta_cp(data)
    solve(
      m$Gplus + gamma * m$Gdelta_plus,
      m$Zplus + gamma * m$Gdelta_plus %*% beta_cp
    )
  }
}

fit_cross_validation <- function(folds, estimator, training, validation) {
  vapply(seq_len(folds), function(i) {
    calc_Rdelta(validation[[i]], estimator(training[[i]]))
  }, numeric(1))
}

validate_cv_data <- function(data, folds) {
  # Check dimensions before constructing any folds.
  required <- c("Xe", "ye", "Xo", "yo")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing data objects: ", paste(missing, collapse = ", "))
  }

  n_e <- nrow(data$Xe)
  n_o <- nrow(data$Xo)
  if (length(data$ye) != n_e) stop("nrow(Xe) must equal length(ye).")
  if (length(data$yo) != n_o) stop("nrow(Xo) must equal length(yo).")
  if (ncol(data$Xe) != ncol(data$Xo)) {
    stop("Xe and Xo must have the same number of predictors.")
  }
  if (folds > min(n_e, n_o)) {
    stop("The number of folds cannot exceed the smaller environment.")
  }
  invisible(c(environment_e = n_e, environment_o = n_o))
}

cross_validation <- function(folds, data, estimators) {
  # Create separate fold assignments for the two environments.
  sample_sizes <- validate_cv_data(data, folds)
  fold_e <- split_idx(folds, sample_sizes[["environment_e"]])
  fold_o <- split_idx(folds, sample_sizes[["environment_o"]])

  training <- validation <- vector("list", folds)
  for (i in seq_len(folds)) {
    training[[i]] <- subset_data(data, fold_e != i, fold_o != i)
    validation[[i]] <- subset_data(data, fold_e == i, fold_o == i)
  }

  # Compare gamma values using the average validation criterion.
  foldwise <- vapply(
    estimators,
    fit_cross_validation,
    folds = folds,
    training = training,
    validation = validation,
    FUN.VALUE = numeric(folds)
  )
  mean_risk <- colMeans(foldwise)
  best <- which.min(mean_risk)

  list(
    best_gamma = names(estimators)[best],
    mean_abs_rdelta = mean_risk,
    foldwise_rdelta = foldwise,
    fold_assignment = list(environment_e = fold_e, environment_o = fold_o),
    sample_sizes = sample_sizes
  )
}
