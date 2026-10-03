# RNA preprocessing, estimation, and prediction risk.

filter_environments <- function(data, environment_column, minimum_cells = 100) {
  environment <- as.character(data[[environment_column]])
  counts <- table(environment)
  keep <- names(counts)[counts >= minimum_cells]
  data[environment %in% setdiff(keep, "excluded"), , drop = FALSE]
}

make_two_environment_data <- function(
    data, environment_column, first_environments, second_environments,
    predictor_genes, target_gene) {
  overlap <- intersect(first_environments, second_environments)
  if (length(overlap) > 0L) {
    stop("The two training environments contain overlapping interventions.")
  }

  available <- unique(as.character(data[[environment_column]]))
  missing <- setdiff(c(first_environments, second_environments), available)
  if (length(missing) > 0L) {
    stop("Training interventions not found: ", paste(missing, collapse = ", "))
  }

  first <- data[
    data[[environment_column]] %in% first_environments, , drop = FALSE
  ]
  second <- data[
    data[[environment_column]] %in% second_environments, , drop = FALSE
  ]
  list(
    data = list(
      Xe = as.matrix(first[, predictor_genes, drop = FALSE]),
      ye = as.numeric(first[[target_gene]]),
      Xo = as.matrix(second[, predictor_genes, drop = FALSE]),
      yo = as.numeric(second[[target_gene]])
    ),
    first = first,
    second = second,
    first_environments = first_environments,
    second_environments = second_environments
  )
}

compute_common_centering <- function(data) {
  x_center <- 0.5 * (
    colMeans(data$Xe) + colMeans(data$Xo)
  )
  y_center <- 0.5 * (mean(data$ye) + mean(data$yo))
  list(
    data = list(
      Xe = sweep(data$Xe, 2, x_center, "-"),
      ye = data$ye - y_center,
      Xo = sweep(data$Xo, 2, x_center, "-"),
      yo = data$yo - y_center
    ),
    x_center = x_center,
    y_center = y_center
  )
}

apply_common_centering <- function(data, centering) {
  list(
    Xe = sweep(data$Xe, 2, centering$x_center, "-"),
    ye = data$ye - centering$y_center,
    Xo = sweep(data$Xo, 2, centering$x_center, "-"),
    yo = data$yo - centering$y_center
  )
}

centered_cross_validation <- function(data, estimators, fold_assignment) {
  fold_e <- fold_assignment$environment_e
  fold_o <- fold_assignment$environment_o
  folds <- length(unique(fold_e))
  gamma_values <- as.numeric(names(estimators))
  foldwise <- matrix(
    NA_real_, nrow = folds, ncol = length(gamma_values),
    dimnames = list(NULL, names(estimators))
  )

  for (fold in seq_len(folds)) {
    raw_training <- subset_data(data, fold_e != fold, fold_o != fold)
    raw_validation <- subset_data(data, fold_e == fold, fold_o == fold)
    centering <- compute_common_centering(raw_training)
    training <- centering$data
    validation <- apply_common_centering(raw_validation, centering)

    training_moments <- build_matrices(training)
    beta_cp <- drop(compute_cd(training_moments))
    beta_path <- vapply(gamma_values, function(gamma) {
      drop(solve(
        training_moments$Gplus +
          gamma * training_moments$Gdelta_plus,
        training_moments$Zplus +
          gamma * training_moments$Gdelta_plus %*% beta_cp
      ))
    }, numeric(training_moments$p))

    validation_moments <- build_matrices(validation)
    validation_cp <- drop(compute_cd(validation_moments))
    beta_difference <- sweep(beta_path, 1, validation_cp, "-")
    foldwise[fold, ] <- colSums(
      beta_difference *
        (validation_moments$Gdelta_plus %*% beta_difference)
    )
  }

  mean_loss <- colMeans(foldwise)
  best <- which.min(mean_loss)
  list(
    best_gamma = names(estimators)[best],
    mean_loss = mean_loss,
    foldwise_loss = foldwise,
    fold_assignment = fold_assignment
  )
}

recover_intercept <- function(beta, centering) {
  drop(centering$y_center - crossprod(centering$x_center, beta))
}

compute_environment_risk <- function(
    data, beta, predictor_genes, target_gene, intercept = 0) {
  X <- as.matrix(data[, predictor_genes, drop = FALSE])
  y <- as.numeric(data[[target_gene]])
  complete <- complete.cases(X, y)
  if (!any(complete)) return(NA_real_)
  prediction <- intercept +
    X[complete, , drop = FALSE] %*% beta[predictor_genes]
  mean((y[complete] - prediction)^2)
}

compute_test_risks <- function(
    data, environments, environment_column, beta_ols, beta_cv,
    predictor_genes, target_gene, intercept_ols = 0, intercept_cv = 0) {
  result <- lapply(environments, function(environment) {
    current <- data[data[[environment_column]] == environment, , drop = FALSE]
    data.frame(
      environment = environment,
      n_cells = nrow(current),
      OLS = compute_environment_risk(
        current, beta_ols, predictor_genes, target_gene, intercept_ols
      ),
      CV = compute_environment_risk(
        current, beta_cv, predictor_genes, target_gene, intercept_cv
      )
    )
  })
  result <- do.call(rbind, result)
  result[is.finite(result$OLS) & is.finite(result$CV), , drop = FALSE]
}

compare_training_risks <- function(data, estimates, gamma) {
  rows <- lapply(names(estimates), function(method) {
    beta <- estimates[[method]]
    r_delta <- calc_Rdelta(data, beta)
    r_plus <- in_sample_risk(data, beta)
    data.frame(
      method = method,
      Rdelta_plus = r_delta,
      Rplus = r_plus,
      Rgamma = 0.5 * r_plus + 0.5 * gamma * r_delta
    )
  })
  do.call(rbind, rows)
}
