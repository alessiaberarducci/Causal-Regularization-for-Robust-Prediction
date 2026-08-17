# Grouped Hotelling T-squared tests for Chamber estimators.

format_chamber_group_test_p_value <- function(value) {
  # Format one p-value for the results table.
  if (is.finite(value)) {
    format.pval(value, digits = 3, eps = 1e-4)
  } else {
    "NA"
  }
}

chamber_covariance_pseudoinverse <- function(covariance) {
  # Compute a stable pseudoinverse of a covariance matrix.
  covariance <- 0.5 * (covariance + t(covariance))
  decomposition <- eigen(covariance, symmetric = TRUE)
  eigenvalue_scale <- max(abs(decomposition$values))
  if (!is.finite(eigenvalue_scale) || eigenvalue_scale == 0) {
    return(list(inverse = matrix(0, nrow(covariance), ncol(covariance)), rank = 0L))
  }
  cutoff <- max(dim(covariance)) * .Machine$double.eps * eigenvalue_scale
  retained <- decomposition$values > cutoff
  if (!any(retained)) {
    return(list(inverse = matrix(0, nrow(covariance), ncol(covariance)), rank = 0L))
  }
  retained_vectors <- decomposition$vectors[, retained, drop = FALSE]
  inverse <- sweep(
    retained_vectors, 2, decomposition$values[retained], "/"
  ) %*% t(retained_vectors)
  list(inverse = inverse, rank = sum(retained))
}

compute_chamber_hotelling_t_squared_tests <- function(
    beta_ols, covariance_ols, beta_cv, covariance_cv,
    coefficient_groups, direct_causes, residual_degrees_of_freedom,
    significance_level = 0.05) {
  # Test whether each predictor's coefficient group is jointly zero.
  parameter_count <- length(beta_ols)
  if (length(beta_cv) != parameter_count ||
      any(dim(covariance_ols) != c(parameter_count, parameter_count)) ||
      any(dim(covariance_cv) != c(parameter_count, parameter_count))) {
    stop("Coefficient and covariance dimensions do not agree.")
  }
  if (is.null(names(coefficient_groups)) || any(!nzchar(names(coefficient_groups)))) {
    stop("coefficient_groups must be a named list of column indices.")
  }
  if (!is.finite(residual_degrees_of_freedom) ||
      residual_degrees_of_freedom <= 0) {
    stop("residual_degrees_of_freedom must be positive.")
  }
  if (!is.finite(significance_level) || significance_level <= 0 ||
      significance_level >= 1) {
    stop("significance_level must be strictly between 0 and 1.")
  }

  estimator_test <- function(beta, covariance, indices) {
    # Run one grouped Hotelling T-squared test.
    coefficient_block <- beta[indices]
    covariance_block <- covariance[indices, indices, drop = FALSE]
    covariance_solution <- chamber_covariance_pseudoinverse(covariance_block)
    effective_terms <- covariance_solution$rank
    if (effective_terms == 0L) {
      return(c(
        statistic = NA_real_, effective_terms = 0,
        f_statistic = NA_real_, p_value = NA_real_
      ))
    }
    statistic <- max(drop(
      t(coefficient_block) %*% covariance_solution$inverse %*%
        coefficient_block
    ), 0)
    f_statistic <- statistic / effective_terms
    c(
      statistic = statistic,
      effective_terms = effective_terms,
      f_statistic = f_statistic,
      p_value = stats::pf(
        f_statistic,
        df1 = effective_terms,
        df2 = residual_degrees_of_freedom,
        lower.tail = FALSE
      )
    )
  }

  result <- do.call(rbind, lapply(names(coefficient_groups), function(variable) {
    indices <- sort(unique(as.integer(coefficient_groups[[variable]])))
    if (!length(indices) || any(indices < 1L | indices > parameter_count)) {
      stop("Invalid coefficient group for ", variable, ".")
    }
    ols <- estimator_test(beta_ols, covariance_ols, indices)
    cv <- estimator_test(beta_cv, covariance_cv, indices)
    data.frame(
      variable = variable,
      terms = length(indices),
      OLS_T2 = ols[["statistic"]],
      OLS_effective_terms = ols[["effective_terms"]],
      OLS_F = ols[["f_statistic"]],
      OLS_p = ols[["p_value"]],
      CV_T2 = cv[["statistic"]],
      CV_effective_terms = cv[["effective_terms"]],
      CV_F = cv[["f_statistic"]],
      CV_p = cv[["p_value"]],
      is_direct_cause = variable %in% direct_causes,
      row.names = NULL
    )
  }))

  result$OLS_reject <- result$OLS_p < significance_level
  result$CV_reject <- result$CV_p < significance_level

  # Order predictors by the strongest evidence against the null hypothesis.
  result <- result[order(result$CV_p, result$OLS_p), , drop = FALSE]
  attr(result, "significance_level") <- significance_level
  attr(result, "residual_degrees_of_freedom") <- residual_degrees_of_freedom
  result
}

write_chamber_hotelling_t_squared_tests_pdf <- function(
    table, file, gamma, model_description) {
  # Write the grouped OLS and CR test results to a PDF table.
  significance_level <- attr(table, "significance_level")
  residual_degrees_of_freedom <- attr(table, "residual_degrees_of_freedom")
  if (is.null(significance_level) || is.null(residual_degrees_of_freedom)) {
    stop("The grouped-test table is missing its test metadata.")
  }
  if (length(model_description) != 1L || !nzchar(model_description)) {
    stop("model_description must be a non-empty string.")
  }

  grDevices::pdf(
    file, width = 14, height = max(8, 0.42 * nrow(table) + 3.2),
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0.8, 0.5, 3.1, 0.5), family = "serif")
  graphics::plot.new()
  row_count <- nrow(table)
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, row_count + 2.5))
  column_x <- c(0.02, 0.15, 0.28, 0.40, 0.53, 0.66, 0.79, 0.97)
  alignment <- c(0, rep(1, length(column_x) - 1L))
  headers <- c(
    "Covariate", "Terms", "OLS T2", "OLS p", "OLS H0",
    "CV T2", "CV p", "CV H0"
  )
  header_y <- row_count + 1.4
  for (column_index in seq_along(column_x)) {
    graphics::text(
      column_x[column_index], header_y, headers[column_index],
      adj = c(alignment[column_index], 0.5), font = 2, cex = 0.68
    )
  }
  graphics::segments(0.015, header_y - 0.48, 0.985, header_y - 0.48)

  for (row_index in seq_len(row_count)) {
    current <- table[row_index, ]
    row_y <- row_count - row_index + 0.55
    row_colour <- if (current$is_direct_cause) "red3" else "black"
    values <- c(
      current$variable,
      current$terms,
      formatC(current$OLS_T2, format = "f", digits = 3),
      format_chamber_group_test_p_value(current$OLS_p),
      if (current$OLS_reject) "reject" else "do not reject",
      formatC(current$CV_T2, format = "f", digits = 3),
      format_chamber_group_test_p_value(current$CV_p),
      if (current$CV_reject) "reject" else "do not reject"
    )
    for (column_index in seq_along(column_x)) {
      graphics::text(
        column_x[column_index], row_y, values[column_index],
        adj = c(alignment[column_index], 0.5),
        col = row_colour, cex = 0.63
      )
    }
  }

  graphics::mtext(
    paste0(
      "Hotelling T squared: OLS versus CV (", model_description,
      ", gamma_CV = ", format(gamma, digits = 6), ")"
    ),
    side = 3, line = 1.45, font = 2, cex = 1.05
  )
  graphics::mtext(
    paste0(
      "T2 = beta[J]' V[J,J]^{-1} beta[J]; F approximation = T2 / rank(V[J,J]); ",
      "df2 = ", residual_degrees_of_freedom,
      "; decisions use alpha = ", format(significance_level, digits = 3),
      "; direct causes in red"
    ),
    side = 3, line = 0.15, cex = 0.68
  )
}
