# Chamber cross-validation and classification metrics.

fit_basis_cross_validation_exact <- function(
    model_data, candidate_indices = seq_along(gamma_grid)) {
  candidate_indices <- sort(unique(as.integer(candidate_indices)))
  if (!length(candidate_indices) ||
      any(candidate_indices < 1L | candidate_indices > length(gamma_grid))) {
    stop("Invalid exact-CV gamma candidate indices.")
  }
  candidate_gamma <- gamma_grid[candidate_indices]
  set.seed(123)
  fold_e <- split_idx(fold_count, nrow(model_data$Xe))
  fold_o <- split_idx(fold_count, nrow(model_data$Xo))
  training_components <- validation_components <- vector("list", fold_count)

  for (fold_index in seq_len(fold_count)) {
    training_data <- subset_data(
      model_data, fold_e != fold_index, fold_o != fold_index
    )
    validation_data <- subset_data(
      model_data, fold_e == fold_index, fold_o == fold_index
    )
    training_matrices <- build_matrices(training_data)
    validation_matrices <- build_matrices(validation_data)
    training_components[[fold_index]] <- list(
      matrices = training_matrices,
      beta_cp = symmetric_pseudoinverse_solve(
        training_matrices$Gdelta, training_matrices$Zdelta
      )
    )
    validation_components[[fold_index]] <- list(
      Gdelta_plus = validation_matrices$Gdelta_plus,
      beta_cp = drop(symmetric_pseudoinverse_solve(
        validation_matrices$Gdelta, validation_matrices$Zdelta
      ))
    )
  }

  mean_validation_risk <- vapply(candidate_gamma, function(gamma) {
    fold_risk <- vapply(seq_len(fold_count), function(fold_index) {
      training_component <- training_components[[fold_index]]
      training_matrices <- training_component$matrices
      beta_hat <- drop(symmetric_pseudoinverse_solve(
        training_matrices$Gplus + gamma * training_matrices$Gdelta_plus,
        training_matrices$Zplus + gamma *
          training_matrices$Gdelta_plus %*% training_component$beta_cp
      ))
      validation_component <- validation_components[[fold_index]]
      beta_difference <- beta_hat - validation_component$beta_cp
      drop(
        t(beta_difference) %*% validation_component$Gdelta_plus %*%
          beta_difference
      )
    }, numeric(1))
    mean(fold_risk)
  }, numeric(1))

  best_index <- which.min(mean_validation_risk)
  list(
    gamma = candidate_gamma[best_index],
    validation_risk = mean_validation_risk[best_index],
    grid_indices = candidate_indices,
    candidate_risk = mean_validation_risk
  )
}

fit_basis_cross_validation_fast <- function(model_data) {
  set.seed(123)
  fold_e <- split_idx(fold_count, nrow(model_data$Xe))
  fold_o <- split_idx(fold_count, nrow(model_data$Xo))
  mean_validation_risk <- numeric(length(gamma_grid))

  for (fold_index in seq_len(fold_count)) {
    training_data <- subset_data(
      model_data, fold_e != fold_index, fold_o != fold_index
    )
    validation_data <- subset_data(
      model_data, fold_e == fold_index, fold_o == fold_index
    )
    training_matrices <- build_matrices(training_data)
    validation_matrices <- build_matrices(validation_data)
    training_beta_cp <- symmetric_pseudoinverse_solve(
      training_matrices$Gdelta, training_matrices$Zdelta
    )
    validation_beta_cp <- drop(symmetric_pseudoinverse_solve(
      validation_matrices$Gdelta, validation_matrices$Zdelta
    ))
    beta_by_gamma <- gamma_solution_matrix(
      training_matrices, training_beta_cp
    )
    beta_difference <- sweep(
      beta_by_gamma, 1L, validation_beta_cp, "-"
    )
    mean_validation_risk <- mean_validation_risk + colSums(
      beta_difference * (validation_matrices$Gdelta_plus %*% beta_difference)
    ) / fold_count
  }

  best_index <- which.min(mean_validation_risk)
  list(
    gamma = gamma_grid[best_index],
    validation_risk = mean_validation_risk[best_index],
    mean_validation_risk = mean_validation_risk
  )
}

fit_basis_cross_validation_hybrid <- function(model_data) {
  fast_result <- fit_basis_cross_validation_fast(model_data)
  coarse_indices <- seq(1L, length(gamma_grid), by = 10L)
  fast_candidate_indices <- head(
    order(fast_result$mean_validation_risk), 20L
  )
  initial_indices <- sort(unique(c(
    coarse_indices, fast_candidate_indices,
    which(gamma_grid == 0), length(gamma_grid)
  )))
  initial_exact <- fit_basis_cross_validation_exact(
    model_data, initial_indices
  )
  refinement_indices <- which(
    abs(gamma_grid - initial_exact$gamma) <= 0.1 + 1e-12
  )
  new_refinement_indices <- setdiff(refinement_indices, initial_indices)

  if (!length(new_refinement_indices)) return(initial_exact)
  refined_exact <- fit_basis_cross_validation_exact(
    model_data, new_refinement_indices
  )
  combined_indices <- c(
    initial_exact$grid_indices, refined_exact$grid_indices
  )
  combined_risk <- c(
    initial_exact$candidate_risk, refined_exact$candidate_risk
  )
  best <- which.min(combined_risk)
  list(
    gamma = gamma_grid[combined_indices[best]],
    validation_risk = combined_risk[best],
    grid_indices = combined_indices,
    candidate_risk = combined_risk
  )
}

use_exact_cv <- tolower(Sys.getenv(
  "CHAMBER_PAIR_EXACT_CV", "false"
)) %in% c("true", "1", "yes")
fit_basis_cross_validation <- if (use_exact_cv) {
  fit_basis_cross_validation_exact
} else {
  fit_basis_cross_validation_hybrid
}

collapse_variables <- function(values) {
  if (length(values)) paste(values, collapse = ";") else ""
}

classification_metrics <- function(true_positives, false_positives) {
  false_negatives <- length(direct_causes) - true_positives
  true_negatives <- length(non_causes) - false_positives
  precision_denominator <- true_positives + false_positives
  mcc_denominator <- sqrt(
    (true_positives + false_positives) *
      (true_positives + false_negatives) *
      (true_negatives + false_positives) *
      (true_negatives + false_negatives)
  )

  sensitivity <- true_positives / length(direct_causes)
  specificity <- true_negatives / length(non_causes)
  precision <- if (precision_denominator > 0) {
    true_positives / precision_denominator
  } else {
    NA_real_
  }
  f1_score <- 2 * true_positives /
    (2 * true_positives + false_positives + false_negatives)
  mcc <- if (mcc_denominator > 0) {
    (true_positives * true_negatives -
       false_positives * false_negatives) / mcc_denominator
  } else {
    NA_real_
  }

  c(
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    F1 = f1_score,
    balanced_accuracy = (sensitivity + specificity) / 2,
    MCC = mcc
  )
}

summarize_tests <- function(
    pair_id, file_1, file_2, gamma, validation_risk, test_table) {
  direct_table <- test_table[test_table$is_direct_cause, , drop = FALSE]
  non_cause_table <- test_table[!test_table$is_direct_cause, , drop = FALSE]
  OLS_true_positives <- sum(direct_table$OLS_reject)
  CV_true_positives <- sum(direct_table$CV_reject)
  OLS_false_positives <- sum(non_cause_table$OLS_reject)
  CV_false_positives <- sum(non_cause_table$CV_reject)
  OLS_metrics <- classification_metrics(
    OLS_true_positives, OLS_false_positives
  )
  CV_metrics <- classification_metrics(
    CV_true_positives, CV_false_positives
  )
  same_sensitivity <- OLS_true_positives == CV_true_positives
  same_direct_decisions <- identical(
    direct_table$OLS_reject, direct_table$CV_reject
  )
  improved_specificity <- CV_false_positives < OLS_false_positives

  data.frame(
    pair_id = pair_id,
    environment_1 = file_1,
    environment_2 = file_2,
    n_1 = NA_integer_,
    n_2 = NA_integer_,
    gamma_CV = gamma,
    CV_validation_Rdelta_plus = validation_risk,
    OLS_direct_rejected = OLS_true_positives,
    CV_direct_rejected = CV_true_positives,
    sensitivity_OLS = unname(OLS_metrics["sensitivity"]),
    sensitivity_CV = unname(CV_metrics["sensitivity"]),
    OLS_noncauses_rejected = OLS_false_positives,
    CV_noncauses_rejected = CV_false_positives,
    specificity_OLS = unname(OLS_metrics["specificity"]),
    specificity_CV = unname(CV_metrics["specificity"]),
    precision_OLS = unname(OLS_metrics["precision"]),
    precision_CV = unname(CV_metrics["precision"]),
    F1_OLS = unname(OLS_metrics["F1"]),
    F1_CV = unname(CV_metrics["F1"]),
    balanced_accuracy_OLS = unname(OLS_metrics["balanced_accuracy"]),
    balanced_accuracy_CV = unname(CV_metrics["balanced_accuracy"]),
    MCC_OLS = unname(OLS_metrics["MCC"]),
    MCC_CV = unname(CV_metrics["MCC"]),
    specificity_gain_count = OLS_false_positives - CV_false_positives,
    same_sensitivity = same_sensitivity,
    same_direct_decisions = same_direct_decisions,
    improved_specificity = improved_specificity,
    desired_pattern = gamma > 0 && same_sensitivity && improved_specificity,
    all_direct_found = OLS_true_positives == length(direct_causes) &&
      CV_true_positives == length(direct_causes),
    OLS_rejected_direct = collapse_variables(
      direct_table$variable[direct_table$OLS_reject]
    ),
    CV_rejected_direct = collapse_variables(
      direct_table$variable[direct_table$CV_reject]
    ),
    OLS_rejected_noncauses = collapse_variables(
      non_cause_table$variable[non_cause_table$OLS_reject]
    ),
    CV_rejected_noncauses = collapse_variables(
      non_cause_table$variable[non_cause_table$CV_reject]
    ),
    OLS_only_direct = collapse_variables(
      direct_table$variable[direct_table$OLS_reject & !direct_table$CV_reject]
    ),
    CV_only_direct = collapse_variables(
      direct_table$variable[!direct_table$OLS_reject & direct_table$CV_reject]
    ),
    OLS_only_noncauses = collapse_variables(
      non_cause_table$variable[
        non_cause_table$OLS_reject & !non_cause_table$CV_reject
      ]
    ),
    CV_only_noncauses = collapse_variables(
      non_cause_table$variable[
        !non_cause_table$OLS_reject & non_cause_table$CV_reject
      ]
    ),
    status = "ok",
    error = "",
    stringsAsFactors = FALSE
  )
}

screen_pair <- function(pair_row, retain_table = FALSE) {
  pair_id <- as.integer(pair_row[[1]])
  file_index_1 <- as.integer(pair_row[[2]])
  file_index_2 <- as.integer(pair_row[[3]])
  file_1 <- environment_files[file_index_1]
  file_2 <- environment_files[file_index_2]

  tryCatch({
    model_data <- prepare_pair(
      environment_data[[file_1]], environment_data[[file_2]]
    )
    cv <- fit_basis_cross_validation(model_data)
    fitted <- fit_estimators_and_tests(model_data, cv$gamma)
    summary <- summarize_tests(
      pair_id, file_1, file_2, cv$gamma, cv$validation_risk,
      fitted$table
    )
    summary$n_1 <- nrow(model_data$Xe)
    summary$n_2 <- nrow(model_data$Xo)
    if (retain_table) {
      list(summary = summary, table = fitted$table)
    } else {
      summary
    }
  }, error = function(condition) {
    data.frame(
      pair_id = pair_id,
      environment_1 = file_1,
      environment_2 = file_2,
      n_1 = NA_integer_, n_2 = NA_integer_, gamma_CV = NA_real_,
      CV_validation_Rdelta_plus = NA_real_,
      OLS_direct_rejected = NA_integer_, CV_direct_rejected = NA_integer_,
      sensitivity_OLS = NA_real_, sensitivity_CV = NA_real_,
      OLS_noncauses_rejected = NA_integer_,
      CV_noncauses_rejected = NA_integer_,
      specificity_OLS = NA_real_, specificity_CV = NA_real_,
      precision_OLS = NA_real_, precision_CV = NA_real_,
      F1_OLS = NA_real_, F1_CV = NA_real_,
      balanced_accuracy_OLS = NA_real_,
      balanced_accuracy_CV = NA_real_,
      MCC_OLS = NA_real_, MCC_CV = NA_real_,
      specificity_gain_count = NA_integer_, same_sensitivity = FALSE,
      same_direct_decisions = FALSE, improved_specificity = FALSE,
      desired_pattern = FALSE, all_direct_found = FALSE,
      OLS_rejected_direct = "", CV_rejected_direct = "",
      OLS_rejected_noncauses = "", CV_rejected_noncauses = "",
      OLS_only_direct = "", CV_only_direct = "",
      OLS_only_noncauses = "", CV_only_noncauses = "",
      status = "error", error = conditionMessage(condition),
      stringsAsFactors = FALSE
    )
  })
}
