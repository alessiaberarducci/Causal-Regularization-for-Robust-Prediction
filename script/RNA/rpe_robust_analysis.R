# RPE robust analysis

# Load estimation, validation and application functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/rpe_robust_functions.R")

set.seed(123)

# Model variables.
environment_column <- "interventions"
observed_genes <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000147604",
  "ENSG00000110700", "ENSG00000172757", "ENSG00000133112",
  "ENSG00000067225", "ENSG00000108518", "ENSG00000125691",
  "ENSG00000173812"
)
target_gene <- "ENSG00000173812"
predictor_genes <- setdiff(observed_genes, target_gene)
control_environment <- "non-targeting"

# Load environments with at least 100 cells.
dat <- read.csv("dataset_rpe.csv", check.names = FALSE)
dat <- filter_environments(dat, environment_column, minimum_cells = 100)
available_environments <- unique(as.character(dat[[environment_column]]))

# Divide interventions on the nine predictors into balanced groups of 4 and 5.
intervention_sizes <- table(as.character(dat[[environment_column]]))[
  predictor_genes
]
candidate_splits <- combn(predictor_genes, 4, simplify = FALSE)
cell_imbalance <- vapply(candidate_splits, function(first_group) {
  second_group <- setdiff(predictor_genes, first_group)
  abs(sum(intervention_sizes[first_group]) -
        sum(intervention_sizes[second_group]))
}, numeric(1))
first_training <- candidate_splits[[which.min(cell_imbalance)]]
second_training <- setdiff(predictor_genes, first_training)
training_interventions <- c(first_training, second_training)

training <- make_two_environment_data(
  dat, environment_column, first_training, second_training,
  predictor_genes, target_gene
)
raw_data_train <- training$data
centering <- compute_common_centering(raw_data_train)
data_train <- centering$data

# All remaining interventions form the candidate test set.
test_environments <- setdiff(
  available_environments,
  c(training_interventions, target_gene, control_environment, "excluded")
)

# OLS and the limiting RL solution do not depend on the fold count.
training_moments <- moments(data_train)
beta_ols <- setNames(drop(compute_ols(training_moments)), predictor_genes)
beta_rl <- setNames(drop(compute_cd(training_moments)), predictor_genes)
intercept_ols <- recover_intercept(beta_ols, centering)
intercept_rl <- recover_intercept(beta_rl, centering)

# Construct the augmented matrices defining C_gamma.
training_xy_1 <- cbind(data_train$Xe, data_train$ye)
training_xy_2 <- cbind(data_train$Xo, data_train$yo)
M_1 <- crossprod(training_xy_1) / nrow(training_xy_1)
M_2 <- crossprod(training_xy_2) / nrow(training_xy_2)
M_plus <- M_1 + M_2
G_delta_plus <- build_matrices(data_train)$Gdelta_plus
G_beta_rl <- drop(G_delta_plus %*% beta_rl)
M_delta_plus <- rbind(
  cbind(G_delta_plus, G_beta_rl),
  c(G_beta_rl, drop(crossprod(beta_rl, G_beta_rl)))
)

compute_environment_moment <- function(environment) {
  current <- dat[dat[[environment_column]] == environment, , drop = FALSE]
  X <- sweep(
    as.matrix(current[, predictor_genes, drop = FALSE]),
    2, centering$x_center, "-"
  )
  y <- as.numeric(current[[target_gene]]) - centering$y_center
  complete <- complete.cases(X, y)
  XY <- cbind(X[complete, , drop = FALSE], y[complete])
  crossprod(XY) / nrow(XY)
}
environment_moments <- setNames(
  lapply(test_environments, compute_environment_moment),
  test_environments
)

# Fit CV and determine C_gamma membership for each requested fold count.
gamma_grid <- seq(0, 40, by = 0.01)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
fold_values <- c(5, 10, 15, 20)

fit_for_folds <- function(folds) {
  # Reset the seed so fold-count comparisons use the same randomization rule.
  set.seed(123)
  sample_sizes <- validate_cv_data(raw_data_train, folds)
  fold_assignment <- list(
    environment_e = split_idx(folds, sample_sizes[["environment_e"]]),
    environment_o = split_idx(folds, sample_sizes[["environment_o"]])
  )
  cv_result <- centered_cross_validation(
    raw_data_train, estimators, fold_assignment
  )
  gamma_cv <- as.numeric(cv_result$best_gamma)
  beta_cv <- setNames(
    drop(make_estimator(gamma_cv)(data_train)), predictor_genes
  )
  intercept_cv <- recover_intercept(beta_cv, centering)

  all_test_risks <- compute_test_risks(
    dat, test_environments, environment_column,
    beta_ols, beta_cv, predictor_genes, target_gene,
    intercept_ols = intercept_ols, intercept_cv = intercept_cv
  )

  M_gamma <- 0.5 * M_plus + 0.5 * gamma_cv * M_delta_plus
  minimum_eigenvalues <- vapply(environment_moments, function(M_environment) {
    min(eigen(M_gamma - M_environment, symmetric = TRUE)$values)
  }, numeric(1))
  membership <- data.frame(
    environment = names(minimum_eigenvalues),
    minimum_eigenvalue = minimum_eigenvalues,
    in_Cgamma = minimum_eigenvalues >= 0,
    row.names = NULL
  )
  selected_risks <- merge(
    membership[membership$in_Cgamma, , drop = FALSE],
    all_test_risks, by = "environment", all = FALSE, sort = FALSE
  )

  training_risks <- compare_training_risks(
    data_train,
    list(OLS = beta_ols, CV = beta_cv, RL = beta_rl),
    gamma_cv
  )
  # Worst risk: R_gamma = R_plus / 2 + gamma R_delta^+ / 2.
  names(training_risks)[names(training_risks) == "Rgamma"] <-
    "worst_risk_Rgamma"

  list(
    folds = folds,
    gamma_cv = gamma_cv,
    beta_cv = beta_cv,
    intercept_cv = intercept_cv,
    all_test_risks = all_test_risks,
    membership = membership,
    selected_risks = selected_risks,
    training_risks = training_risks,
    fold_assignment = fold_assignment,
    cv_result = cv_result
  )
}

fold_results <- setNames(lapply(fold_values, fit_for_folds), fold_values)
main_result <- fold_results[["10"]]

# Sensitivity analysis: no centering and no intercept.
set.seed(123)
cv_no_centering <- cross_validation(10, raw_data_train, estimators)
gamma_no_centering <- as.numeric(cv_no_centering$best_gamma)
moments_no_centering <- moments(raw_data_train)
beta_ols_no_centering <- setNames(
  drop(compute_ols(moments_no_centering)), predictor_genes
)
beta_cv_no_centering <- setNames(
  drop(make_estimator(gamma_no_centering)(raw_data_train)), predictor_genes
)
beta_rl_no_centering <- setNames(
  drop(compute_cd(moments_no_centering)), predictor_genes
)
all_risks_no_centering <- compute_test_risks(
  dat, test_environments, environment_column,
  beta_ols_no_centering, beta_cv_no_centering,
  predictor_genes, target_gene,
  intercept_ols = 0, intercept_cv = 0
)

# Recompute C_gamma using raw, uncentred augmented moments.
training_xy_1_no_centering <- cbind(raw_data_train$Xe, raw_data_train$ye)
training_xy_2_no_centering <- cbind(raw_data_train$Xo, raw_data_train$yo)
M_plus_no_centering <-
  crossprod(training_xy_1_no_centering) / nrow(training_xy_1_no_centering) +
  crossprod(training_xy_2_no_centering) / nrow(training_xy_2_no_centering)
G_delta_plus_no_centering <- build_matrices(raw_data_train)$Gdelta_plus
G_beta_rl_no_centering <- drop(
  G_delta_plus_no_centering %*% beta_rl_no_centering
)
M_delta_plus_no_centering <- rbind(
  cbind(G_delta_plus_no_centering, G_beta_rl_no_centering),
  c(
    G_beta_rl_no_centering,
    drop(crossprod(beta_rl_no_centering, G_beta_rl_no_centering))
  )
)
M_gamma_no_centering <- 0.5 * M_plus_no_centering +
  0.5 * gamma_no_centering * M_delta_plus_no_centering
minimum_eigenvalues_no_centering <- vapply(
  test_environments,
  function(environment) {
    current <- dat[dat[[environment_column]] == environment, , drop = FALSE]
    X <- as.matrix(current[, predictor_genes, drop = FALSE])
    y <- as.numeric(current[[target_gene]])
    complete <- complete.cases(X, y)
    XY <- cbind(X[complete, , drop = FALSE], y[complete])
    M_environment <- crossprod(XY) / nrow(XY)
    min(eigen(
      M_gamma_no_centering - M_environment, symmetric = TRUE
    )$values)
  },
  numeric(1)
)
membership_no_centering <- data.frame(
  environment = test_environments,
  minimum_eigenvalue = minimum_eigenvalues_no_centering,
  in_Cgamma = minimum_eigenvalues_no_centering >= 0
)
selected_risks_no_centering <- merge(
  membership_no_centering[membership_no_centering$in_Cgamma, ],
  all_risks_no_centering,
  by = "environment", all = FALSE, sort = FALSE
)
worst_risks_no_centering <- compare_training_risks(
  raw_data_train,
  list(
    OLS = beta_ols_no_centering,
    CV = beta_cv_no_centering,
    RL = beta_rl_no_centering
  ),
  gamma_no_centering
)
names(worst_risks_no_centering)[
  names(worst_risks_no_centering) == "Rgamma"
] <- "worst_risk_Rgamma"

# Numerical summaries.
summarise_risks <- function(risks, folds, gamma_cv, subset_name) {
  data.frame(
    subset = subset_name,
    folds = folds,
    gamma_cv = gamma_cv,
    environments = nrow(risks),
    mean_OLS = mean(risks$OLS),
    mean_CV = mean(risks$CV),
    median_OLS = median(risks$OLS),
    median_CV = median(risks$CV),
    maximum_observed_MSE_OLS = max(risks$OLS),
    maximum_observed_MSE_CV = max(risks$CV),
    proportion_CV_better = mean(risks$CV < risks$OLS)
  )
}

all_test_summary <- summarise_risks(
  main_result$all_test_risks, 10, main_result$gamma_cv,
  "all_test_environments"
)
cgamma_summary <- do.call(rbind, lapply(fold_results, function(result) {
  summarise_risks(
    result$selected_risks, result$folds, result$gamma_cv,
    "Cgamma"
  )
}))
centering_comparison_summary <- rbind(
  transform(
    summarise_risks(
      main_result$all_test_risks, 10, main_result$gamma_cv, "all_test"
    ),
    specification = "common_centering_with_intercept"
  ),
  transform(
    summarise_risks(
      main_result$selected_risks, 10, main_result$gamma_cv, "Cgamma"
    ),
    specification = "common_centering_with_intercept"
  ),
  transform(
    summarise_risks(
      all_risks_no_centering, 10, gamma_no_centering, "all_test"
    ),
    specification = "no_centering_no_intercept"
  ),
  transform(
    summarise_risks(
      selected_risks_no_centering, 10, gamma_no_centering, "Cgamma"
    ),
    specification = "no_centering_no_intercept"
  )
)
centering_comparison_summary <- centering_comparison_summary[, c(
  "specification", setdiff(
    names(centering_comparison_summary), "specification"
  )
)]
coefficient_table <- data.frame(
  gene = predictor_genes,
  OLS = beta_ols,
  CV = main_result$beta_cv,
  RL = beta_rl,
  row.names = NULL
)
training_setup <- data.frame(
  training_environment = c("environment_1", "environment_2"),
  interventions = c(
    paste(first_training, collapse = "; "),
    paste(second_training, collapse = "; ")
  ),
  cells = c(nrow(data_train$Xe), nrow(data_train$Xo))
)

# Relate test performance to second-moment directions within C_gamma.
absolute_symmetric_part <- function(matrix) {
  decomposition <- eigen(matrix, symmetric = TRUE)
  decomposition$vectors %*%
    diag(abs(decomposition$values), nrow = length(decomposition$values)) %*%
    t(decomposition$vectors)
}
reference_second_moment <- 0.5 * (
  training_moments$XXe + training_moments$XXo
)
Gdelta_frobenius <- sqrt(sum(G_delta_plus^2))
cgamma_direction_diagnostics <- do.call(rbind, lapply(
  main_result$selected_risks$environment,
  function(environment) {
    current <- dat[dat[[environment_column]] == environment, , drop = FALSE]
    X_raw <- as.matrix(current[, predictor_genes, drop = FALSE])
    y_raw <- as.numeric(current[[target_gene]])
    complete <- complete.cases(X_raw, y_raw)
    X_raw <- X_raw[complete, , drop = FALSE]
    y_raw <- y_raw[complete]
    X <- sweep(X_raw, 2, centering$x_center, "-")
    test_second_moment <- crossprod(X) / nrow(X)
    test_shift_plus <- absolute_symmetric_part(
      test_second_moment - reference_second_moment
    )
    denominator <- Gdelta_frobenius * sqrt(sum(test_shift_plus^2))
    alignment <- if (denominator > 0) {
      sum(G_delta_plus * test_shift_plus) / denominator
    } else {
      NA_real_
    }
    risk_row <- main_result$selected_risks[
      main_result$selected_risks$environment == environment, , drop = FALSE
    ]
    prediction_ols <- intercept_ols + drop(X_raw %*% beta_ols)
    prediction_cv <- main_result$intercept_cv +
      drop(X_raw %*% main_result$beta_cv)
    residual_ols <- y_raw - prediction_ols
    prediction_change <- prediction_cv - prediction_ols
    quadratic_cost <- mean(prediction_change^2)
    residual_correction <- 2 * mean(residual_ols * prediction_change)
    data.frame(
      environment = environment,
      second_moment_alignment = alignment,
      MSE_OLS = risk_row$OLS,
      MSE_CV = risk_row$CV,
      MSE_CV_minus_OLS = risk_row$CV - risk_row$OLS,
      relative_MSE_difference = (risk_row$CV - risk_row$OLS) / risk_row$OLS,
      quadratic_cost = quadratic_cost,
      residual_correction = residual_correction,
      reconstructed_MSE_difference = quadratic_cost - residual_correction,
      CV_better = risk_row$CV < risk_row$OLS
    )
  }
))
second_moment_association <- cor.test(
  cgamma_direction_diagnostics$second_moment_alignment,
  cgamma_direction_diagnostics$MSE_CV_minus_OLS,
  method = "spearman", exact = FALSE
)
cgamma_direction_summary <- data.frame(
  environments = nrow(cgamma_direction_diagnostics),
  environments_CV_better = sum(cgamma_direction_diagnostics$CV_better),
  proportion_CV_better = mean(cgamma_direction_diagnostics$CV_better),
  mean_relative_MSE_difference = mean(
    cgamma_direction_diagnostics$relative_MSE_difference
  ),
  median_relative_MSE_difference = median(
    cgamma_direction_diagnostics$relative_MSE_difference
  ),
  spearman_rho = unname(second_moment_association$estimate),
  p_value = second_moment_association$p.value,
  maximum_decomposition_error = max(abs(
    cgamma_direction_diagnostics$reconstructed_MSE_difference -
      cgamma_direction_diagnostics$MSE_CV_minus_OLS
  ))
)

# Save the compact set of tables.
rna_output_directory <- "visualization/rna_output"
dir.create(rna_output_directory, recursive = TRUE, showWarnings = FALSE)
rna_output_file <- function(filename) file.path(rna_output_directory, filename)
write.csv(training_setup, rna_output_file("rpe_training_setup.csv"), row.names = FALSE)
write.csv(coefficient_table, rna_output_file("rpe_coefficients_10fold.csv"), row.names = FALSE)
write.csv(all_test_summary, rna_output_file("rpe_all_test_summary_10fold.csv"), row.names = FALSE)
write.csv(
  main_result$all_test_risks,
  rna_output_file("rpe_all_test_risks_10fold.csv"), row.names = FALSE
)
write.csv(cgamma_summary, rna_output_file("rpe_Cgamma_fold_summary.csv"), row.names = FALSE)
write.csv(
  main_result$training_risks,
  rna_output_file("rpe_worst_case_risk_10fold.csv"), row.names = FALSE
)
write.csv(
  main_result$membership,
  rna_output_file("rpe_Cgamma_membership_10fold.csv"), row.names = FALSE
)
write.csv(
  main_result$selected_risks,
  rna_output_file("rpe_Cgamma_test_risks_10fold.csv"), row.names = FALSE
)
write.csv(
  cgamma_direction_diagnostics,
  rna_output_file("rpe_Cgamma_second_moment_diagnostics_10fold.csv"),
  row.names = FALSE
)
write.csv(
  cgamma_direction_summary,
  rna_output_file("rpe_Cgamma_second_moment_summary_10fold.csv"),
  row.names = FALSE
)
write.csv(
  centering_comparison_summary,
  rna_output_file("rpe_centering_intercept_comparison_10fold.csv"),
  row.names = FALSE
)
write.csv(
  worst_risks_no_centering,
  rna_output_file("rpe_no_centering_no_intercept_worst_risk_10fold.csv"),
  row.names = FALSE
)

# Common plotting function: every observed environment remains visible.
plot_risk_boxplot <- function(risks, file = NULL, title = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 6, height = 4.5, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 6, height = 4.5)
  }
  graphics::layout(1)
  old_par <- graphics::par(mfrow = c(1, 1), mar = c(4.2, 4.6, 2, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  values <- risks[, c("OLS", "CV")]
  limits <- range(unlist(values), finite = TRUE)
  padding <- max(0.02, 0.05 * diff(limits))
  graphics::boxplot(
    values, outline = FALSE,
    ylim = c(max(0, limits[1] - padding), limits[2] + padding),
    xlab = "", ylab = expression(R[e](hat(beta))), main = title
  )
  graphics::stripchart(
    values, vertical = TRUE, method = "jitter", add = TRUE,
    pch = 16, cex = 0.5,
    col = grDevices::adjustcolor("black", alpha.f = 0.35)
  )
}

# Main 10-fold comparisons.
plot_risk_boxplot(
  main_result$all_test_risks,
  rna_output_file("rpe_all_test_OLS_CV_boxplot_10fold.pdf")
)
plot_risk_boxplot(
  main_result$selected_risks,
  rna_output_file("rpe_Cgamma_OLS_CV_boxplot_10fold.pdf")
)

# C_gamma sensitivity to 5, 10, 15 and 20 folds.
plot_cgamma_fold_comparison <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 10, height = 7, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 10, height = 7)
  }
  graphics::layout(matrix(1:4, nrow = 2, byrow = TRUE))
  old_par <- graphics::par(mar = c(3.8, 4.2, 2.4, 0.8))
  on.exit(graphics::par(old_par), add = TRUE)

  common_limits <- range(unlist(lapply(fold_results, function(result) {
    unlist(result$selected_risks[, c("OLS", "CV")])
  })), finite = TRUE)
  padding <- max(0.02, 0.05 * diff(common_limits))
  common_limits <- c(
    max(0, common_limits[1] - padding), common_limits[2] + padding
  )

  for (result in fold_results) {
    values <- result$selected_risks[, c("OLS", "CV")]
    graphics::boxplot(
      values, outline = FALSE, ylim = common_limits,
      xlab = "", ylab = expression(R[e](hat(beta))),
      main = paste0(result$folds, " folds")
    )
    graphics::stripchart(
      values, vertical = TRUE, method = "jitter", add = TRUE,
      pch = 16, cex = 0.4,
      col = grDevices::adjustcolor("black", alpha.f = 0.3)
    )
  }
}

plot_cgamma_fold_comparison(
  rna_output_file("rpe_Cgamma_fold_comparison_boxplots.pdf")
)

# Show when CV improves within C_gamma as a function of shift direction.
plot_second_moment_diagnostics <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 6.5, height = 4.8, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 6.5, height = 4.8)
  }
  graphics::layout(1)
  old_par <- graphics::par(mfrow = c(1, 1), mar = c(4.6, 5, 1.2, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  x <- cgamma_direction_diagnostics$second_moment_alignment
  y <- cgamma_direction_diagnostics$MSE_CV_minus_OLS
  point_colours <- ifelse(
    cgamma_direction_diagnostics$CV_better, "#0072B2", "grey55"
  )
  graphics::plot(
    x, y, pch = 16,
    col = grDevices::adjustcolor(point_colours, alpha.f = 0.7),
    xlab = expression(A[e]^moment),
    ylab = expression(MSE[e](hat(beta)[CV]) - MSE[e](hat(beta)[OLS]))
  )
  graphics::abline(h = 0, lty = 2, col = "grey40")
  smooth <- stats::lowess(x, y, f = 2 / 3)
  graphics::lines(smooth, lwd = 2, col = "#0072B2")
  graphics::legend(
    "topright",
    legend = sprintf(
      "CV better: %d/%d (%.1f%%)",
      sum(cgamma_direction_diagnostics$CV_better),
      nrow(cgamma_direction_diagnostics),
      100 * mean(cgamma_direction_diagnostics$CV_better)
    ),
    bty = "n", cex = 0.85
  )
}

plot_second_moment_diagnostics(
  rna_output_file("rpe_Cgamma_second_moment_diagnostics_10fold.pdf")
)

# Paired MSE comparison: CV improves below the diagonal.
plot_paired_mse <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 5.8, height = 5.2, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 5.8, height = 5.2)
  }
  graphics::layout(1)
  old_par <- graphics::par(mfrow = c(1, 1), mar = c(4.8, 5, 1.2, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  limits <- range(
    cgamma_direction_diagnostics$MSE_OLS,
    cgamma_direction_diagnostics$MSE_CV,
    finite = TRUE
  )
  padding <- 0.05 * diff(limits)
  limits <- limits + c(-padding, padding)
  point_colours <- ifelse(
    cgamma_direction_diagnostics$CV_better, "#0072B2", "grey55"
  )
  graphics::plot(
    cgamma_direction_diagnostics$MSE_OLS,
    cgamma_direction_diagnostics$MSE_CV,
    pch = 16,
    col = grDevices::adjustcolor(point_colours, alpha.f = 0.7),
    xlim = limits, ylim = limits,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CV]))
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
  graphics::legend(
    "topleft",
    legend = sprintf(
      "CV better: %d/%d (%.1f%%)",
      sum(cgamma_direction_diagnostics$CV_better),
      nrow(cgamma_direction_diagnostics),
      100 * mean(cgamma_direction_diagnostics$CV_better)
    ),
    bty = "n", cex = 0.85
  )
}

plot_paired_mse(rna_output_file("rpe_Cgamma_paired_MSE_10fold.pdf"))

# Squared-loss decomposition: CV improves above the diagonal.
plot_squared_loss_decomposition <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 5.8, height = 5.2, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 5.8, height = 5.2)
  }
  graphics::layout(1)
  old_par <- graphics::par(mfrow = c(1, 1), mar = c(4.8, 5, 1.2, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  x_limits <- range(cgamma_direction_diagnostics$quadratic_cost, finite = TRUE)
  y_limits <- range(cgamma_direction_diagnostics$residual_correction, finite = TRUE)
  x_padding <- max(0.0005, 0.05 * diff(x_limits))
  y_padding <- max(0.0005, 0.05 * diff(y_limits))
  x_limits <- c(max(0, x_limits[1] - x_padding), x_limits[2] + x_padding)
  y_limits <- y_limits + c(-y_padding, y_padding)
  point_colours <- ifelse(
    cgamma_direction_diagnostics$CV_better, "#0072B2", "grey55"
  )
  graphics::plot(
    cgamma_direction_diagnostics$quadratic_cost,
    cgamma_direction_diagnostics$residual_correction,
    pch = 16,
    col = grDevices::adjustcolor(point_colours, alpha.f = 0.7),
    xlim = x_limits, ylim = y_limits,
    xlab = "Quadratic cost",
    ylab = "OLS-residual correction"
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
  graphics::legend(
    "topright", legend = "CV improves above the diagonal",
    bty = "n", cex = 0.8
  )
}

plot_squared_loss_decomposition(
  rna_output_file("rpe_Cgamma_squared_loss_decomposition_10fold.pdf")
)

# Compare the two specifications on all tests and their own C_gamma subsets.
plot_centering_comparison <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 10, height = 7, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 10, height = 7)
  }
  graphics::layout(matrix(1:4, nrow = 2, byrow = TRUE))
  old_par <- graphics::par(mar = c(3.8, 4.2, 2.5, 0.8))
  on.exit(graphics::par(old_par), add = TRUE)
  panels <- list(
    "Centred: all tests" = main_result$all_test_risks,
    "No intercept: all tests" = all_risks_no_centering,
    "Centred: Cgamma" = main_result$selected_risks,
    "No intercept: Cgamma" = selected_risks_no_centering
  )
  for (panel_name in names(panels)) {
    values <- panels[[panel_name]][, c("OLS", "CV")]
    limits <- range(unlist(values), finite = TRUE)
    padding <- max(0.02, 0.05 * diff(limits))
    limits <- c(max(0, limits[1] - padding), limits[2] + padding)
    graphics::boxplot(
      values, outline = FALSE, ylim = limits,
      xlab = "", ylab = expression(R[e](hat(beta))), main = panel_name
    )
    graphics::stripchart(
      values, vertical = TRUE, method = "jitter", add = TRUE,
      pch = 16, cex = 0.4,
      col = grDevices::adjustcolor("black", alpha.f = 0.3)
    )
  }
}

plot_centering_comparison(
  rna_output_file("rpe_centering_intercept_comparison_10fold.pdf")
)

# Display figures when the script is run interactively.
if (interactive()) {
  plot_risk_boxplot(
    main_result$all_test_risks,
    title = "All test environments: 10 folds"
  )
  plot_risk_boxplot(
    main_result$selected_risks,
    title = expression(C[gamma[CV]] * ": 10 folds")
  )
  plot_cgamma_fold_comparison()
  plot_second_moment_diagnostics()
  plot_paired_mse()
  plot_squared_loss_decomposition()
  plot_centering_comparison()
}

cat("Training environment 1:\n")
print(first_training)
cat("Training environment 2:\n")
print(second_training)
print(training_setup, row.names = FALSE)
cat("\nAll test environments, 10-fold CV:\n")
print(all_test_summary, row.names = FALSE)
cat("\nC_gamma-selected environments:\n")
print(cgamma_summary, row.names = FALSE)
cat("\nEstimated worst-case risk, 10-fold CV:\n")
print(main_result$training_risks, row.names = FALSE)
cat("\nSecond-moment diagnostic within C_gamma:\n")
print(cgamma_direction_summary, row.names = FALSE)
cat("\nSensitivity to removing centering and the intercept:\n")
print(centering_comparison_summary, row.names = FALSE)
cat("\nWorst risk without centering or intercept:\n")
print(worst_risks_no_centering, row.names = FALSE)
