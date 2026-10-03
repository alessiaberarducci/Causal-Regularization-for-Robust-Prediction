# RNA prediction across training sample sizes.

source("code/functions/datautil.R")
source("code/functions/cd.R")
source("code/functions/ols.R")
source("code/functions/measurements.R")
source("code/functions/cross_valid_fun_n.R")
source("code/functions/rpe_robust_functions.R")

base_seed <- 123L
sample_sizes <- as.integer(strsplit(
  Sys.getenv("RPE_SAMPLE_SIZES", "650,500,250,150,100"), ",",
  fixed = TRUE
)[[1]])
if (any(!is.finite(sample_sizes)) || any(sample_sizes < 1L)) {
  stop("RPE_SAMPLE_SIZES must be a comma-separated list of positive integers.")
}
folds <- 10L
output_suffix <- Sys.getenv("RPE_OUTPUT_SUFFIX", "")
membership_tolerance <- as.numeric(
  Sys.getenv("RPE_MEMBERSHIP_TOLERANCE", "0.05")
)
if (!is.finite(membership_tolerance) || membership_tolerance < 0) {
  stop("RPE_MEMBERSHIP_TOLERANCE must be a nonnegative finite number.")
}

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

dat <- read.csv("data/dataset_rpe.csv", check.names = FALSE)
dat <- filter_environments(dat, environment_column, minimum_cells = 100)
available_environments <- unique(as.character(dat[[environment_column]]))

all_environment_sizes <- table(as.character(dat[[environment_column]]))
read_training_group <- function(variable) {
  value <- Sys.getenv(variable, "")
  if (!nzchar(value)) {
    stop(variable, " must be provided explicitly by the calling script.")
  }
  group <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(group) == 0L || any(!nzchar(group))) {
    stop(variable, " must contain a comma-separated list of environments.")
  }
  group
}
first_training <- read_training_group("RPE_FIRST_TRAINING")
second_training <- read_training_group("RPE_SECOND_TRAINING")
training_imbalance <- abs(
  sum(all_environment_sizes[first_training]) -
    sum(all_environment_sizes[second_training])
)
training_interventions <- c(first_training, second_training)

training <- make_two_environment_data(
  dat, environment_column, first_training, second_training,
  predictor_genes, target_gene
)
full_training <- training$data
if (max(sample_sizes) > min(nrow(full_training$Xe), nrow(full_training$Xo))) {
  stop("The requested sample size exceeds one of the training groups.")
}

test_environments <- setdiff(
  available_environments,
  c(training_interventions, target_gene, control_environment, "excluded")
)

set.seed(base_seed)
row_order_e <- sample(seq_len(nrow(full_training$Xe)))
row_order_o <- sample(seq_len(nrow(full_training$Xo)))

gamma_grid <- seq(0, 40, by = 0.01)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)

run_subsample_experiment <- function(sample_size = NULL, full_sample = FALSE) {
  if (full_sample) {
    rows_e <- seq_len(nrow(full_training$Xe))
    rows_o <- seq_len(nrow(full_training$Xo))
    sample_size_e <- length(rows_e)
    sample_size_o <- length(rows_o)
    sample_label <- sprintf(
      "Full sample: n1 = %d, n2 = %d", sample_size_e, sample_size_o
    )
    experiment_seed <- base_seed
  } else {
    rows_e <- row_order_e[seq_len(sample_size)]
    rows_o <- row_order_o[seq_len(sample_size)]
    sample_size_e <- sample_size
    sample_size_o <- sample_size
    sample_label <- sprintf("n = %d", sample_size)
    experiment_seed <- base_seed + sample_size
  }
  raw_training <- subset_data(
    full_training, rows_e = rows_e, rows_o = rows_o
  )
  centering <- compute_common_centering(raw_training)
  centered_training <- centering$data

  set.seed(experiment_seed)
  cv_sizes <- validate_cv_data(raw_training, folds)
  fold_assignment <- list(
    environment_e = split_idx(folds, cv_sizes[["environment_e"]]),
    environment_o = split_idx(folds, cv_sizes[["environment_o"]])
  )
  cv_result <- centered_cross_validation(
    raw_training, estimators, fold_assignment
  )
  gamma_cv <- as.numeric(cv_result$best_gamma)

  training_moments <- moments(centered_training)
  beta_ols <- setNames(
    drop(compute_ols(training_moments)), predictor_genes
  )
  beta_cv <- setNames(
    drop(make_estimator(gamma_cv)(centered_training)), predictor_genes
  )
  beta_rl <- setNames(
    drop(compute_cd(training_moments)), predictor_genes
  )
  intercept_ols <- recover_intercept(beta_ols, centering)
  intercept_cv <- recover_intercept(beta_cv, centering)

  worst_risk <- compare_training_risks(
    centered_training,
    list(OLS = beta_ols, CV = beta_cv),
    gamma_cv
  )
  names(worst_risk)[names(worst_risk) == "Rgamma"] <-
    "worst_risk_Rgamma"
  worst_risk$sample_size <- if (full_sample) NA_integer_ else sample_size
  worst_risk$sample_size_e <- sample_size_e
  worst_risk$sample_size_o <- sample_size_o
  worst_risk$sample_label <- sample_label
  worst_risk$gamma_cv <- gamma_cv

  risks <- compute_test_risks(
    dat, test_environments, environment_column,
    beta_ols, beta_cv, predictor_genes, target_gene,
    intercept_ols = intercept_ols, intercept_cv = intercept_cv
  )

  training_xy_1 <- cbind(centered_training$Xe, centered_training$ye)
  training_xy_2 <- cbind(centered_training$Xo, centered_training$yo)
  M_1 <- crossprod(training_xy_1) / nrow(training_xy_1)
  M_2 <- crossprod(training_xy_2) / nrow(training_xy_2)
  M_0 <- 0.5 * (M_1 + M_2)

  G_delta_plus <- build_matrices(centered_training)$Gdelta_plus
  G_beta_rl <- drop(G_delta_plus %*% beta_rl)
  M_delta_plus <- rbind(
    cbind(G_delta_plus, G_beta_rl),
    c(G_beta_rl, drop(crossprod(beta_rl, G_beta_rl)))
  )
  M_gamma <- M_0 + 0.5 * gamma_cv * M_delta_plus

  compute_test_moment <- function(environment) {
    current <- dat[
      as.character(dat[[environment_column]]) == environment,
      , drop = FALSE
    ]
    X <- sweep(
      as.matrix(current[, predictor_genes, drop = FALSE]),
      2, centering$x_center, "-"
    )
    y <- as.numeric(current[[target_gene]]) - centering$y_center
    complete <- complete.cases(X, y)
    XY <- cbind(X[complete, , drop = FALSE], y[complete])
    crossprod(XY) / nrow(XY)
  }
  test_moments <- setNames(
    lapply(test_environments, compute_test_moment),
    test_environments
  )
  min_eigen_Cgamma <- vapply(test_moments, function(M_test) {
    min(eigen(
      M_gamma - M_test, symmetric = TRUE, only.values = TRUE
    )$values)
  }, numeric(1))
  min_eigen_C0 <- vapply(test_moments, function(M_test) {
    min(eigen(
      M_0 - M_test, symmetric = TRUE, only.values = TRUE
    )$values)
  }, numeric(1))

  membership <- data.frame(
    environment = names(test_moments),
    minimum_eigenvalue_Cgamma = min_eigen_Cgamma,
    minimum_eigenvalue_C0 = min_eigen_C0,
    in_Cgamma = min_eigen_Cgamma >= -membership_tolerance,
    in_C0 = min_eigen_C0 >= -membership_tolerance,
    row.names = NULL
  )
  result <- merge(risks, membership, by = "environment", sort = FALSE)
  result$sample_size <- if (full_sample) NA_integer_ else sample_size
  result$sample_size_e <- sample_size_e
  result$sample_size_o <- sample_size_o
  result$sample_label <- sample_label
  result$gamma_cv <- gamma_cv
  result$membership_tolerance <- membership_tolerance
  result$class <- ifelse(
    result$in_C0, "C0",
    ifelse(result$in_Cgamma, "Cgamma", "outside")
  )
  result$class <- factor(
    result$class, levels = c("outside", "Cgamma", "C0")
  )

  list(
    sample_size = if (full_sample) NA_integer_ else sample_size,
    sample_size_e = sample_size_e,
    sample_size_o = sample_size_o,
    sample_label = sample_label,
    experiment_seed = experiment_seed,
    full_sample = full_sample,
    gamma_cv = gamma_cv,
    worst_risk = worst_risk,
    data = result
  )
}

experiment_results <- c(
  list(full = run_subsample_experiment(full_sample = TRUE)),
  setNames(lapply(sample_sizes, run_subsample_experiment), sample_sizes)
)
plot_columns <- 2L
plot_rows <- ceiling(length(experiment_results) / plot_columns)
panel_ids <- c(
  seq_along(experiment_results),
  rep(0L, plot_rows * plot_columns - length(experiment_results))
)
panel_layout <- matrix(panel_ids, ncol = plot_columns, byrow = TRUE)
combined_results <- do.call(
  rbind, lapply(experiment_results, function(result) result$data)
)
row.names(combined_results) <- NULL
combined_worst_risks <- do.call(
  rbind, lapply(experiment_results, function(result) result$worst_risk)
)
row.names(combined_worst_risks) <- NULL

rna_output_directory <- "output/rna"
dir.create(rna_output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  combined_results,
  file.path(rna_output_directory, paste0(
    "rpe_subsample_Cgamma_C0_results",
    output_suffix, ".csv"
  )),
  row.names = FALSE
)
write.csv(
  combined_worst_risks,
  file.path(rna_output_directory, paste0(
    "rpe_subsample_worst_risk_comparison",
    output_suffix, ".csv"
  )),
  row.names = FALSE
)

scatter_colours <- c(
  outside = "#009E73",
  Cgamma = "#E69F00",
  C0 = "#0072B2"
)
boxplot_colours <- c(
  outside = "grey65",
  Cgamma = "#0072B2",
  C0 = "#D55E00"
)

scatter_file <- file.path(rna_output_directory, paste0(
  "rpe_subsample_OLS_CV_scatterplots",
  output_suffix, ".pdf"
))
grDevices::pdf(
  scatter_file, width = 10.5, height = 5.2 * plot_rows,
  family = "serif", useDingbats = FALSE
)
graphics::layout(panel_layout)
old_par <- graphics::par(
  family = "serif", mar = c(4.7, 4.7, 2.5, 0.8),
  lend = "round"
)
scatter_limits <- range(
  combined_results$OLS, combined_results$CV, finite = TRUE
)
scatter_padding <- max(0.01, 0.02 * diff(scatter_limits))
scatter_limits <- c(
  max(0, scatter_limits[1] - scatter_padding),
  scatter_limits[2] + scatter_padding
)

for (i in seq_along(experiment_results)) {
  result <- experiment_results[[i]]
  current <- result$data
  counts <- table(current$class)
  graphics::plot(
    current$OLS, current$CV,
    type = "n", xlim = scatter_limits, ylim = scatter_limits, asp = 1,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CV])),
    main = bquote(.(result$sample_label) ~ "," ~ gamma[CV] == .(
      result$gamma_cv
    ))
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey35")
  for (current_class in levels(current$class)) {
    selected <- current$class == current_class
    graphics::points(
      current$OLS[selected], current$CV[selected],
      pch = 16, cex = 0.58,
      col = grDevices::adjustcolor(
        scatter_colours[[current_class]], alpha.f = 0.72
      )
    )
  }
  graphics::legend(
    "topleft",
    legend = as.expression(list(
      bquote("Outside " * C[gamma] * ": " * .(counts[["outside"]])),
      bquote(C[gamma] - C[0] * ": " * .(counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(counts[["C0"]]))
    )),
    col = scatter_colours, pch = 16,
    pt.cex = 0.7, bty = "n", cex = 0.67
  )
}
graphics::par(old_par)
grDevices::dev.off()

boxplot_file <- file.path(rna_output_directory, paste0(
  "rpe_subsample_OLS_CV_boxplots",
  output_suffix, ".pdf"
))
grDevices::pdf(
  boxplot_file, width = 10.5, height = 5.5 * plot_rows,
  family = "serif", useDingbats = FALSE
)
graphics::layout(panel_layout)
old_par <- graphics::par(
  family = "serif", mar = c(4.4, 4.7, 3.3, 0.8),
  lend = "round"
)
boxplot_values <- c(combined_results$OLS, combined_results$CV)
boxplot_values <- boxplot_values[is.finite(boxplot_values)]
boxplot_upper_quantile <- 0.90
boxplot_upper <- unname(stats::quantile(
  boxplot_values, probs = boxplot_upper_quantile, names = FALSE
))
boxplot_limits <- c(0, 1.05 * boxplot_upper)

for (i in seq_along(experiment_results)) {
  result <- experiment_results[[i]]
  current <- result$data
  counts <- table(current$class)
  graphics::boxplot(
    current[, c("OLS", "CV")],
    names = c("OLS", "CV"), outline = FALSE,
    ylim = boxplot_limits,
    ylab = expression(MSE[e](hat(beta))),
    main = bquote(.(result$sample_label) ~ "," ~ gamma[CV] == .(
      result$gamma_cv
    ))
  )
  set.seed(result$experiment_seed)
  jitter_offset <- stats::runif(nrow(current), -0.12, 0.12)
  point_colours <- grDevices::adjustcolor(
    unname(boxplot_colours[as.character(current$class)]),
    alpha.f = 0.62
  )
  graphics::points(
    1 + jitter_offset, current$OLS,
    pch = 16, cex = 0.45, col = point_colours
  )
  graphics::points(
    2 + jitter_offset, current$CV,
    pch = 16, cex = 0.45, col = point_colours
  )
  clipped_ols <- sum(current$OLS > boxplot_limits[2], na.rm = TRUE)
  clipped_cv <- sum(current$CV > boxplot_limits[2], na.rm = TRUE)
  graphics::text(
    x = 0.58, y = 0.985 * boxplot_limits[2],
    labels = sprintf("Above limit: OLS %d, CV %d", clipped_ols, clipped_cv),
    adj = c(0, 1), cex = 0.64
  )
  graphics::legend(
    "topright",
    legend = as.expression(list(
      bquote("Outside " * C[gamma] * ": " * .(counts[["outside"]])),
      bquote(C[gamma] - C[0] * ": " * .(counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(counts[["C0"]]))
    )),
    col = boxplot_colours, pch = 16,
    pt.cex = 0.7, bty = "n", cex = 0.67
  )
}
graphics::par(old_par)
grDevices::dev.off()

summary_table <- do.call(rbind, lapply(experiment_results, function(result) {
  current <- result$data
  worst_risk_ols <- result$worst_risk$worst_risk_Rgamma[
    result$worst_risk$method == "OLS"
  ]
  worst_risk_cv <- result$worst_risk$worst_risk_Rgamma[
    result$worst_risk$method == "CV"
  ]
  data.frame(
    sample = result$sample_label,
    sample_size_per_training_group = result$sample_size,
    sample_size_training_group_1 = result$sample_size_e,
    sample_size_training_group_2 = result$sample_size_o,
    gamma_cv = result$gamma_cv,
    test_environments = nrow(current),
    in_Cgamma = sum(current$in_Cgamma),
    in_C0 = sum(current$in_C0),
    proportion_CV_better = mean(current$CV < current$OLS),
    worst_risk_OLS = worst_risk_ols,
    worst_risk_CV = worst_risk_cv,
    worst_risk_CV_minus_OLS = worst_risk_cv - worst_risk_ols
  )
}))
write.csv(
  summary_table,
  file.path(rna_output_directory, paste0(
    "rpe_subsample_Cgamma_C0_summary",
    output_suffix, ".csv"
  )),
  row.names = FALSE
)

cat("Training group imbalance:", training_imbalance, "\n")
cat("Training group 1:", paste(first_training, collapse = ", "), "\n")
cat("Training group 2:", paste(second_training, collapse = ", "), "\n")
cat("Membership tolerance:", membership_tolerance, "\n")
print(summary_table, row.names = FALSE)
cat("Worst-risk comparison at each gamma_CV:\n")
print(
  combined_worst_risks[, c(
    "sample_label", "gamma_cv", "method",
    "Rdelta_plus", "Rplus", "worst_risk_Rgamma"
  )],
  row.names = FALSE
)
cat("Scatterplots:", scatter_file, "\n")
cat("Boxplots:", boxplot_file, "\n")
