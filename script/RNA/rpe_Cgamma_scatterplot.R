# RPE paired OLS/CV risks with C_gamma and C_0 membership

source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/rpe_robust_functions.R")

set.seed(123)

# Allow small empirical violations of the positive-semidefinite constraint.
membership_tolerance <- as.numeric(
  Sys.getenv("RPE_MEMBERSHIP_TOLERANCE", "0.05")
)
if (!is.finite(membership_tolerance) || membership_tolerance < 0) {
  stop("RPE_MEMBERSHIP_TOLERANCE must be a nonnegative finite number.")
}

# Use the same variables and training-environment construction as
# rpe_robust_analysis.R.
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

dat <- read.csv("dataset_rpe.csv", check.names = FALSE)
dat <- filter_environments(dat, environment_column, minimum_cells = 100)
available_environments <- unique(as.character(dat[[environment_column]]))

# Explicit full-sample training groups, matching alternative n0.
first_training <- c(
  "ENSG00000147604", "ENSG00000110700",
  "ENSG00000133112", "ENSG00000125691"
)
second_training <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000172757",
  "ENSG00000067225", "ENSG00000108518"
)
training_interventions <- c(first_training, second_training)

training <- make_two_environment_data(
  dat, environment_column, first_training, second_training,
  predictor_genes, target_gene
)
raw_data_train <- training$data
centering <- compute_common_centering(raw_data_train)
data_train <- centering$data

# Every remaining eligible intervention is used as a test environment.
test_environments <- setdiff(
  available_environments,
  c(training_interventions, target_gene, control_environment, "excluded")
)

# Select gamma using the same centered 10-fold CV specification.
folds <- 10L
gamma_grid <- seq(0, 40, by = 0.01)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
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

# Refit pooled OLS and the selected CV estimator on all training observations.
training_moments <- moments(data_train)
beta_ols <- setNames(drop(compute_ols(training_moments)), predictor_genes)
beta_cv <- setNames(
  drop(make_estimator(gamma_cv)(data_train)), predictor_genes
)
beta_rl <- setNames(drop(compute_cd(training_moments)), predictor_genes)
intercept_ols <- recover_intercept(beta_ols, centering)
intercept_cv <- recover_intercept(beta_cv, centering)

test_risks <- compute_test_risks(
  dat, test_environments, environment_column,
  beta_ols, beta_cv, predictor_genes, target_gene,
  intercept_ols = intercept_ols, intercept_cv = intercept_cv
)

# Construct the augmented second-moment matrices defining C_gamma.
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

M_0 <- 0.5 * M_plus
M_gamma <- M_0 + 0.5 * gamma_cv * M_delta_plus

compute_environment_moment <- function(environment) {
  current <- dat[
    as.character(dat[[environment_column]]) == environment, , drop = FALSE
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

environment_moments <- setNames(
  lapply(test_environments, compute_environment_moment),
  test_environments
)
minimum_eigen_Cgamma <- vapply(environment_moments, function(M_test) {
  min(eigen(M_gamma - M_test, symmetric = TRUE, only.values = TRUE)$values)
}, numeric(1))
minimum_eigen_C0 <- vapply(environment_moments, function(M_test) {
  min(eigen(M_0 - M_test, symmetric = TRUE, only.values = TRUE)$values)
}, numeric(1))

membership <- data.frame(
  environment = names(environment_moments),
  minimum_eigenvalue_Cgamma = minimum_eigen_Cgamma,
  minimum_eigenvalue_C0 = minimum_eigen_C0,
  membership_tolerance = membership_tolerance,
  in_Cgamma = minimum_eigen_Cgamma >= -membership_tolerance,
  in_C0 = minimum_eigen_C0 >= -membership_tolerance,
  row.names = NULL
)

scatter_data <- merge(
  test_risks, membership,
  by = "environment", all = FALSE, sort = FALSE
)
scatter_data$class <- ifelse(
  scatter_data$in_C0, "C0",
  ifelse(scatter_data$in_Cgamma, "Cgamma", "outside")
)
scatter_data$class <- factor(
  scatter_data$class,
  levels = c("outside", "Cgamma", "C0")
)

dir.create("visualization", showWarnings = FALSE)
write.csv(
  scatter_data,
  "visualization/rpe_OLS_CV_Cgamma_C0_scatter_data.csv",
  row.names = FALSE
)

plot_file <- "visualization/rpe_OLS_CV_Cgamma_C0_scatterplot.pdf"
grDevices::pdf(
  plot_file, width = 6.2, height = 5.5,
  family = "serif", useDingbats = FALSE
)
old_par <- graphics::par(
  family = "serif", mar = c(4.8, 5, 1.2, 1),
  lend = "round"
)

limits <- range(scatter_data$OLS, scatter_data$CV, finite = TRUE)
padding <- max(0.01, 0.04 * diff(limits))
limits <- c(max(0, limits[1] - padding), limits[2] + padding)
scatter_colours <- c(
  outside = "black",
  Cgamma = "#0072B2",
  C0 = "#D55E00"
)
boxplot_colours <- c(
  outside = "grey65",
  Cgamma = "#0072B2",
  C0 = "#D55E00"
)

graphics::plot(
  scatter_data$OLS, scatter_data$CV,
  type = "n", xlim = limits, ylim = limits, asp = 1,
  xlab = expression(MSE[e](hat(beta)[OLS])),
  ylab = expression(MSE[e](hat(beta)[CV]))
)
graphics::abline(a = 0, b = 1, lty = 2, col = "grey35")

# Draw nested classes in this order so C_0 remains visible above C_gamma.
for (current_class in levels(scatter_data$class)) {
  selected <- scatter_data$class == current_class
  graphics::points(
    scatter_data$OLS[selected], scatter_data$CV[selected],
    pch = 16, cex = 0.65,
    col = grDevices::adjustcolor(
      scatter_colours[[current_class]], alpha.f = 0.75
    )
  )
}

class_counts <- table(scatter_data$class)
graphics::legend(
  "topleft",
  legend = c(
    sprintf("Outside Cgamma (n = %d)", class_counts[["outside"]]),
    sprintf("Cgamma without C0 (n = %d)", class_counts[["Cgamma"]]),
    sprintf("C0 (n = %d)", class_counts[["C0"]])
  ),
  col = scatter_colours,
  pch = 16, pt.cex = 0.8, bty = "n", cex = 0.78
)
graphics::mtext(
  bquote(gamma[CV] == .(gamma_cv)),
  side = 3, adj = 1, line = 0.1, cex = 0.8
)

graphics::par(old_par)
grDevices::dev.off()

# Boxplot over all test environments, retaining the same membership colours.
boxplot_file <- "visualization/rpe_OLS_CV_Cgamma_C0_boxplot.pdf"
grDevices::pdf(
  boxplot_file, width = 6.2, height = 5.2,
  family = "serif", useDingbats = FALSE
)
old_par <- graphics::par(
  family = "serif", mar = c(4.5, 5, 1.2, 1),
  lend = "round"
)

boxplot_values <- scatter_data[, c("OLS", "CV")]
boxplot_limits <- range(boxplot_values, finite = TRUE)
boxplot_padding <- max(0.01, 0.04 * diff(boxplot_limits))
boxplot_limits <- c(
  max(0, boxplot_limits[1] - boxplot_padding),
  boxplot_limits[2] + boxplot_padding
)
graphics::boxplot(
  boxplot_values,
  names = c("OLS", "CV"),
  outline = FALSE, ylim = boxplot_limits,
  ylab = expression(MSE[e](hat(beta)))
)

set.seed(123)
jitter_offset <- stats::runif(nrow(scatter_data), -0.12, 0.12)
point_colours <- grDevices::adjustcolor(
  unname(boxplot_colours[as.character(scatter_data$class)]),
  alpha.f = 0.65
)
graphics::points(
  1 + jitter_offset, scatter_data$OLS,
  pch = 16, cex = 0.5, col = point_colours
)
graphics::points(
  2 + jitter_offset, scatter_data$CV,
  pch = 16, cex = 0.5, col = point_colours
)
graphics::legend(
  "topright",
  legend = c(
    sprintf("Outside Cgamma (n = %d)", class_counts[["outside"]]),
    sprintf("Cgamma without C0 (n = %d)", class_counts[["Cgamma"]]),
    sprintf("C0 (n = %d)", class_counts[["C0"]])
  ),
  col = boxplot_colours,
  pch = 16, pt.cex = 0.8, bty = "n", cex = 0.76
)

graphics::par(old_par)
grDevices::dev.off()

cat("Training group 1:", paste(first_training, collapse = ", "), "\n")
cat("Training group 2:", paste(second_training, collapse = ", "), "\n")
cat("Selected gamma:", gamma_cv, "\n")
cat("Membership tolerance:", membership_tolerance, "\n")
cat("Test environments:", nrow(scatter_data), "\n")
cat("In Cgamma:", sum(scatter_data$in_Cgamma), "\n")
cat("In C0:", sum(scatter_data$in_C0), "\n")
cat("Scatterplot:", plot_file, "\n")
cat("Boxplot:", boxplot_file, "\n")
