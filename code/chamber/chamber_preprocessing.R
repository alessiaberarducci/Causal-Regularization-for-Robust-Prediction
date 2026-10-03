# Chamber preprocessing and spline design matrices.

Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

source("code/functions/datautil.R")
source("code/functions/cross_valid_fun_n.R")
source("code/functions/chamber_hotelling_t_squared_tests.R")

data_directory <- "data/lt_interventions_standard_v1"
fixed_knot_specification <- Sys.getenv(
  "CHAMBER_PAIR_FIXED_KNOTS", "-1,0,1"
)
fixed_knots <- suppressWarnings(as.numeric(strsplit(
  fixed_knot_specification, ",", fixed = TRUE
)[[1]]))
if (!length(fixed_knots) || any(!is.finite(fixed_knots))) {
  stop("CHAMBER_PAIR_FIXED_KNOTS must be a comma-separated numeric vector.")
}
fixed_knots <- sort(unique(fixed_knots))
supported_knots <- identical(fixed_knots, c(-1, 0, 1)) ||
  identical(fixed_knots, 0)
if (!supported_knots) {
  stop("CHAMBER_PAIR_FIXED_KNOTS must be either -1,0,1 or 0.")
}
fixed_knot_text <- paste(format(
  fixed_knots, scientific = FALSE, trim = TRUE
), collapse = ", ")
fixed_knot_slug <- paste(vapply(fixed_knots, function(knot) {
  token <- gsub(
    ".", "p",
    format(abs(knot), scientific = FALSE, trim = TRUE),
    fixed = TRUE
  )
  if (knot < 0) paste0("minus", token) else token
}, character(1)), collapse = "_")
spline_variant <- paste0("fixed_knots_", fixed_knot_slug)
exclude_ir2_related_environments <- tolower(Sys.getenv(
  "CHAMBER_PAIR_EXCLUDE_IR2_RELATED", "true"
)) %in% c("true", "1", "yes")
default_output_parent <- file.path(
  "output", "chamber_output", "quadratic_splines",
  spline_variant
)
if (exclude_ir2_related_environments) {
  default_output_parent <- file.path(
    default_output_parent, "exclude_ir2_related"
  )
}
default_output_directory <- file.path(
  default_output_parent, "all_pairs_gamma_grid_0_3_step_0.01"
)
output_directory <- Sys.getenv(
  "CHAMBER_PAIR_OUTPUT_DIRECTORY", default_output_directory
)
csv_directory <- file.path(output_directory, "csv")
pdf_directory <- file.path(output_directory, "pdf")
dir.create(csv_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_directory, recursive = TRUE, showWarnings = FALSE)

target <- "ir_2"

predictors <- c(
  "red", "green", "blue", "current", "pol_1", "pol_2",
  "ir_1", "vis_1", "vis_2", "ir_3", "vis_3",
  "l_11", "l_12", "l_21", "l_22", "l_31", "l_32"
)
direct_causes <- c("red", "green", "blue", "l_21", "l_22")
non_causes <- setdiff(predictors, direct_causes)
required_columns <- c(predictors, target)
gamma_grid <- seq(0, 3, by = 0.01)
terms_per_predictor <-
  2L + length(fixed_knots)
fold_count <- 5L
significance_level <- suppressWarnings(as.numeric(Sys.getenv(
  "CHAMBER_PAIR_HOTELLING_ALPHA", "0.05"
)))
if (!is.finite(significance_level) ||
    significance_level <= 0 || significance_level >= 1) {
  stop("CHAMBER_PAIR_HOTELLING_ALPHA must be strictly between 0 and 1.")
}
pseudoinverse_tolerance <- 1e-8

environment_files <- sort(list.files(
  data_directory,
  pattern = "^uniform_.*\\.csv$",
  full.names = FALSE
))

environment_files <- setdiff(environment_files, "uniform_reference.csv")
excluded_environment_files <- character(0)
if (exclude_ir2_related_environments) {
  excluded_environment_files <- grep(
    "^uniform_(diode|t)_ir_2_.*\\.csv$",
    environment_files,
    value = TRUE
  )
  environment_files <- setdiff(
    environment_files, excluded_environment_files
  )
}
if (length(environment_files) < 2L) {
  stop("At least two intervened Chamber environments are required.")
}

read_environment <- function(file) {
  current <- read.csv(
    file.path(data_directory, file),
    check.names = FALSE
  )
  missing_columns <- setdiff(required_columns, names(current))
  if (length(missing_columns)) {
    stop(file, " is missing: ", paste(missing_columns, collapse = ", "))
  }
  current <- current[, required_columns, drop = FALSE]
  if (!all(vapply(current, is.numeric, logical(1)))) {
    stop(file, " contains a non-numeric required column.")
  }
  current
}

if (length(excluded_environment_files)) {
  cat(
    "Excluding ", length(excluded_environment_files),
    " ir_2-related environments: ",
    paste(excluded_environment_files, collapse = ", "), "\n",
    sep = ""
  )
}
cat("Loading ", length(environment_files), " intervened environments...\n", sep = "")
environment_data <- setNames(
  lapply(environment_files, read_environment),
  environment_files
)

pair_indices <- t(utils::combn(seq_along(environment_files), 2L))
requested_environment_1 <- Sys.getenv("CHAMBER_PAIR_ENVIRONMENT_1", "")
requested_environment_2 <- Sys.getenv("CHAMBER_PAIR_ENVIRONMENT_2", "")
if (nzchar(requested_environment_1) || nzchar(requested_environment_2)) {
  if (!nzchar(requested_environment_1) || !nzchar(requested_environment_2)) {
    stop(
      "CHAMBER_PAIR_ENVIRONMENT_1 and CHAMBER_PAIR_ENVIRONMENT_2 must ",
      "be supplied together."
    )
  }
  requested_indices <- match(
    c(requested_environment_1, requested_environment_2),
    environment_files
  )
  if (anyNA(requested_indices) || requested_indices[1L] == requested_indices[2L]) {
    stop("The requested validation pair is invalid.")
  }
  pair_indices <- matrix(sort(requested_indices), nrow = 1L)
}
pair_limit <- suppressWarnings(as.integer(Sys.getenv(
  "CHAMBER_PAIR_LIMIT", as.character(nrow(pair_indices))
)))
if (is.na(pair_limit) || pair_limit < 1L) {
  stop("CHAMBER_PAIR_LIMIT must be a positive integer.")
}
pair_indices <- pair_indices[seq_len(min(pair_limit, nrow(pair_indices))), , drop = FALSE]

default_workers <- if (.Platform$OS.type == "unix") 8L else 1L
worker_count <- suppressWarnings(as.integer(Sys.getenv(
  "CHAMBER_PAIR_WORKERS", as.character(default_workers)
)))
if (is.na(worker_count) || worker_count < 1L) {
  stop("CHAMBER_PAIR_WORKERS must be a positive integer.")
}
worker_count <- min(worker_count, nrow(pair_indices))

symmetric_pseudoinverse_solve <- function(matrix, right_hand_side) {
  matrix <- 0.5 * (matrix + t(matrix))
  decomposition <- eigen(matrix, symmetric = TRUE)
  eigenvalue_scale <- max(abs(decomposition$values))
  if (!is.finite(eigenvalue_scale) || eigenvalue_scale == 0) {
    stop("A moment matrix has numerical rank zero.")
  }
  retained <- abs(decomposition$values) >
    pseudoinverse_tolerance * eigenvalue_scale
  if (!any(retained)) stop("A moment matrix has numerical rank zero.")
  retained_vectors <- decomposition$vectors[, retained, drop = FALSE]
  projected_right_hand_side <-
    t(retained_vectors) %*% right_hand_side
  scaled_projection <- sweep(
    projected_right_hand_side,
    1L,
    decomposition$values[retained],
    "/"
  )
  retained_vectors %*% scaled_projection
}

raw_truncated_power_spline_basis <- function(X, knots_by_predictor) {
  expanded <- lapply(seq_len(ncol(X)), function(column_index) {
    x <- X[, column_index]
    polynomial_terms <- vapply(
      1:2,
      function(power) x^power,
      numeric(length(x))
    )
    hinge_terms <- vapply(
      knots_by_predictor[[column_index]],
      function(knot) pmax(x - knot, 0)^2,
      numeric(length(x))
    )
    cbind(polynomial_terms, hinge_terms)
  })
  do.call(cbind, expanded)
}

prepare_pair <- function(environment_1, environment_2) {

  complete_1 <- complete.cases(environment_1[, required_columns, drop = FALSE])
  complete_2 <- complete.cases(environment_2[, required_columns, drop = FALSE])
  environment_1 <- environment_1[complete_1, , drop = FALSE]
  environment_2 <- environment_2[complete_2, , drop = FALSE]

  pooled_predictors <- rbind(
    environment_1[, predictors, drop = FALSE],
    environment_2[, predictors, drop = FALSE]
  )

  predictor_center <- colMeans(pooled_predictors)
  predictor_scale <- vapply(pooled_predictors, stats::sd, numeric(1))
  if (any(!is.finite(predictor_scale) | predictor_scale <= 1e-12)) {
    stop("At least one original predictor has zero numerical variance.")
  }
  response_center <- mean(c(
    environment_1[[target]], environment_2[[target]]
  ))

  standardize_predictors <- function(current) {
    sweep(
      sweep(as.matrix(current[, predictors, drop = FALSE]),
            2L, predictor_center, "-"),
      2L, predictor_scale, "/"
    )
  }
  standardized_1 <- standardize_predictors(environment_1)
  standardized_2 <- standardize_predictors(environment_2)

  knots_by_predictor <- rep(
    list(fixed_knots), length(predictors)
  )
  raw_basis_1 <- raw_truncated_power_spline_basis(
    standardized_1, knots_by_predictor
  )
  raw_basis_2 <- raw_truncated_power_spline_basis(
    standardized_2, knots_by_predictor
  )
  pooled_raw_basis <- rbind(raw_basis_1, raw_basis_2)

  basis_center <- colMeans(pooled_raw_basis)
  basis_scale <- apply(pooled_raw_basis, 2L, stats::sd)
  if (any(!is.finite(basis_scale) | basis_scale <= 1e-12)) {
    stop("At least one spline basis column has zero numerical variance.")
  }
  standardize_basis <- function(raw_basis) {
    sweep(
      sweep(raw_basis, 2L, basis_center, "-"),
      2L, basis_scale, "/"
    )
  }

  model_data <- list(
    Xe = standardize_basis(raw_basis_1),
    ye = environment_1[[target]] - response_center,
    Xo = standardize_basis(raw_basis_2),
    yo = environment_2[[target]] - response_center
  )
  term_counts <- 2L + lengths(knots_by_predictor)

  group_ends <- cumsum(term_counts)
  group_starts <- c(1L, head(group_ends, -1L) + 1L)
  attr(model_data, "coefficient_groups") <- setNames(
    Map(seq.int, group_starts, group_ends), predictors
  )
  model_data
}
