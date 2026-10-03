# Simulation experiments for Figures 3, 5, 7, and 9.

source("code/functions/datautil.R")
source("code/functions/cd.R")
source("code/functions/ols.R")
source("code/functions/cross_valid_fun_n.R")
source("code/functions/visualizations.R")
source("code/functions/confidence_intervals.R")

run_paper_example <- function(example) {
  settings <- list(
    A = list(figure = "03", n = 100L, seed = 123L, ci_seed = 2026L,
             max_gamma = 30, path_max = 30, cv_max = 15),
    B = list(figure = "05", n = 100L, seed = 134L, ci_seed = 2027L,
             max_gamma = 300, path_max = 50, cv_max = 15),
    C = list(figure = "07", n = 1000L, seed = 123L, ci_seed = 2027L,
             max_gamma = 30, path_max = 30, cv_max = 30),
    D = list(figure = "09", n = 1000L, seed = 123L, ci_seed = 2028L,
             max_gamma = 30, path_max = 30, cv_max = 30)
  )
  config <- settings[[example]]
  if (is.null(config)) stop("Choose example A, B, C or D.")
  simulate <- get(paste0("simulate_example_", tolower(example), "_data"))
  set.seed(config$seed)
  data <- simulate(config$n)
  colnames(data$Xe) <- colnames(data$Xo) <- c("X1", "X2")
  gamma_grid <- seq(0, config$max_gamma, by = 0.1)
  estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
  if (example == "D") set.seed(456L)
  cv <- cross_validation(5L, data, estimators)
  gamma_cv <- as.numeric(cv$best_gamma)
  m <- build_matrices(data)
  beta_rl <- drop(compute_cd(m))
  beta_ols <- drop(compute_ols(m))
  path <- estimate_beta_path_ci(data, gamma_grid)
  repetitions <- as.integer(Sys.getenv("PAPER_EXAMPLE_CI_REPS", "500"))
  if (is.na(repetitions) || repetitions < 2L) {
    stop("PAPER_EXAMPLE_CI_REPS must be an integer of at least two.")
  }
  intervals <- beta_path_confidence_intervals(
    simulate, config$n, gamma_grid, repetitions = repetitions,
    seed = config$ci_seed
  )
  prefix <- file.path("output", "simulations", paste0("figure_", config$figure))
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  truth <- if (example == "B") NULL else c(1, 0)
  plot_beta_path(
    path, gamma_grid, beta_rl, gamma_cv, beta_true = truth,
    file = paste0(prefix, "_a.pdf"), x_max = config$path_max,
    beta_lower = intervals$beta_lower, beta_upper = intervals$beta_upper,
    reference_label = if (example %in% c("C", "D")) "CP" else "RL",
    truth_label = "CP"
  )
  plot_cv_risk(cv, gamma_grid, file = paste0(prefix, "_b.pdf"),
               x_max = config$cv_max)
  chosen <- which.min(abs(gamma_grid - gamma_cv))
  estimates <- data.frame(
    coefficient = c("X1", "X2"), OLS = beta_ols,
    CR = path[, chosen], RL = beta_rl, gamma_cv = gamma_cv
  )
  write.csv(estimates, paste0(prefix, "_estimates.csv"), row.names = FALSE)
  write.csv(data.frame(
    gamma = gamma_grid, beta_1 = path[1, ], beta_2 = path[2, ],
    beta_1_lower = intervals$beta_lower[1, ],
    beta_1_upper = intervals$beta_upper[1, ],
    beta_2_lower = intervals$beta_lower[2, ],
    beta_2_upper = intervals$beta_upper[2, ]
  ), paste0(prefix, "_path.csv"), row.names = FALSE)
  write.csv(data.frame(gamma = gamma_grid, t(cv$foldwise_rdelta),
                       average = colMeans(cv$foldwise_rdelta)),
            paste0(prefix, "_cv.csv"), row.names = FALSE)
  saveRDS(list(data = data, cv = cv, config = config, estimates = estimates,
               intervals = intervals, gamma_grid = gamma_grid, path = path),
          paste0(prefix, "_results.rds"))
  cat("Example", example, "gamma_CV =", gamma_cv, "\n")
  print(estimates, row.names = FALSE)
  invisible(estimates)
}
