# Chamber screening and reporting: run all pairs and write final outputs.

write_screening_summary_pdf <- function(results, file) {
  # Create a PDF summary of all successful environment-pair results.
  successful <- results[results$status == "ok", , drop = FALSE]
  desired <- successful[successful$desired_pattern, , drop = FALSE]
  grDevices::pdf(file, width = 14, height = 9, family = "serif",
                 useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0.6, 0.6, 3.8, 0.6), family = "serif")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::title(
    "All-pairs Chamber screening: quadratic spline",
    line = 2.2
  )
  graphics::mtext(
    paste0(
      "degree 2; fixed knots (", fixed_knot_text, ")",
      "; gamma grid 0--3 by 0.01; alpha = ",
      format(significance_level, scientific = FALSE, trim = TRUE)
    ),
    side = 3, line = 0.6, cex = 0.76
  )
  overview <- c(
    paste("Intervened environments:", length(environment_files)),
    if (length(excluded_environment_files)) {
      paste(
        "Excluded ir_2-related environments:",
        length(excluded_environment_files)
      )
    },
    paste("Unordered pairs evaluated:", nrow(results)),
    paste("Successful fits:", nrow(successful)),
    paste("Pairs with gamma_CV > 0:", sum(successful$gamma_CV > 0)),
    paste(
      "Same sensitivity and higher CV specificity:",
      sum(successful$desired_pattern)
    ),
    paste(
      "Same sensitivity, higher CV specificity, all five causes found:",
      sum(successful$desired_pattern & successful$all_direct_found)
    )
  )
  graphics::text(0.04, 0.88, paste(overview, collapse = "\n"),
                 adj = c(0, 1), cex = 1.05)
  graphics::text(
    0.04, 0.55,
    paste0(
      "Desired pattern: gamma_CV > 0, equal sensitivity = TP/5, and ",
      "specificity_CV > specificity_OLS over 12 non-causes."
    ),
    adj = c(0, 1), cex = 0.9
  )

  if (!nrow(desired)) return(invisible(NULL))
  desired <- desired[order(
    -desired$all_direct_found,
    -desired$specificity_gain_count,
    -desired$CV_direct_rejected,
    desired$CV_noncauses_rejected,
    desired$pair_id
  ), , drop = FALSE]
  display_count <- min(25L, nrow(desired))
  display <- desired[seq_len(display_count), , drop = FALSE]

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, display_count + 2.7))
  graphics::title(
    paste0("Top ", display_count, " pairs with the desired pattern"),
    line = 2.2
  )
  graphics::mtext(
    "TP = rejected direct causes; FP = rejected non-causes",
    side = 3, line = 0.6, cex = 0.76
  )
  x <- c(0.01, 0.31, 0.64, 0.73, 0.81, 0.90, 0.98)
  headers <- c("Environment 1", "Environment 2", "gamma", "OLS TP",
               "CV TP", "OLS FP", "CV FP")
  align <- c(0, 0, rep(1, 5))
  header_y <- display_count + 1.3
  for (column in seq_along(x)) {
    graphics::text(x[column], header_y, headers[column],
                   adj = c(align[column], 0.5), font = 2, cex = 0.72)
  }
  graphics::segments(0.01, header_y - 0.45, 0.99, header_y - 0.45)
  shorten <- function(value) {
    # Shorten an environment filename for display in the PDF.
    value <- sub("^uniform_", "", value)
    sub("\\.csv$", "", value)
  }
  for (row in seq_len(display_count)) {
    y <- display_count - row + 0.55
    values <- c(
      shorten(display$environment_1[row]),
      shorten(display$environment_2[row]),
      format(display$gamma_CV[row], digits = 4),
      display$OLS_direct_rejected[row], display$CV_direct_rejected[row],
      display$OLS_noncauses_rejected[row],
      display$CV_noncauses_rejected[row]
    )
    for (column in seq_along(x)) {
      graphics::text(x[column], y, values[column],
                     adj = c(align[column], 0.5), cex = 0.65)
    }
  }
}

build_classification_recap <- function(results) {
  # Aggregate OLS and CR classification metrics across all pairs.
  successful <- results[results$status == "ok", , drop = FALSE]
  gamma_zero <- successful[successful$gamma_CV == 0, , drop = FALSE]
  metric_names <- c(
    Sensitivity = "sensitivity",
    Specificity = "specificity",
    Precision = "precision",
    F1 = "F1",
    `Balanced accuracy` = "balanced_accuracy",
    MCC = "MCC"
  )

  mean_metrics <- do.call(rbind, lapply(names(metric_names), function(label) {
    metric <- unname(metric_names[label])
    data.frame(
      Metric = label,
      OLS = mean(successful[[paste0(metric, "_OLS")]], na.rm = TRUE),
      CR = mean(successful[[paste0(metric, "_CV")]], na.rm = TRUE),
      check.names = FALSE
    )
  }))
  row.names(mean_metrics) <- NULL

  comparison_metrics <- metric_names[c(
    "Precision", "F1", "Balanced accuracy", "MCC"
  )]
  pairwise_comparison <- do.call(
    rbind,
    lapply(names(comparison_metrics), function(label) {
      metric <- unname(comparison_metrics[label])
      difference <- successful[[paste0(metric, "_CV")]] -
        successful[[paste0(metric, "_OLS")]]
      tolerance <- 1e-12
      data.frame(
        Metric = label,
        `CR better` = sum(difference > tolerance, na.rm = TRUE),
        Equal = sum(abs(difference) <= tolerance, na.rm = TRUE),
        `OLS better` = sum(difference < -tolerance, na.rm = TRUE),
        check.names = FALSE
      )
    })
  )
  row.names(pairwise_comparison) <- NULL

  gamma_zero_metrics <- do.call(
    rbind,
    lapply(names(metric_names), function(label) {
      metric <- unname(metric_names[label])
      OLS_values <- gamma_zero[[paste0(metric, "_OLS")]]
      CR_values <- gamma_zero[[paste0(metric, "_CV")]]
      finite_or_missing_equal <-
        (is.na(OLS_values) & is.na(CR_values)) |
        (!is.na(OLS_values) & !is.na(CR_values) &
           abs(OLS_values - CR_values) <= 1e-12)
      data.frame(
        Metric = label,
        `OLS = CV` = if (all(finite_or_missing_equal)) {
          mean(OLS_values, na.rm = TRUE)
        } else {
          NA_real_
        },
        check.names = FALSE
      )
    })
  )
  row.names(gamma_zero_metrics) <- NULL

  desired <- successful[successful$desired_pattern, , drop = FALSE]
  best_pair <- NULL
  if (nrow(desired)) {
    best <- desired[1L, , drop = FALSE]
    best_pair <- data.frame(
      Metric = names(metric_names),
      OLS = vapply(metric_names, function(metric) {
        best[[paste0(metric, "_OLS")]]
      }, numeric(1)),
      CR = vapply(metric_names, function(metric) {
        best[[paste0(metric, "_CV")]]
      }, numeric(1)),
      check.names = FALSE
    )
    attr(best_pair, "environment_1") <- best$environment_1
    attr(best_pair, "environment_2") <- best$environment_2
    attr(best_pair, "gamma_CV") <- best$gamma_CV
  }

  list(
    knots = fixed_knots,
    pairs_evaluated = nrow(results),
    successful_pairs = nrow(successful),
    positive_gamma_pairs = sum(successful$gamma_CV > 0),
    gamma_zero_pairs = nrow(gamma_zero),
    desired_pairs = sum(successful$desired_pattern),
    desired_all_direct_pairs = sum(
      successful$desired_pattern & successful$all_direct_found
    ),
    mean_metrics = mean_metrics,
    gamma_zero_metrics = gamma_zero_metrics,
    gamma_zero_results = gamma_zero,
    pairwise_comparison = pairwise_comparison,
    best_pair = best_pair
  )
}

print_classification_recap <- function(recap, digits = 4L) {
  # Print the aggregate classification results to the console.
  cat("\nCLASSIFICATION RECAP\n")
  cat("Quadratic knots:", paste(recap$knots, collapse = ", "), "\n")
  cat(
    "Successful pairs:", recap$successful_pairs, "/",
    recap$pairs_evaluated, "\n"
  )
  cat("Pairs with gamma_CV > 0:", recap$positive_gamma_pairs, "\n")
  cat("Desired-pattern pairs:", recap$desired_pairs, "\n")
  cat(
    "Desired-pattern pairs finding all direct causes:",
    recap$desired_all_direct_pairs, "\n\n"
  )

  cat("Mean metrics over all successful pairs\n")
  mean_table <- recap$mean_metrics
  mean_table[c("OLS", "CR")] <- round(mean_table[c("OLS", "CR")], digits)
  print(mean_table, row.names = FALSE)

  cat("\nPairwise comparison over all successful pairs\n")
  print(recap$pairwise_comparison, row.names = FALSE)

  cat(
    "\nMean metrics for the ", recap$gamma_zero_pairs,
    " pairs with gamma_CV = 0\n", sep = ""
  )
  gamma_zero_table <- recap$gamma_zero_metrics
  gamma_zero_table[["OLS = CV"]] <- round(
    gamma_zero_table[["OLS = CV"]], digits
  )
  print(gamma_zero_table, row.names = FALSE)

  if (!is.null(recap$best_pair)) {
    cat(
      "\nHighest-ranked desired pair:\n",
      attr(recap$best_pair, "environment_1"), " vs ",
      attr(recap$best_pair, "environment_2"),
      "; gamma_CV = ", attr(recap$best_pair, "gamma_CV"), "\n",
      sep = ""
    )
    best_table <- recap$best_pair
    best_table[c("OLS", "CR")] <- round(best_table[c("OLS", "CR")], digits)
    print(best_table, row.names = FALSE)
  }

  invisible(recap)
}

# Screen the environment pairs in parallel batches.
cat(
  "Screening ", nrow(pair_indices), " unordered pairs with ", worker_count,
  " worker(s)...\n", sep = ""
)
numbered_pairs <- cbind(pair_id = seq_len(nrow(pair_indices)), pair_indices)
batch_size <- max(8L, worker_count * 5L)
result_batches <- vector("list", ceiling(nrow(numbered_pairs) / batch_size))
batch_number <- 0L
start_time <- proc.time()[["elapsed"]]

for (batch_start in seq(1L, nrow(numbered_pairs), by = batch_size)) {
  batch_number <- batch_number + 1L
  batch_end <- min(nrow(numbered_pairs), batch_start + batch_size - 1L)
  current_pairs <- split(
    numbered_pairs[batch_start:batch_end, , drop = FALSE],
    seq_len(batch_end - batch_start + 1L)
  )
  current_results <- if (.Platform$OS.type == "unix" && worker_count > 1L) {
    parallel::mclapply(
      current_pairs, screen_pair,
      mc.cores = worker_count, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    lapply(current_pairs, screen_pair)
  }
  result_batches[[batch_number]] <- do.call(rbind, current_results)
  elapsed <- proc.time()[["elapsed"]] - start_time
  cat(
    "SCREEN_PROGRESS ", batch_end, "/", nrow(numbered_pairs),
    " elapsed_seconds=", format(round(elapsed, 1), nsmall = 1), "\n",
    sep = ""
  )
  flush.console()
}

screening_results <- do.call(rbind, result_batches)

# Sort the pairwise results and write the main CSV files.
screening_results$OLS_equals_CR <- screening_results$gamma_CV == 0
screening_results$OLS_equals_CR[screening_results$status != "ok"] <- NA
screening_results <- screening_results[order(
  screening_results$status != "ok",
  -screening_results$desired_pattern,
  -screening_results$all_direct_found,
  -screening_results$specificity_gain_count,
  -screening_results$CV_direct_rejected,
  screening_results$CV_noncauses_rejected,
  screening_results$pair_id
), , drop = FALSE]
row.names(screening_results) <- NULL

screening_csv <- file.path(
  csv_directory,
  "chamber_quadratic_all_pairs_screening.csv"
)
utils::write.csv(screening_results, screening_csv, row.names = FALSE)

classification_columns <- c(
  "pair_id", "environment_1", "environment_2", "gamma_CV",
  "sensitivity_OLS", "sensitivity_CV",
  "specificity_OLS", "specificity_CV",
  "precision_OLS", "precision_CV",
  "F1_OLS", "F1_CV",
  "balanced_accuracy_OLS", "balanced_accuracy_CV",
  "MCC_OLS", "MCC_CV", "OLS_equals_CR", "status"
)
classification_metrics_csv <- file.path(
  csv_directory,
  "chamber_quadratic_all_pairs_classification_metrics.csv"
)
utils::write.csv(
  screening_results[, classification_columns, drop = FALSE],
  classification_metrics_csv,
  row.names = FALSE
)

positive_gamma_results <- screening_results[
  screening_results$status == "ok" & screening_results$gamma_CV > 0,
  , drop = FALSE
]
positive_gamma_csv <- file.path(
  csv_directory,
  "chamber_quadratic_all_pairs_positive_gamma.csv"
)
utils::write.csv(positive_gamma_results, positive_gamma_csv, row.names = FALSE)

desired_pattern_results <- positive_gamma_results[
  positive_gamma_results$desired_pattern,
  , drop = FALSE
]
desired_pattern_csv <- file.path(
  csv_directory,
  paste0(
    "chamber_quadratic_same_sensitivity_higher_specificity.csv"
  )
)
utils::write.csv(desired_pattern_results, desired_pattern_csv, row.names = FALSE)

summary_pdf <- file.path(
  pdf_directory,
  "chamber_quadratic_all_pairs_screening_summary.pdf"
)

# Write a visual screening summary.
write_screening_summary_pdf(screening_results, summary_pdf)

desired_results <- screening_results[
  screening_results$status == "ok" & screening_results$desired_pattern,
  , drop = FALSE
]
if (nrow(desired_results)) {
  # Refit the best pair and write its detailed Hotelling test table.
  best <- desired_results[1L, , drop = FALSE]
  best_pair_row <- c(
    best$pair_id,
    match(best$environment_1, environment_files),
    match(best$environment_2, environment_files)
  )
  best_fit <- screen_pair(best_pair_row, retain_table = TRUE)
  best_stem <- paste0(
    "chamber_quadratic_best_pair_",
    tools::file_path_sans_ext(sub("^uniform_", "", best$environment_1)),
    "__",
    tools::file_path_sans_ext(sub("^uniform_", "", best$environment_2)),
    "_hotelling_t_squared.pdf"
  )
  best_pdf <- file.path(pdf_directory, best_stem)
  write_chamber_hotelling_t_squared_tests_pdf(
    best_fit$table,
    best_pdf,
    gamma = best$gamma_CV,
    model_description = paste0(
      "quadratic spline degree 2; knots ", fixed_knot_text,
      "; ",
      tools::file_path_sans_ext(sub("^uniform_", "", best$environment_1)),
      " vs ",
      tools::file_path_sans_ext(sub("^uniform_", "", best$environment_2))
    )
  )
  cat("BEST_PAIR_PDF=", best_pdf, "\n", sep = "")
}

successful_count <- sum(screening_results$status == "ok")
desired_count <- sum(screening_results$desired_pattern)
complete_desired_count <- sum(
  screening_results$desired_pattern & screening_results$all_direct_found
)
classification_recap <- build_classification_recap(screening_results)

# Print the final recap and all output locations.
gamma_zero_results <- classification_recap$gamma_zero_results
gamma_zero_metrics <- classification_recap$gamma_zero_metrics
print_classification_recap(classification_recap)

cat("SCREENING_CSV=", screening_csv, "\n", sep = "")
cat(
  "CLASSIFICATION_METRICS_CSV=", classification_metrics_csv,
  "\n", sep = ""
)
cat("POSITIVE_GAMMA_CSV=", positive_gamma_csv, "\n", sep = "")
cat("DESIRED_PATTERN_CSV=", desired_pattern_csv, "\n", sep = "")
cat("SUMMARY_PDF=", summary_pdf, "\n", sep = "")
cat("SUCCESSFUL_PAIRS=", successful_count, "\n", sep = "")
cat("DESIRED_PAIRS=", desired_count, "\n", sep = "")
cat("DESIRED_ALL_DIRECT_PAIRS=", complete_desired_count, "\n", sep = "")
