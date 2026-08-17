# Repeat the subsample experiment with explicitly specified training groups.

first_training <- "non-targeting"
second_training <- c(
  "ENSG00000133112", "ENSG00000075624",
  "ENSG00000187514", "ENSG00000067225"
)

Sys.setenv(
  RPE_FIRST_TRAINING = paste(first_training, collapse = ","),
  RPE_SECOND_TRAINING = paste(second_training, collapse = ","),
  RPE_OUTPUT_SUFFIX = "_alternative_n2",
  RPE_SAMPLE_SIZES = "850,650,500,250,150,100"
)
source("script/RNA/rpe_Cgamma_subsample_comparison.R")

# Save the full-sample and n = 100 scatterplots as separate article figures.
selected_scatter_results <- list(
  full_sample = experiment_results[["full"]],
  n100 = experiment_results[["100"]]
)
selected_scatter_values <- unlist(lapply(
  selected_scatter_results,
  function(result) c(result$data$OLS, result$data$CV)
))
selected_scatter_limits <- range(selected_scatter_values, finite = TRUE)
selected_scatter_padding <- max(
  0.01, 0.02 * diff(selected_scatter_limits)
)
selected_scatter_limits <- c(
  max(0, selected_scatter_limits[1] - selected_scatter_padding),
  selected_scatter_limits[2] + selected_scatter_padding
)

write_single_scatterplot <- function(result, file) {
  grDevices::pdf(
    file, width = 6.3, height = 6.0,
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  old_par <- graphics::par(
    family = "serif", mar = c(4.7, 4.7, 2.5, 0.8), lend = "round"
  )
  on.exit(graphics::par(old_par), add = TRUE)
  current_data <- result$data
  class_counts <- table(current_data$class)
  graphics::plot(
    current_data$OLS, current_data$CV,
    type = "n", xlim = selected_scatter_limits,
    ylim = selected_scatter_limits, asp = 1,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CV])),
    main = bquote(.(result$sample_label) ~ "," ~ gamma[CV] == .(
      result$gamma_cv
    ))
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey35")
  for (current_class in levels(current_data$class)) {
    selected <- current_data$class == current_class
    graphics::points(
      current_data$OLS[selected], current_data$CV[selected],
      pch = 16, cex = 0.65,
      col = grDevices::adjustcolor(
        scatter_colours[[current_class]], alpha.f = 0.72
      )
    )
  }
  graphics::legend(
    "topleft",
    legend = as.expression(list(
      bquote(
        "Outside " * C[gamma] * ": " * .(class_counts[["outside"]])
      ),
      bquote(C[gamma] - C[0] * ": " * .(class_counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(class_counts[["C0"]]))
    )),
    col = scatter_colours, pch = 16,
    pt.cex = 0.75, bty = "n", cex = 0.72
  )
}

full_sample_scatter_file <- file.path(
  rna_output_directory,
  "rpe_subsample_OLS_CV_scatterplot_full_sample_alternative_n2.pdf"
)
n100_scatter_file <- file.path(
  rna_output_directory,
  "rpe_subsample_OLS_CV_scatterplot_n100_alternative_n2.pdf"
)
write_single_scatterplot(
  selected_scatter_results$full_sample, full_sample_scatter_file
)
write_single_scatterplot(selected_scatter_results$n100, n100_scatter_file)

# Log-log versions retain every environment while expanding the dense region.
if (any(selected_scatter_values <= 0)) {
  stop("Log-scale scatterplots require strictly positive MSE values.")
}
selected_log_limits <- range(log10(selected_scatter_values), finite = TRUE)
selected_log_padding <- max(0.01, 0.02 * diff(selected_log_limits))
selected_log_limits <- 10^(selected_log_limits + c(-1, 1) * selected_log_padding)

write_single_log_scatterplot <- function(result, file) {
  grDevices::pdf(
    file, width = 6.3, height = 6.0,
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  old_par <- graphics::par(
    family = "serif", mar = c(4.7, 4.7, 2.5, 0.8), lend = "round"
  )
  on.exit(graphics::par(old_par), add = TRUE)
  current_data <- result$data
  class_counts <- table(current_data$class)
  graphics::plot(
    current_data$OLS, current_data$CV,
    type = "n", log = "xy",
    xlim = selected_log_limits, ylim = selected_log_limits, asp = 1,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CV])),
    main = bquote(.(result$sample_label) ~ "," ~ gamma[CV] == .(
      result$gamma_cv
    ))
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey35")
  for (current_class in levels(current_data$class)) {
    selected <- current_data$class == current_class
    graphics::points(
      current_data$OLS[selected], current_data$CV[selected],
      pch = 16, cex = 0.65,
      col = grDevices::adjustcolor(
        scatter_colours[[current_class]], alpha.f = 0.72
      )
    )
  }
  graphics::legend(
    "topleft",
    legend = as.expression(list(
      bquote(
        "Outside " * C[gamma] * ": " * .(class_counts[["outside"]])
      ),
      bquote(C[gamma] - C[0] * ": " * .(class_counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(class_counts[["C0"]]))
    )),
    col = scatter_colours, pch = 16,
    pt.cex = 0.75, bty = "n", cex = 0.72
  )
}

full_sample_log_scatter_file <- file.path(
  rna_output_directory,
  "rpe_subsample_OLS_CV_scatterplot_log_full_sample_alternative_n2.pdf"
)
n100_log_scatter_file <- file.path(
  rna_output_directory,
  "rpe_subsample_OLS_CV_scatterplot_log_n100_alternative_n2.pdf"
)
write_single_log_scatterplot(
  selected_scatter_results$full_sample, full_sample_log_scatter_file
)
write_single_log_scatterplot(
  selected_scatter_results$n100, n100_log_scatter_file
)

# Plot the paired MSE differences for the same two sample sizes.
selected_mse_differences <- lapply(
  selected_scatter_results,
  function(result) result$data$OLS - result$data$CV
)

write_mse_difference_histogram <- function(differences, file) {
  differences <- differences[is.finite(differences)]
  observed_limits <- range(differences)
  displayed_span <- diff(observed_limits)
  histogram_padding <- max(0.001, 0.04 * displayed_span)
  histogram_limits <- c(
    observed_limits[1] - histogram_padding,
    observed_limits[2] + histogram_padding
  )
  histogram_data <- graphics::hist(
    differences, breaks = 80, plot = FALSE
  )
  grDevices::pdf(
    file, width = 6.3, height = 5.2,
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  old_par <- graphics::par(
    family = "serif", mar = c(4.8, 4.8, 0.8, 0.8), lend = "round"
  )
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::plot(
    histogram_data,
    xlim = histogram_limits,
    xaxt = "n",
    main = "",
    xlab = expression(
      MSE[e](hat(beta)[OLS]) - MSE[e](hat(beta)[CV])
    ),
    ylab = "Frequency",
    col = "grey80", border = "white"
  )
  x_ticks <- pretty(histogram_limits, n = 5)
  if (observed_limits[1] < 0) {
    negative_tick <- floor(100 * observed_limits[1]) / 100
    x_ticks <- c(x_ticks, negative_tick)
  }
  x_ticks <- sort(unique(x_ticks[
    x_ticks >= histogram_limits[1] & x_ticks <= histogram_limits[2]
  ]))
  graphics::axis(1, at = x_ticks)
  graphics::abline(v = 0, lty = 2, lwd = 1.2, col = "black")
  graphics::rug(
    differences, side = 1, ticksize = 0.025,
    col = grDevices::adjustcolor("black", alpha.f = 0.55)
  )
}

full_sample_histogram_file <- file.path(
  rna_output_directory,
  "rpe_OLS_minus_CV_MSE_histogram_full_sample_alternative_n2.pdf"
)
n100_histogram_file <- file.path(
  rna_output_directory,
  "rpe_OLS_minus_CV_MSE_histogram_n100_alternative_n2.pdf"
)
write_mse_difference_histogram(
  selected_mse_differences$full_sample, full_sample_histogram_file
)
write_mse_difference_histogram(
  selected_mse_differences$n100, n100_histogram_file
)

cat("Full-sample scatterplot:", full_sample_scatter_file, "\n")
cat("n = 100 scatterplot:", n100_scatter_file, "\n")
cat("Full-sample log scatterplot:", full_sample_log_scatter_file, "\n")
cat("n = 100 log scatterplot:", n100_log_scatter_file, "\n")
cat("Full-sample MSE-difference histogram:", full_sample_histogram_file, "\n")
cat("n = 100 MSE-difference histogram:", n100_histogram_file, "\n")
