# Figures 11 and 12: RNA prediction experiments.

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
source("code/RNA/rpe_Cgamma_subsample_comparison.R")

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
    family = "serif", mar = c(5.0, 5.0, 2.7, 0.8), lend = "round",
    cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.10
  )
  on.exit(graphics::par(old_par), add = TRUE)
  current_data <- result$data
  class_counts <- table(current_data$class)
  graphics::plot(
    current_data$OLS, current_data$CV,
    type = "n", xlim = selected_scatter_limits,
    ylim = selected_scatter_limits, asp = 1,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CR])),
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
      bquote(C[gamma] * "\\" * C[0] * ": " * .(class_counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(class_counts[["C0"]]))
    )),
    col = scatter_colours, pch = 16,
    pt.cex = 1.0, bty = "n", cex = 1.0
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
    family = "serif", mar = c(5.0, 5.0, 2.7, 0.8), lend = "round",
    cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.10
  )
  on.exit(graphics::par(old_par), add = TRUE)
  current_data <- result$data
  class_counts <- table(current_data$class)
  graphics::plot(
    current_data$OLS, current_data$CV,
    type = "n", log = "xy",
    xlim = selected_log_limits, ylim = selected_log_limits, asp = 1,
    xlab = expression(MSE[e](hat(beta)[OLS])),
    ylab = expression(MSE[e](hat(beta)[CR])),
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
      bquote(C[gamma] * "\\" * C[0] * ": " * .(class_counts[["Cgamma"]])),
      bquote(C[0] * ": " * .(class_counts[["C0"]]))
    )),
    col = scatter_colours, pch = 16,
    pt.cex = 1.0, bty = "n", cex = 1.0
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

write_mse_difference_histogram <- function(result, file) {
  current_data <- result$data
  mse_difference <- current_data$OLS - current_data$CV
  complete <- is.finite(mse_difference) & !is.na(current_data$class)
  current_data <- current_data[complete, , drop = FALSE]
  mse_difference <- mse_difference[complete]

  class_order <- c("outside", "Cgamma", "C0")
  histogram_colours <- c(
    outside = "#009E73",
    Cgamma = "#E69F00",
    C0 = "#0072B2"
  )
  data_range <- range(mse_difference)
  data_padding <- max(0.001, 0.02 * diff(data_range))
  plot_limits <- data_range + c(-1, 1) * data_padding
  common_breaks <- seq(plot_limits[1], plot_limits[2], length.out = 81L)
  class_counts <- vapply(class_order, function(current_class) {
    graphics::hist(
      mse_difference[current_data$class == current_class],
      breaks = common_breaks, plot = FALSE
    )$counts
  }, numeric(80L))
  total_counts <- rowSums(class_counts)

  grDevices::pdf(
    file, width = 6.3, height = 5.2,
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  old_par <- graphics::par(
    family = "serif", mar = c(5.0, 5.0, 1.2, 0.8), bty = "l",
    cex.axis = 1.15, cex.lab = 1.20
  )
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot(
    NA, NA, type = "n",
    xlim = range(common_breaks), ylim = c(0, max(total_counts) * 1.08),
    xlab = expression(
      MSE[e](hat(beta)[OLS]) - MSE[e](hat(beta)[CR])
    ),
    ylab = "Frequency", main = "", bty = "l", xaxt = "n"
  )
  cumulative_counts <- rep(0, length(total_counts))
  for (current_class in class_order) {
    current_counts <- class_counts[, current_class]
    graphics::rect(
      xleft = common_breaks[-length(common_breaks)],
      ybottom = cumulative_counts,
      xright = common_breaks[-1L],
      ytop = cumulative_counts + current_counts,
      col = histogram_colours[[current_class]],
      border = "grey25", lwd = 0.25
    )
    cumulative_counts <- cumulative_counts + current_counts
  }
  graphics::abline(v = 0, lty = 2, lwd = 2.2, col = "black")
  x_ticks <- sort(unique(c(pretty(plot_limits), 0)))
  graphics::axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE))
  graphics::box(bty = "l")
  graphics::legend(
    "topright",
    legend = as.expression(list(
      bquote("Outside " * C[gamma] * ": " * .(
        sum(current_data$class == "outside")
      )),
      bquote(C[gamma] * "\\" * C[0] * ": " * .(
        sum(current_data$class == "Cgamma")
      )),
      bquote(C[0] * ": " * .(sum(current_data$class == "C0")))
    )),
    fill = histogram_colours[class_order],
    border = histogram_colours[class_order],
    bty = "n"
  )
}

full_sample_difference_histogram_file <- file.path(
  rna_output_directory,
  "rpe_OLS_minus_CR_MSE_histogram_full_sample_alternative_n2.pdf"
)
n100_difference_histogram_file <- file.path(
  rna_output_directory,
  "rpe_OLS_minus_CR_MSE_histogram_n100_alternative_n2.pdf"
)
write_mse_difference_histogram(
  selected_scatter_results$full_sample,
  full_sample_difference_histogram_file
)
write_mse_difference_histogram(
  selected_scatter_results$n100,
  n100_difference_histogram_file
)
