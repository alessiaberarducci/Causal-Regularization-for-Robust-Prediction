# Figure 10: causal regularization and anchor regression paths.

source("code/functions/datautil.R")
source("code/functions/cd.R")
source("code/functions/ols.R")
source("code/functions/cross_valid_fun_n.R")
source("code/functions/confidence_intervals.R")
source("code/functions/anchor_infinity.R")
source("code/functions/anchor_regression.R")

set.seed(123L)
data <- simulate_example_d_data(1000L)
X <- rbind(data$Xe, data$Xo)
y <- c(data$ye, data$yo)
anchor <- matrix(c(rep(-0.5, nrow(data$Xe)),
                   rep(0.5, nrow(data$Xo))), ncol = 1L)
gamma_grid <- seq(0, 30, by = 0.1)
cr <- estimate_beta_path_ci(data, gamma_grid)
ar <- vapply(gamma_grid, function(gamma) {
  fit_anchor_regression(X, y, anchor, gamma)
}, numeric(2L))
ar_limit <- fit_anchor_infinity(X, y, anchor, intercept = FALSE)
output_directory <- "output/simulations"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
results <- data.frame(
  gamma = gamma_grid, beta_1_CR = cr[1, ], beta_2_CR = cr[2, ],
  beta_1_AR = ar[1, ], beta_2_AR = ar[2, ],
  beta_1_AR_infinity = ar_limit$coefficients[1],
  beta_2_AR_infinity = ar_limit$coefficients[2]
)
write.csv(results, file.path(output_directory, "figure_10_path.csv"),
          row.names = FALSE)
saveRDS(list(data = data, X = X, y = y, anchor = anchor, results = results,
             ar_limit = ar_limit),
        file.path(output_directory, "figure_10_results.rds"))

pdf(file.path(output_directory, "figure_10.pdf"), width = 7, height = 5,
    family = "serif", useDingbats = FALSE)
par(mar = c(8.0, 5.2, 1.2, 2.6), mgp = c(2.7, 0.8, 0),
    cex.axis = 1.1, cex.lab = 1.1, lend = "round")
matplot(gamma_grid, t(cr), type = "l", lty = 1, lwd = 2,
        col = "black", xlab = expression(gamma), ylab = expression(beta(gamma)),
        ylim = range(cr, ar, ar_limit$coefficients, c(1, 0)) + c(-0.08, 0.16),
        xlim = c(0, 30), axes = FALSE, xaxs = "i")
axis(1); axis(2, las = 1); box(bty = "l")
matlines(gamma_grid, t(ar), lty = 1, lwd = 2, col = "grey55")
abline(h = c(1, 0), lty = 2, col = "black", lwd = 1.5)
abline(h = ar_limit$coefficients, lty = 2, col = "grey55", lwd = 1.5)
legend("bottom", inset = c(0, -0.48), xpd = TRUE, bty = "n", ncol = 4,
       legend = c(expression(hat(beta)[1](gamma)~CR),
                  expression(hat(beta)[2](gamma)~CR),
                  expression(hat(beta)[1](gamma)~AR),
                  expression(hat(beta)[2](gamma)~AR),
                  expression(beta[1*",CP"]), expression(beta[2*",CP"]),
                  expression(beta[1*",AR,"*infinity]),
                  expression(beta[2*",AR,"*infinity])),
       col = c("black", "black", "grey55", "grey55",
               "black", "black", "grey55", "grey55"),
       lty = c(rep(1, 4), rep(2, 4)), lwd = 2, cex = 0.9)
dev.off()
cat("Anchor limit:", ar_limit$coefficients, "\n")
