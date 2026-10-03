# Table 1: Chamber classification metrics and Welch tests.

source("code/chamber/chamber_hotelling_quadratic_knot_0.R")

if (any(screening_results$status != "ok")) {
  stop("At least one Chamber pair failed; inspect per-pair results.")
}
positive <- screening_results[screening_results$gamma_CV > 0, , drop = FALSE]
zero <- screening_results[screening_results$gamma_CV == 0, , drop = FALSE]
if (!nrow(positive) || !nrow(zero)) stop("Both gamma subsets are required.")
metrics <- c(Sensitivity = "sensitivity", Specificity = "specificity",
             Precision = "precision", `F1 score` = "F1",
             `Balanced accuracy` = "balanced_accuracy", MCC = "MCC")
table_1 <- data.frame(
  Metric = names(metrics),
  OLS = vapply(metrics, function(m) mean(positive[[paste0(m, "_OLS")]]), 0),
  CR = vapply(metrics, function(m) mean(positive[[paste0(m, "_CV")]]), 0),
  OLS_equals_CR = vapply(metrics, function(m) mean(zero[[paste0(m, "_OLS")]]), 0),
  row.names = NULL
)
reported_tests <- do.call(rbind, lapply(c("sensitivity", "specificity"), function(m) {
  test <- t.test(positive[[paste0(m, "_CV")]],
                 positive[[paste0(m, "_OLS")]], paired = FALSE)
  data.frame(metric = m, method = test$method, statistic = unname(test$statistic),
             df = unname(test$parameter), p_value = test$p.value)
}))
write.csv(table_1, "output/table_1.csv", row.names = FALSE)
write.csv(reported_tests, "output/table_1_tests.csv", row.names = FALSE)
write.csv(data.frame(total_pairs = nrow(screening_results),
                     positive_gamma_pairs = nrow(positive), zero_gamma_pairs = nrow(zero)),
          "output/table_1_pair_counts.csv", row.names = FALSE)
print(table_1, digits = 7, row.names = FALSE)
print(reported_tests, row.names = FALSE)
