# Risk measurements

# Squared error for each observation.
risk_vector <- function(X, y, beta) {
  drop(y - X %*% beta)^2
}

risk <- function(X, y, beta) {
  mean(risk_vector(X, y, beta))
}

in_sample_risk <- function(data, beta) {
  # Rplus is the sum of the risks in the two environments.
  risk(data$Xe, data$ye, beta) + risk(data$Xo, data$yo, beta)
}

difference <- function(data, beta) {
  risk(data$Xe, data$ye, beta) - risk(data$Xo, data$yo, beta)
}

abs_difference <- function(data, beta) {
  abs(difference(data, beta))
}

measure <- function(data, beta, metric) {
  # Apply the selected metric to every row of a coefficient matrix.
  apply(beta, 1L, function(b) metric(data, matrix(b, ncol = 1L)))
}
