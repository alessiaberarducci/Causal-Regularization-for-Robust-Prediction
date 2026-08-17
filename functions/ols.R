# Pooled OLS estimator

# Gplus and Zplus contain the moments from both environments.
compute_ols <- function(m) {
  qr.solve(m$Gplus, m$Zplus)
}
