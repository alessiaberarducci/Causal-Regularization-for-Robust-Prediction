# Pooled ordinary least squares.

compute_ols <- function(m) {
  qr.solve(m$Gplus, m$Zplus)
}
