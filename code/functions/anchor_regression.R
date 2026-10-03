# No-intercept anchor regression with centered anchors.

fit_anchor_regression <- function(X, y, anchor, gamma) {
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("gamma must be one finite nonnegative number.")
  }
  X <- as.matrix(X)
  anchor <- as.matrix(anchor)
  if (nrow(X) != length(y) || nrow(anchor) != length(y)) {
    stop("X, y and anchor must have the same number of observations.")
  }
  Q <- anchor_column_basis(anchor)
  projected_X <- Q %*% crossprod(Q, X)
  projected_y <- drop(Q %*% crossprod(Q, y))
  weighted_X <- X + (sqrt(gamma) - 1) * projected_X
  weighted_y <- y + (sqrt(gamma) - 1) * projected_y
  drop(qr.solve(weighted_X, weighted_y))
}
