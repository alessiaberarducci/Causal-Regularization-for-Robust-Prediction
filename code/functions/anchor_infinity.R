# The limiting anchor regression solution as gamma tends to infinity.

symmetric_pinv_solve <- function(matrix, rhs, tolerance = 1e-10) {
  matrix <- 0.5 * (matrix + t(matrix))
  rhs <- as.matrix(rhs)
  decomposition <- eigen(matrix, symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(decomposition$values)))
  keep <- abs(decomposition$values) > cutoff

  if (!any(keep)) {
    return(matrix(0, nrow = nrow(matrix), ncol = ncol(rhs)))
  }

  vectors <- decomposition$vectors[, keep, drop = FALSE]
  vectors %*% sweep(
    crossprod(vectors, rhs), 1L, decomposition$values[keep], "/"
  )
}

anchor_column_basis <- function(anchor, tolerance = 1e-10) {
  decomposition <- svd(anchor, nu = min(dim(anchor)), nv = 0)
  if (!length(decomposition$d)) {
    return(matrix(0, nrow = nrow(anchor), ncol = 0L))
  }

  cutoff <- tolerance * max(1, max(decomposition$d))
  keep <- decomposition$d > cutoff
  decomposition$u[, keep, drop = FALSE]
}

fit_anchor_infinity <- function(X, y, anchor, intercept, tolerance = 1e-10) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  anchor <- as.matrix(anchor)

  if (nrow(X) != length(y) || nrow(anchor) != length(y)) {
    stop("X, y, and anchor must contain the same number of observations.")
  }

  if (intercept) {
    X_work <- sweep(X, 2L, colMeans(X), "-")
    y_work <- y - mean(y)
    anchor_work <- sweep(anchor, 2L, colMeans(anchor), "-")
  } else {
    X_work <- X
    y_work <- y
    anchor_work <- anchor
  }

  Q <- anchor_column_basis(anchor_work, tolerance)
  projected_X <- crossprod(Q, X_work)
  projected_y <- crossprod(Q, y_work)

  projected_gram <- crossprod(projected_X)
  decomposition <- eigen(
    0.5 * (projected_gram + t(projected_gram)), symmetric = TRUE
  )
  cutoff <- tolerance * max(1, max(abs(decomposition$values)))
  identified <- abs(decomposition$values) > cutoff
  beta_first <- drop(symmetric_pinv_solve(
    projected_gram, crossprod(projected_X, projected_y), tolerance
  ))

  null_basis <- decomposition$vectors[, !identified, drop = FALSE]
  beta <- beta_first
  if (ncol(null_basis)) {
    X_null <- X_work %*% null_basis
    eta <- drop(symmetric_pinv_solve(
      crossprod(X_null),
      crossprod(X_null, y_work - X_work %*% beta_first),
      tolerance
    ))
    beta <- beta_first + drop(null_basis %*% eta)
  }

  residual <- y_work - drop(X_work %*% beta)
  projected_residual <- drop(crossprod(Q, residual))
  projected_residual_squared <- sum(projected_residual^2)
  zero_tolerance <- tolerance * max(1, sum(y_work^2))

  list(
    coefficients = beta,
    intercept = if (intercept) {
      mean(y) - drop(colMeans(X) %*% beta)
    } else {
      0
    },
    projected_residual_squared = projected_residual_squared,
    projected_residual_zero = projected_residual_squared <= zero_tolerance,
    anchor_rank = ncol(Q)
  )
}
