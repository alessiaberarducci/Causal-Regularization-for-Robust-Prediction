# Empirical moments

# Build the sum and difference moments used by all estimators.
moments <- function(data) {
  n_e <- nrow(data$Xe)
  n_o <- nrow(data$Xo)
  XYe <- crossprod(data$Xe, data$ye) / n_e
  XXe <- crossprod(data$Xe) / n_e
  XYo <- crossprod(data$Xo, data$yo) / n_o
  XXo <- crossprod(data$Xo) / n_o

  # Keep the environment-specific moments for later comparisons.
  list(
    XYe = XYe,
    XXe = XXe,
    ne = n_e,
    XYo = XYo,
    XXo = XXo,
    no = n_o,
    Zdelta = XYe - XYo,
    Gdelta = XXe - XXo,
    Gplus = XXe + XXo,
    Zplus = XYe + XYo,
    p = ncol(data$Xe)
  )
}
