# Causal Dantzig estimator

# Solve the empirical moment equation for beta.
compute_cd <- function(m) {
  qr.solve(m$Gdelta, m$Zdelta)
}
