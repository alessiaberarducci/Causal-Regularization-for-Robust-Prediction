# Causal Dantzig estimation.

compute_cd <- function(m) {
  qr.solve(m$Gdelta, m$Zdelta)
}
