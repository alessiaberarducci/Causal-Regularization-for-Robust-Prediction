# Chamber quadratic splines with a single knot at zero.

Sys.setenv(
  CHAMBER_PAIR_FIXED_KNOTS = "0",
  CHAMBER_PAIR_EXCLUDE_IR2_RELATED = "true",
  CHAMBER_PAIR_HOTELLING_ALPHA = "0.05",
  CHAMBER_PAIR_OUTPUT_DIRECTORY = Sys.getenv(
    "CHAMBER_PAIR_OUTPUT_DIRECTORY_OVERRIDE",
    file.path(
      "output", "chamber_output", "quadratic_splines",
      "fixed_knots_0", "exclude_ir2_related", "alpha_0.05",
      "all_pairs_gamma_grid_0_3_step_0.01"
    )
  )
)

source("code/chamber/chamber_quadratic_all_pairs_screening.R")
