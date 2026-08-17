# Hotelling analysis with a quadratic spline and the single knot 0

# Set the spline, test, and output options for this analysis.
Sys.setenv(
  CHAMBER_PAIR_FIXED_KNOTS = "0",
  CHAMBER_PAIR_EXCLUDE_IR2_RELATED = "true",
  CHAMBER_PAIR_HOTELLING_ALPHA = "0.05",
  CHAMBER_PAIR_OUTPUT_DIRECTORY = Sys.getenv(
    "CHAMBER_PAIR_OUTPUT_DIRECTORY_OVERRIDE",
    file.path(
      "visualization", "chamber_output", "quadratic_splines",
      "fixed_knots_0", "exclude_ir2_related", "alpha_0.05",
      "all_pairs_gamma_grid_0_3_step_0.01"
    )
  )
)

# Run the common all-pairs Chamber screening workflow.
source("script/chamber/chamber_quadratic_all_pairs_screening.R")
