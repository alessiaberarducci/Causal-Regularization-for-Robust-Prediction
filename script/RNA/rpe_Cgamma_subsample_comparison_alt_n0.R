# Run the subsample experiment with the explicitly specified original groups.

first_training <- c(
  "ENSG00000147604", "ENSG00000110700",
  "ENSG00000133112", "ENSG00000125691"
)
second_training <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000172757",
  "ENSG00000067225", "ENSG00000108518"
)

Sys.setenv(
  RPE_FIRST_TRAINING = paste(first_training, collapse = ","),
  RPE_SECOND_TRAINING = paste(second_training, collapse = ","),
  RPE_OUTPUT_SUFFIX = "",
  RPE_SAMPLE_SIZES = "650,500,250,150,100"
)
source("script/RNA/rpe_Cgamma_subsample_comparison.R")
