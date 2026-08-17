# Repeat the subsample experiment with explicitly specified training groups.

first_training <- c(
  "ENSG00000147604", "ENSG00000133112", "ENSG00000108518",
  "ENSG00000125691", "ENSG00000006715", "ENSG00000079999"
)
second_training <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000110700",
  "ENSG00000172757", "ENSG00000067225", "ENSG00000012174",
  "ENSG00000149357"
)

Sys.setenv(
  RPE_FIRST_TRAINING = paste(first_training, collapse = ","),
  RPE_SECOND_TRAINING = paste(second_training, collapse = ","),
  RPE_OUTPUT_SUFFIX = "_alternative_n3",
  RPE_SAMPLE_SIZES = "850,650,500,250,150,100"
)
source("script/RNA/rpe_Cgamma_subsample_comparison.R")
hist((current$OLS-current$CV))
