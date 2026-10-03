# Download the RNA and Causal Chamber datasets.

download_inputs <- function() {
  dir.create("data", showWarnings = FALSE)
  cache <- tempfile("causal_regularization_data_")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)

  rna_url <- paste0("https://raw.githubusercontent.com/xwshen51/DRIG/",
                    "278d3ea74ddd1f37302a06ed609c7169d8772564/data/dataset_rpe1.csv")
  rna_file <- file.path(cache, "dataset_rpe.csv")
  download.file(rna_url, rna_file, mode = "wb", method = "libcurl")
  if (unname(tools::md5sum(rna_file)) != "79ded4e226c6613d40d0999d38f46e8e") {
    stop("RNA download is incomplete or differs from the published dataset.")
  }
  if (!file.copy(rna_file, "data/dataset_rpe.csv", overwrite = TRUE)) {
    stop("Could not save RNA data.")
  }

  chamber_url <- paste0("https://causalchamber.s3.eu-central-1.amazonaws.com/",
                       "downloadables/lt_interventions_standard_v1.zip")
  archive <- file.path(cache, "chamber.zip")
  download.file(chamber_url, archive, mode = "wb", method = "libcurl")
  if (unname(tools::md5sum(archive)) != "476664d024f88e8b7640998bb5e9ee33") {
    stop("Chamber download is incomplete or differs from the published dataset.")
  }
  entries <- unzip(archive, list = TRUE)$Name
  entries <- entries[grepl("[.]csv$", entries)]
  if (length(entries) != 59L || anyDuplicated(basename(entries)) ||
      any(grepl("(^/|[.][.]/)", entries))) {
    stop("Unexpected Chamber archive contents.")
  }
  unzip(archive, files = entries, exdir = cache)
  target_dir <- "data/lt_interventions_standard_v1"
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  for (entry in entries) {
    if (!file.copy(file.path(cache, entry), file.path(target_dir, basename(entry)),
                   overwrite = TRUE)) stop("Could not save ", basename(entry))
  }
  cat("Downloaded RNA data and 59 Chamber CSV files.\n")
}
download_inputs()
