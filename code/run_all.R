# Reproduce all computational experiments.

arguments <- commandArgs(trailingOnly = TRUE)
scope <- if (length(arguments)) arguments[1] else "all"
if (!scope %in% c("all", "simulations", "rna", "chamber")) {
  stop("Usage: Rscript code/run_all.R [all|simulations|rna|chamber]")
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument)) {
  script_path <- normalizePath(sub("^--file=", "", script_argument[1]))
  setwd(dirname(dirname(script_path)))
}
if (!file.exists("code/functions/datautil.R")) {
  stop("Run from the repository root or invoke code/run_all.R by its path.")
}
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)
inherited <- names(Sys.getenv())
clear <- inherited[grepl("^(RPE_|CHAMBER_PAIR_)", inherited)]
clear <- setdiff(clear, "CHAMBER_PAIR_WORKERS")
Sys.unsetenv(clear)
scripts <- character(0)
if (scope %in% c("all", "simulations")) {
  scripts <- c(scripts, paste0("code/simulations/example_", LETTERS[1:4], ".R"),
               "code/simulations/anchor_comparison.R")
}
if (scope %in% c("all", "rna")) scripts <- c(scripts, "code/RNA/paper_RNA.R")
if (scope %in% c("all", "chamber")) scripts <- c(scripts, "code/chamber/paper_chamber.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (script in scripts) {
  log <- file.path("output/logs", paste0(tools::file_path_sans_ext(basename(script)), ".log"))
  cat("Running", script, "...\n")
  status <- system2(rscript, c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  if (status != 0L) stop(script, " failed. See ", log)
}
cat("Completed:", scope, "\n")
