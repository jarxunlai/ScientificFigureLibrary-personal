# Generate the bundled synthetic inputs if they are absent, then run the plot.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE)) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
needed <- file.path(root, "data", c("input_deg_heatmap.csv", "input_cell_annotation.csv", "input_go_enrichment.csv"))
if (!all(file.exists(needed))) {
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(root)
  source(file.path("code", "generate-example-data.R"), encoding = "UTF-8")
}
source(file.path(script_dir, "organized.R"), encoding = "UTF-8")
