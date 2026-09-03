# 绘图层入口：与 organized.R 相同（linkET qcorrplot + geom_couple）。
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(script_dir, "organized.R"), encoding = "UTF-8")
