# Rebuild the analysis from original downloads into a new directory.
# Rscript reproduce.R <exercise_directory> <new_output_directory> <quarto_executable>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply exercise directory, new output directory and Quarto executable.")
exercise <- normalizePath(args[1], mustWork = TRUE)
destination <- args[2]
if (dir.exists(destination) || file.exists(destination)) stop("Choose a new output directory. Existing outputs are never removed.")
quarto <- normalizePath(args[3], mustWork = TRUE)
scripts <- file.path(exercise, "data/Ketapang")
required <- c("sf", "dplyr", "ggplot2", "spatstat.geom", "spatstat.explore", "knitr", "rmarkdown")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
dir.create(destination, recursive = TRUE)
destination <- normalizePath(destination, mustWork = TRUE)
dir.create(file.path(destination, "data"))
dir.create(file.path(destination, "data/NASA FIRMS"))
nasa_inputs <- list.files(file.path(exercise, "data/NASA FIRMS"), pattern = "[.]csv$", full.names = TRUE)
stopifnot(length(nasa_inputs) == 2L, all(file.copy(nasa_inputs, file.path(destination, "data/NASA FIRMS"))))
stopifnot(identical(unname(tools::md5sum(nasa_inputs)), unname(tools::md5sum(file.path(destination, "data/NASA FIRMS", basename(nasa_inputs))))))
dir.create(file.path(destination, "logs"))
root <- file.path(destination, "data/Ketapang")
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, arguments, label) {
  cat("RUNNING:", label, "\n"); flush.console()
  result <- system2(rscript, args = c(shQuote(script), vapply(arguments, shQuote, character(1))),
    stdout = file.path(destination, "logs", paste0(label, ".log")), stderr = "")
  if (!identical(result, 0L)) stop("Step failed: ", label, ". Inspect logs; no existing project outputs were changed.")
}
run(file.path(scripts, "filter_ketapang.R"), c(file.path(exercise, "data"), root), "01-extract")
d <- read.csv(file.path(root, "ketapang_modis_2026_jan_aug.csv"), colClasses = "character", check.names = FALSE)
confidence <- as.numeric(d$confidence)
stopifnot(!anyNA(confidence), nrow(d) == 4738L)
for (threshold in c(30, 80)) {
  filtered <- d[confidence >= threshold, ]
  original <- read.csv(file.path(scripts, paste0("ketapang_modis_2026_conf", threshold, ".csv")), colClasses = "character", check.names = FALSE)
  rownames(filtered) <- NULL; rownames(original) <- NULL
  stopifnot(identical(filtered, original))
  write.csv(filtered, file.path(root, paste0("ketapang_modis_2026_conf", threshold, ".csv")), row.names = FALSE, na = "")
}
months <- sprintf("2026-%02d", 1:8)
monthly <- do.call(rbind, lapply(months, function(m) {
  hit <- substr(d$acq_date, 1, 7) == m
  data.frame(month = m, all_detections = sum(hit), confidence_ge30 = sum(hit & confidence >= 30),
    confidence_ge80 = sum(hit & confidence >= 80), excluded_below30 = sum(hit & confidence < 30),
    excluded_below80 = sum(hit & confidence < 80),
    ge30_standard = sum(hit & confidence >= 30 & d$data_status == "standard"),
    ge30_nrt = sum(hit & confidence >= 30 & d$data_status == "NRT"),
    ge80_standard = sum(hit & confidence >= 80 & d$data_status == "standard"),
    ge80_nrt = sum(hit & confidence >= 80 & d$data_status == "NRT"))
}))
write.csv(monthly, file.path(root, "monthly_confidence_comparison.csv"), row.names = FALSE)
file.copy(list.files(scripts, pattern = "[.]R$", full.names = TRUE), root)
dir.create(file.path(root, "learning"))
file.copy(file.path(scripts, "learning/kde_learning.R"), file.path(root, "learning"))
sequence <- c("plot_data_checks.R", "learning/kde_learning.R", "second_order_diagnostics.R",
  "inhomogeneous_diagnostics.R", "weight_leverage_audit.R", "bandwidth_selector_checks.R",
  "temporal_repeat_checks.R", "time_window_comparison.R", "time_window_second_order.R",
  "conditional_permutation.R", "day_block_permutation.R", "spatial_recurrence.R",
  "robustness_extension.R", "verify_revision.R", "verify_robustness.R")
for (i in seq_along(sequence)) run(file.path(root, sequence[i]), root, sprintf("%02d-%s", i + 1, basename(sequence[i])))
stopifnot(all(file.copy(file.path(exercise, c("Take-Home_Ex01.qmd", "executive_summary.qmd", "summary.css", "report.css", "reproduce.R", "recurrence_map.R", "robustness_figures.R")), destination)))
file.copy(file.path(scripts, "README.md"), root)
writeLines(c("project:", "  type: default", "  render:", "    - Take-Home_Ex01.qmd", "    - executive_summary.qmd"), file.path(destination, "_quarto.yml"))
# Give Quarto an output-local temporary directory on Windows. This avoids
# a native-encoding failure when the account's default TEMP path is non-ASCII.
quarto_temp <- file.path(destination, "logs", "quarto-temp")
dir.create(quarto_temp, showWarnings = FALSE)
Sys.setenv(TMP = quarto_temp, TEMP = quarto_temp, TMPDIR = quarto_temp)
cat("RENDERING REPORT AND SUMMARY\n"); flush.console()
status <- system2(quarto, c("render", shQuote(destination)), stdout = file.path(destination, "logs", "render.log"), stderr = "")
if (!identical(status, 0L)) stop("Quarto render failed; inspect render.log.")
stopifnot(all(file.exists(file.path(destination, c("Take-Home_Ex01.html", "executive_summary.html")))))
capture.output(sessionInfo(), file = file.path(destination, "logs", "session-info.txt"))
cat("COMPLETE:", destination, "\n")
