# Run from the exercise directory: Rscript data/download-data.R
# Downloads the instructor's pinned course-data files without modifying them.
commit <- "966eea6df498b38f10bc209d07a922f4d87b95d9"
base_url <- paste0("https://raw.githubusercontent.com/tskam/ISSS626-AY2026-27Aug/",
                   commit, "/lesson/Lesson04/data/")
files <- c("aspatial/Hunan_2012.csv", paste0("geospatial/Hunan.",
            c("dbf", "prj", "qpj", "shp", "shx")))
dir.create("data/aspatial", recursive = TRUE, showWarnings = FALSE)
dir.create("data/geospatial", recursive = TRUE, showWarnings = FALSE)
for (file in files) {
  destination <- file.path("data", file)
  if (file.exists(destination)) stop("Refusing to overwrite: ", destination)
  download.file(paste0(base_url, file), destination, mode = "wb", quiet = FALSE)
  stopifnot(file.info(destination)$size > 0)
}
capture.output({
  cat("Source: instructor's ISSS626-AY2026-27Aug repository\n")
  cat("Commit:", commit, "\nRetrieved:", as.character(Sys.Date()), "\n")
  cat("Indicator reference year: 2012 (not the download year).\n")
  cat("Raw file MD5 checksums:\n")
  print(tools::md5sum(file.path("data", files)))
  cat("\nSource URLs:\n", paste0(base_url, files, collapse = "\n"), "\n")
}, file = "data/download-manifest.txt")
