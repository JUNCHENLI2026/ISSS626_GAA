# Run from the exercise directory: Rscript data/download-data.R
# The same pinned Hunan inputs support both parts of Exercise 5.
commit <- "966eea6df498b38f10bc209d07a922f4d87b95d9"
base <- paste0("https://raw.githubusercontent.com/tskam/ISSS626-AY2026-27Aug/",
               commit, "/lesson/Lesson04/data/")
files <- c("aspatial/Hunan_2012.csv",
           paste0("geospatial/Hunan.", c("dbf", "prj", "qpj", "shp", "shx")))
dest <- file.path("data", files)
if (any(file.exists(dest))) stop("Input files already exist; refusing to overwrite them.")
dir.create("data/aspatial", recursive = TRUE, showWarnings = FALSE)
dir.create("data/geospatial", recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(files)) {
  download.file(paste0(base, files[i]), dest[i], mode = "wb")
  stopifnot(file.info(dest[i])$size > 0)
}
capture.output({
  cat("Hunan county boundaries and 2012 indicators\n")
  cat("Repository commit:", commit, "\nRetrieval date:", as.character(Sys.Date()), "\n")
  cat("Indicator year: 2012. Download date is not the indicator year.\n")
  cat(paste0(base, files, collapse = "\n"), "\n")
  print(tools::md5sum(dest))
}, file = "data/download-manifest.txt")
