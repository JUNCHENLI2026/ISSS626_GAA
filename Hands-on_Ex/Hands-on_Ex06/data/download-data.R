# Run from the exercise directory: Rscript data/download-data.R
commit <- "acaa1f9995d3d778a02a09c6bbdd5ef4f4eae773"
base <- paste0("https://raw.githubusercontent.com/tskam/IS415/", commit,
               "/Hands-on_Ex/Hands-on_Ex08/data/")
files <- c("aspatial/Shan-ICT.csv", paste0("geospatial/myanmar_township_boundaries.",
                                         c("cst", "dbf", "prj", "shp", "shx")))
paths <- file.path("data", files)
if (any(file.exists(paths))) stop("Existing raw files found; refusing overwrite.")
dir.create("data/aspatial", recursive = TRUE, showWarnings = FALSE)
dir.create("data/geospatial", recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(paths)) {
  download.file(paste0(base, files[i]), paths[i], mode = "wb")
  stopifnot(file.info(paths[i])$size > 0)
}
capture.output({
  cat("Source: instructor tskam/IS415, pinned commit", commit, "\n")
  cat("Retrieved:", as.character(Sys.Date()), "\nCensus indicator year: 2014\n")
  cat("Original provenance described by Chapter 12: MIMU.\n")
  cat("Boundary effective vintage is not independently established by this download.\n")
  cat(paste0(base, files, collapse = "\n"), "\n")
  print(tools::md5sum(paths))
}, file = "data/download-manifest.txt")
