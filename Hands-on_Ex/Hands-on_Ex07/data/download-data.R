# Run from the Hands-on_Ex07 directory. Existing files are never overwritten.
base_url <- "https://raw.githubusercontent.com/tskam/IS415/acaa1f9995d3d778a02a09c6bbdd5ef4f4eae773/Hands-on_Ex/Hands-on_EX11/data/"
files <- c("aspatial/OD_Matrix.csv",
  paste0("geospatial/MP14_SUBZONE_NO_SEA_PL.", c("shp", "shx", "dbf", "prj")),
  paste0("geospatial/hexagons.", c("shp", "shx", "dbf", "prj", "cpg")),
  paste0("geospatial/ELDERCARE.", c("shp", "shx", "dbf", "prj", "cpg")))
for (f in files) {
  target <- file.path("data", f)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(target)) download.file(paste0(base_url, f), target, mode = "wb")
}
manifest <- data.frame(file = files, source = paste0(base_url, files),
                       md5 = unname(tools::md5sum(file.path("data", files))))
write.table(manifest, "data/download-manifest.txt", sep = "\t", row.names = FALSE)
