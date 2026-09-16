# Geographic extraction helper. Review and understand before using in coursework.
# Raw sources are never edited. Existing outputs are not overwritten.
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "E:/JUNCHENLI/ISSS626_GAA/Take-Home_Ex/Take-Home_Ex01/data"
out <- file.path(root, "Ketapang")
if (dir.exists(out)) stop("Output folder already exists; inspect it before rerunning.")
shp <- list.files(file.path(root, "Indonesia-geospasial"), pattern="[.]shp$", recursive=TRUE, full.names=TRUE)
stopifnot(length(shp) == 1)
layer <- tools::file_path_sans_ext(basename(shp))
sql <- sprintf("SELECT * FROM \"%s\" WHERE WADMPR = 'Kalimantan Barat' AND WADMKK = 'Ketapang'", layer)
cat("Reading Ketapang village polygons...\n")
villages_raw <- st_read(shp, query=sql, quiet=TRUE)
stopifnot(nrow(villages_raw) > 0, all(villages_raw$KDPKAB == "61.04", na.rm=TRUE))
# Use XY only. UTM 49S is a working CRS for topology operations, not a final
# selection of the projection/parameters for the student's statistical analysis.
villages <- st_transform(st_zm(villages_raw), 32749)
valid_before <- st_is_valid(villages)
cat("Village records:", nrow(villages), "; invalid before repair:", sum(!valid_before), "\n")
villages <- st_make_valid(villages)
stopifnot(all(st_is_valid(villages)), !any(st_is_empty(villages)))
window <- st_sf(province="Kalimantan Barat", regency="Ketapang", regency_code="61.04",
                geometry=st_union(st_geometry(villages)))
stopifnot(all(st_is_valid(window)), !any(st_is_empty(window)))

files <- list.files(file.path(root,"NASA FIRMS"), pattern="[.]csv$", full.names=TRUE)
stopifnot(length(files) == 2)
records <- lapply(files, function(f) {
  d <- read.csv(f, colClasses="character", check.names=FALSE, stringsAsFactors=FALSE)
  d$source_file <- basename(f)
  d$source_row <- seq_len(nrow(d))
  d$data_status <- if (grepl("fire_archive_",basename(f))) "standard" else "NRT"
  d
})
all <- bind_rows(records)
stopifnot(!anyNA(as.Date(all$acq_date)))
lon <- as.numeric(all$longitude); lat <- as.numeric(all$latitude)
stopifnot(all(is.finite(lon)), all(is.finite(lat)), all(abs(lon)<=180), all(abs(lat)<=90))
stopifnot(all(as.Date(all$acq_date)>=as.Date("2026-01-01")), all(as.Date(all$acq_date)<=as.Date("2026-08-31")))
points <- st_as_sf(all, coords=c("longitude","latitude"), crs=4326, remove=FALSE) |> st_transform(32749)
cat("Matching all fire detections to the derived Ketapang boundary...\n")
# Include points exactly on the supplied boundary; do not buffer or snap.
inside <- lengths(st_intersects(points, window)) > 0
selected <- all[inside,,drop=FALSE]
selected <- selected[order(selected$acq_date, selected$acq_time, selected$source_file, selected$source_row),]
strict <- lengths(st_within(points[inside,], window)) > 0
key <- c("latitude","longitude","acq_date","acq_time","satellite","instrument")
months <- data.frame(month=sprintf("2026-%02d",1:8))
monthly <- selected |> mutate(month=substr(acq_date,1,7)) |> count(month,data_status,name="detections")
monthly <- merge(expand.grid(month=months$month,data_status=c("standard","NRT")),monthly,all.x=TRUE)
monthly$detections[is.na(monthly$detections)] <- 0L

# Diagnose the unnamed West Kalimantan polygon without assigning it to Ketapang.
unknown_sql <- sprintf("SELECT * FROM \"%s\" WHERE WADMPR = 'Kalimantan Barat' AND (WADMKK IS NULL OR WADMKK = '')",layer)
unknown <- st_read(shp,query=unknown_sql,quiet=TRUE)
unknown_points <- 0L
unknown_overlap_km2 <- 0
if(nrow(unknown)) {
  unknown <- st_make_valid(st_transform(st_zm(unknown),32749))
  unknown_points <- sum(lengths(st_intersects(points,unknown))>0)
  overlap <- suppressWarnings(st_intersection(st_union(st_geometry(unknown)), st_geometry(window)))
  unknown_overlap_km2 <- sum(as.numeric(st_area(overlap)))/1e6
}
stopifnot(nrow(selected)>0)
dir.create(out,recursive=TRUE)
write.csv(selected,file.path(out,"ketapang_modis_2026_jan_aug.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(monthly,file.path(out,"monthly_record_counts.csv"),row.names=FALSE)
st_write(st_transform(window,4326),file.path(out,"ketapang_boundary.gpkg"),layer="study_window",quiet=TRUE)
st_write(st_transform(villages,4326),file.path(out,"ketapang_boundary.gpkg"),layer="villages",quiet=TRUE)
checks <- data.frame(metric=c("input_records","selected_records","excluded_records","village_records","invalid_villages_before_repair","invalid_villages_after_repair","boundary_only_records","duplicate_detection_keys","unnamed_west_kalimantan_polygons","records_in_unnamed_polygons","unnamed_polygon_overlap_km2"),
 value=c(nrow(all),nrow(selected),sum(!inside),nrow(villages),sum(!valid_before),sum(!st_is_valid(villages)),sum(!strict),sum(duplicated(selected[,key])),nrow(unknown),unknown_points,unknown_overlap_km2))
write.csv(checks,file.path(out,"filter_checks.csv"),row.names=FALSE)
inputs <- c(list.files(dirname(shp),full.names=TRUE),files)
hashes <- tools::md5sum(inputs)
write.csv(data.frame(file=basename(inputs),md5=unname(hashes)),file.path(out,"source_checksums.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"session-info.txt"))
# Round-trip checks, with text preserved so HHMM retains leading zeros.
saved <- read.csv(file.path(out,"ketapang_modis_2026_jan_aug.csv"),colClasses="character")
stopifnot(nrow(saved)==sum(inside),identical(saved$acq_time,selected$acq_time))
saved_window <- st_read(file.path(out,"ketapang_boundary.gpkg"),layer="study_window",quiet=TRUE) |> st_transform(32749)
stopifnot(all(lengths(st_intersects(points[inside,],saved_window))>0))
cat("\nCHECKS\n");print(checks,row.names=FALSE)
cat("\nMONTHLY\n");print(xtabs(detections~month+data_status,monthly))
cat("\nDATE RANGE\n");print(range(selected$acq_date))
cat("\nTYPE\n");print(table(selected$type,useNA="ifany"))
cat("\nOutput:",out,"\n")
