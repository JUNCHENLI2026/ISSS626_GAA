# Descriptive spatial recurrence separates many pixels from many recorded dates.
# Run: Rscript spatial_recurrence.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "spatial_recurrence")
if (dir.exists(out)) stop("Output directory exists; preserve it before a deliberate rerun.")
inputs <- file.path(root, c("ketapang_boundary.gpkg", "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
hashes <- tools::md5sum(inputs)
window <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
villages <- st_transform(st_read(inputs[1], layer = "villages", quiet = TRUE), 32749)
stopifnot(!anyNA(villages$WADMKC), all(nzchar(villages$WADMKC)))
districts <- villages |> group_by(WADMKC) |> summarise(do_union = TRUE, .groups = "drop")
districts$area_km2 <- as.numeric(st_area(districts)) / 1e6
tables <- summaries <- assignments <- admin_tables <- list()
settings <- expand.grid(cell_km = c(5, 10), shift_fraction = c(0, .5))
labels <- c("Confidence >=30", "Confidence >=80")
for (s in seq_along(labels)) {
  d <- read.csv(inputs[s + 1L], colClasses = "character")
  pt <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), 32749)
  xy <- st_coordinates(pt) / 1000
  hit <- st_intersects(pt, districts)
  stopifnot(all(lengths(hit) == 1L))
  d$district <- districts$WADMKC[unlist(hit)]
  d$date <- as.Date(d$acq_date)
  admin <- d |> group_by(district) |> summarise(records = n(), recorded_dates = n_distinct(date),
    august_records = sum(format(date, "%m") == "08"), .groups = "drop")
  admin <- merge(st_drop_geometry(districts)[, c("WADMKC", "area_km2")], admin,
    by.x = "WADMKC", by.y = "district", all.x = TRUE)
  for (col in c("records", "recorded_dates", "august_records")) admin[[col]][is.na(admin[[col]])] <- 0L
  admin$records_per_100_km2 <- 100 * admin$records / admin$area_km2
  admin$subset <- labels[s]
  stopifnot(sum(admin$records) == nrow(d))
  admin_tables[[s]] <- admin
  for (k in seq_len(nrow(settings))) {
    size <- settings$cell_km[k]; shift <- size * settings$shift_fraction[k]
    gx <- floor((xy[, 1] - shift) / size); gy <- floor((xy[, 2] - shift) / size)
    cell <- paste(gx, gy, sep = "_")
    rows <- data.frame(cell = cell, gx = gx, gy = gy, date = d$date, district = d$district)
    stats <- rows |> group_by(cell, gx, gy) |> summarise(records = n(), recorded_dates = n_distinct(date),
      recorded_months = n_distinct(format(date, "%Y-%m")),
      largest_daily_count = max(table(date)), .groups = "drop")
    stats$max_daily_fraction <- stats$largest_daily_count / stats$records
    dominant <- rows |> count(cell, district) |> arrange(cell, desc(n), district) |>
      group_by(cell) |> slice_head(n = 1) |> ungroup() |> select(cell, district)
    stats <- left_join(stats, dominant, by = "cell")
    stats$cell_km <- size; stats$shift_fraction <- settings$shift_fraction[k]; stats$subset <- labels[s]
    stats$rank_records <- rank(-stats$records, ties.method = "min")
    stats$rank_dates <- rank(-stats$recorded_dates, ties.method = "min")
    key <- paste(labels[s], size, settings$shift_fraction[k], sep = " | ")
    tables[[key]] <- stats
    assignments[[key]] <- data.frame(source_file = d$source_file, source_row = d$source_row, cell = cell)
    summaries[[key]] <- data.frame(subset = labels[s], cell_km = size, shift_fraction = settings$shift_fraction[k],
      occupied_cells = nrow(stats), cells_one_date = sum(stats$recorded_dates == 1),
      median_recorded_dates = median(stats$recorded_dates), max_recorded_dates = max(stats$recorded_dates),
      spearman_records_dates = cor(stats$records, stats$recorded_dates, method = "spearman"))
    stopifnot(sum(stats$records) == nrow(d), all(stats$recorded_dates <= stats$records))
  }
}
base <- tables[["Confidence >=30 | 5 | 0"]]
polygons <- lapply(seq_len(nrow(base)), function(i) {
  x <- base$gx[i] * 5000; y <- base$gy[i] * 5000
  st_polygon(list(matrix(c(x,y,x+5000,y,x+5000,y+5000,x,y+5000,x,y), ncol = 2, byrow = TRUE)))
})
cells <- st_sf(base, geometry = st_sfc(polygons, crs = 32749))
clipped <- suppressWarnings(st_intersection(cells, st_geometry(window)))
stopifnot(nrow(clipped) == nrow(base))
clipped$clipped_area_km2 <- as.numeric(st_area(clipped)) / 1e6
top <- base |> arrange(desc(recorded_dates), desc(records), cell) |> slice_head(n = 10)
high <- tables[["Confidence >=80 | 5 | 0"]]
top$high_confidence_dates <- high$recorded_dates[match(top$cell, high$cell)]
top$high_confidence_dates[is.na(top$high_confidence_dates)] <- 0L
top$cell_label <- paste0("E", top$gx * 5, "_N", top$gy * 5)
map_data <- rbind(transform(clipped, metric = "Recorded dates", value = recorded_dates),
  transform(clipped, metric = "Detection records", value = records))
map_data$metric <- factor(map_data$metric, levels = c("Detection records", "Recorded dates"))
# Separate panels use metric-specific scales by producing separate scientific maps.
make_map <- function(field, title, legend) ggplot() +
  geom_sf(data = window, fill = "#EEF2F4", colour = "#687681", linewidth = .25) +
  geom_sf(data = clipped, aes(fill = .data[[field]]), colour = NA) +
  geom_sf(data = districts, fill = NA, colour = "#687681", linewidth = .2) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = legend) +
  coord_sf(datum = st_crs(4326)) +
  labs(title = title, subtitle = "Ketapang, January-August 2026, confidence >=30, 5 km grid",
    caption = "Grey = no retained record. Grid origin: UTM 49S (0, 0). Border cells are clipped. Values are observed counts, not risk.") +
  theme_minimal(base_size = 12) + theme(plot.caption = element_text(size = 9, hjust = 0))
dir.create(out)
ggsave(file.path(out, "recorded_dates_map.png"), make_map("recorded_dates", "Where detections recur across dates", "Dates"), width = 9, height = 8, dpi = 160, bg = "white")
ggsave(file.path(out, "record_count_map.png"), make_map("records", "Where many detection records accumulate", "Records"), width = 9, height = 8, dpi = 160, bg = "white")
scatter <- ggplot(base, aes(records, recorded_dates, colour = max_daily_fraction)) +
  geom_point(size = 2, alpha = .75) + scale_x_log10() +
  scale_colour_viridis_c(option = "D", name = "Largest day's\nrecord share", limits = c(0, 1)) +
  labs(x = "Records per occupied 5 km cell (log scale)", y = "Distinct UTC dates with records",
    title = "Record volume and temporal recurrence answer different questions",
    subtitle = "Confidence >=30, January-August 2026",
    caption = "Each point is an occupied cell. No-record days have unknown observation opportunity. No significance test is attached to this plot.") +
  theme_minimal(base_size = 12) + theme(plot.caption = element_text(hjust = 0, size = 9))
ggsave(file.path(out, "volume_recurrence.png"), scatter, width = 11, height = 6, dpi = 160, bg = "white")
stopifnot(identical(hashes, tools::md5sum(inputs)))
summary <- do.call(rbind, summaries); admin <- do.call(rbind, admin_tables)
saveRDS(list(cells = tables, assignments = assignments, summary = summary, top = top,
  main_grid = clipped, districts = districts, administration = admin, source_md5 = hashes,
  settings = list(grid = settings, crs = 32749, temporal_unit = "UTC calendar date", period = "2026-01-01 through 2026-08-31")),
  file.path(out, "spatial_recurrence.rds"))
write.csv(top, file.path(out, "recurrent_cells.csv"), row.names = FALSE)
write.csv(admin, file.path(out, "subdistrict_context.csv"), row.names = FALSE)
capture.output({print(summary); print(top); print(sessionInfo())}, file = file.path(out, "checks.txt"))
print(summary); print(top[, c("cell_label", "district", "records", "recorded_dates", "high_confidence_dates", "max_daily_fraction")])
