# Trace high reciprocal-intensity weights to source detections without deletion.
# Usage: Rscript weight_leverage_audit.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "weight_audit")
if (dir.exists(out)) stop("Audit output already exists; review before rerunning.")
estimate_file <- file.path(root, "inhomogeneous_checks", "intensity_reweighted_summaries.rds")
input_paths <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv",
  "ketapang_modis_2026_jan_aug.csv"))
raw_names <- c("fire_archive_M-C61_803220.csv", "fire_nrt_M-C61_803220.csv")
raw_paths <- file.path(dirname(root), "NASA FIRMS", raw_names)
all_paths <- c(input_paths, estimate_file, raw_paths)
hashes <- tools::md5sum(all_paths)
estimates <- readRDS(estimate_file)
stopifnot(identical(unname(hashes[1:3]), unname(estimates$source_md5)))
read_text <- function(path) read.csv(path, colClasses = "character", check.names = FALSE)
raw_tables <- setNames(lapply(raw_paths, read_text), raw_names)
sets <- lapply(input_paths[2:3], read_text)
all_conf <- read_text(input_paths[4])
boundary <- st_transform(st_read(input_paths[1], layer = "study_window", quiet = TRUE), 32749)
window <- rescale(as.owin(boundary), 1000, "km")
coordinates <- function(d) st_coordinates(st_transform(
  st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), 32749)) / 1000
key <- function(d) paste(d$source_file, d$source_row, sep = ":")
all_xy <- coordinates(all_conf)
all_keys <- key(all_conf)
stopifnot(!anyDuplicated(all_keys))
audit_rows <- list()
rank_curves <- list()
map_points <- list()
highlight_points <- list()
for (i in 1:2) {
  d <- sets[[i]]
  threshold <- c(30, 80)[i]
  label <- paste0("Confidence >=", threshold)
  xy <- coordinates(d)
  pattern <- ppp(xy[, 1], xy[, 2], window = window)
  stopifnot(npoints(pattern) == nrow(d), all(inside.owin(xy[, 1], xy[, 2], window)))
  border_distance <- bdist.points(pattern)
  # Verify every source column against its original source row, not just keys.
  for (filename in unique(d$source_file)) {
    positions <- which(d$source_file == filename)
    source <- raw_tables[[filename]]
    stopifnot(!is.null(source))
    rows <- as.integer(d$source_row[positions])
    stopifnot(all(rows >= 1 & rows <= nrow(source)))
    for (field in names(source)) {
      stopifnot(identical(d[[field]][positions], source[[field]][rows]))
    }
  }
  map_points[[i]] <- data.frame(x = xy[, 1], y = xy[, 2], subset = label)
  for (sigma in c(5, 10, 20)) {
    result <- estimates$estimates[[paste0("conf", threshold, "_sigma", sigma)]]
    stopifnot(identical(result$source_file, d$source_file),
              identical(result$source_row, d$source_row))
    lambda <- as.numeric(result$lambda_at_points)
    stopifnot(all(is.finite(lambda)), all(lambda > 0))
    shares <- (1 / lambda) / sum(1 / lambda)
    ranked <- order(shares, decreasing = TRUE)
    rank_curves[[length(rank_curves) + 1L]] <- data.frame(
      rank = 1:20, cumulative_percent = 100 * cumsum(shares[ranked])[1:20],
      subset = label, bandwidth = paste0(sigma, " km"))
    for (rank in 1:5) {
      j <- ranked[rank]
      distances <- sqrt((xy[, 1] - xy[j, 1])^2 + (xy[, 2] - xy[j, 2])^2)
      distances[j] <- Inf
      nearest <- which.min(distances)
      distances_all <- sqrt((all_xy[, 1] - xy[j, 1])^2 + (all_xy[, 2] - xy[j, 2])^2)
      distances_all[all_keys == key(d)[j]] <- Inf
      nearest_all <- which.min(distances_all)
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        subset = label, sigma_km = sigma, rank = rank,
        source_file = d$source_file[j], source_row = as.integer(d$source_row[j]),
        latitude = as.numeric(d$latitude[j]), longitude = as.numeric(d$longitude[j]),
        acq_date = d$acq_date[j], acq_time = d$acq_time[j],
        confidence = as.numeric(d$confidence[j]), satellite = d$satellite[j],
        data_status = d$data_status[j], supplied_type = d$type[j],
        intensity = lambda[j], reciprocal_share_percent = 100 * shares[j],
        distance_to_boundary_km = border_distance[j],
        nearest_retained_km = distances[nearest],
        retained_neighbours_within_10km = sum(distances <= 10),
        nearest_retained_date = d$acq_date[nearest],
        nearest_retained_day_gap = abs(as.integer(as.Date(d$acq_date[j]) - as.Date(d$acq_date[nearest]))),
        nearest_all_confidence_km = distances_all[nearest_all],
        nearest_all_confidence_value = as.numeric(all_conf$confidence[nearest_all]),
        source_fields_match = TRUE, inside_supplied_window = TRUE)
    }
    if (sigma == 5) {
      highlights <- ranked[1:5]
      highlight_points[[i]] <- data.frame(x = xy[highlights, 1], y = xy[highlights, 2],
                                         subset = label, rank = 1:5)
    }
  }
}
audit <- do.call(rbind, audit_rows)
cumulative <- do.call(rbind, rank_curves)
cumulative$bandwidth <- factor(cumulative$bandwidth, levels = c("5 km", "10 km", "20 km"))
all_points <- do.call(rbind, map_points)
highlights <- do.call(rbind, highlight_points)
shape <- st_coordinates(boundary)
outline <- data.frame(x = shape[, "X"] / 1000, y = shape[, "Y"] / 1000,
  ring = apply(shape[, grep("^L", colnames(shape)), drop = FALSE], 1, paste, collapse = "_"))
plot_theme <- theme_minimal(base_size = 12) + theme(
  plot.title = element_text(face = "bold", size = 16),
  strip.text = element_text(face = "bold"), legend.position = "bottom",
  plot.caption = element_text(hjust = 0, size = 9), plot.caption.position = "plot",
  plot.margin = margin(16, 16, 16, 16))
p_weights <- ggplot(cumulative, aes(rank, cumulative_percent, colour = bandwidth)) +
  geom_line(linewidth = 0.85) + geom_point(size = 1.3) +
  facet_wrap(~subset, nrow = 1) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Concentration of reciprocal-intensity weights",
    subtitle = "Ketapang MODIS | Highest-weight records ranked separately for each bandwidth and subset",
    x = "Number of highest-weight records", y = "Cumulative reciprocal-weight share (%)",
    colour = "Gaussian bandwidth",
    caption = paste("Weights are proportional to 1 / leave-one-out intensity. These are NOT shares of all K-function pair contributions.",
      "The audit does not delete records, select a final bandwidth, or test clustering. Sources: NASA FIRMS; supplied boundary.", sep = "\n")) + plot_theme
p_map <- ggplot() +
  geom_path(data = outline, aes(x, y, group = ring), colour = "#64748B", linewidth = 0.3) +
  geom_point(data = all_points, aes(x, y), colour = "#94A3B8", size = 0.35, alpha = 0.45) +
  geom_point(data = highlights, aes(x, y), shape = 21, fill = "#D55E00", colour = "white", size = 2.6) +
  facet_wrap(~subset, nrow = 1) + coord_equal() +
  labs(title = "Locations flagged for reciprocal-weight review",
    subtitle = "Orange circles: five highest-weight records in each subset at sigma = 5 km",
    x = "UTM zone 49S easting (km)", y = "UTM zone 49S northing (km)",
    caption = paste("All retained detections remain in grey. Flags identify numerical leverage candidates, not confirmed data errors.",
      "Source-row details and distances are recorded in checks.txt. Sources: NASA FIRMS; supplied Indonesia Geospatial boundary.", sep = "\n")) + plot_theme
stopifnot(identical(hashes, tools::md5sum(all_paths)), nrow(audit) == 30,
          all(audit$source_fields_match), all(audit$inside_supplied_window))
dir.create(out, recursive = TRUE)
ggsave(file.path(out, "reciprocal_weight_concentration.png"), p_weights,
       width = 12, height = 7, dpi = 160, bg = "white")
ggsave(file.path(out, "flagged_record_locations.png"), p_map,
       width = 12, height = 9, dpi = 160, bg = "white")
saveRDS(list(record_audit = audit, cumulative_weights = cumulative,
             source_md5 = hashes), file.path(out, "record_audit.rds"))
capture.output({
  cat("RECIPROCAL-WEIGHT AUDIT | ", as.character(Sys.Date()), "\n", sep = "")
  cat("All source columns of both confidence subsets match their original NASA source rows.\n")
  cat("All retained points are inside the supplied study window. No records were removed or changed.\n")
  cat("Source rows are one-based data rows, excluding the CSV header. Times retain HHMM formatting.\n")
  cat("Distances use EPSG:32749 in kilometres; full-period spatial nearest neighbours, not time-constrained neighbours.\n")
  cat("All-confidence nearest neighbours are searched only within the existing Ketapang geographic subset.\n")
  cat("These checks cannot independently establish fire type, coordinate accuracy or boundary currency.\n\n")
  for (label in unique(audit$subset)) {
    for (sigma in c(5, 10, 20)) {
      cat("\n", label, " | sigma = ", sigma, " km\n", sep = "")
      a <- audit[audit$subset == label & audit$sigma_km == sigma, ]
      print(a[, c("rank", "source_file", "source_row", "acq_date", "acq_time",
                  "latitude", "longitude", "confidence", "satellite", "data_status", "supplied_type")],
            row.names = FALSE)
      print(a[, c("rank", "intensity", "reciprocal_share_percent", "distance_to_boundary_km",
                  "nearest_retained_km", "retained_neighbours_within_10km",
                  "nearest_retained_date", "nearest_retained_day_gap",
                  "nearest_all_confidence_km", "nearest_all_confidence_value")], row.names = FALSE, digits = 6)
    }
  }
  cat("\nA large reciprocal share alone is not the same as measured influence on the entire K curve.\n")
  cat("For border correction, centre eligibility also depends on distance to the boundary at each r.\n")
  cat("Records can affect intensity normalization even if they have no nearby pairs in the evaluated range.\n")
  cat("No evidence-based exclusion rule or final bandwidth was selected by this audit.\n")
  cat("\nUnchanged input MD5:\n"); print(hashes)
  print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(audit[audit$rank == 1, c("subset", "sigma_km", "source_row", "acq_date", "confidence",
      "reciprocal_share_percent", "distance_to_boundary_km", "nearest_retained_km",
      "retained_neighbours_within_10km", "nearest_all_confidence_km")], row.names = FALSE, digits = 5)
cat("Completed. Original files and prior results unchanged. Outputs:", out, "\n")
