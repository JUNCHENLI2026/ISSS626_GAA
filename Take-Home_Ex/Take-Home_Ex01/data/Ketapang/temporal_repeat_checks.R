# Exploratory space-time proximity checks. No event merging or record deletion.
# Usage: Rscript temporal_repeat_checks.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "temporal_checks")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
previous_file <- file.path(root, "bandwidth_checks", "bandwidth_diagnostics.rds")
previous_hash <- tools::md5sum(previous_file)
previous <- readRDS(previous_file)
stopifnot(identical(unname(before), unname(previous$source_md5)))
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
end_exclusive <- as.POSIXct("2026-09-01 00:00:00", tz = "UTC")
dates <- seq(as.Date("2026-01-01"), as.Date("2026-08-31"), by = "day")
spatial_limits <- c(0.5, 1, 2)
time_limits <- c(6, 24, 72)
rules <- c("Any strictly earlier record", "At least 1 hour earlier")
labels <- c("Confidence >=30", "Confidence >=80")
results <- list()
summary_frames <- sensitivity_frames <- daily_frames <- acquisition_frames <- list()

for (s in 1:2) {
  d <- read.csv(inputs[s + 1L], colClasses = "character")
  stopifnot(nrow(d) == c(4472L, 1726L)[s], all(d$instrument == "MODIS"),
    all(d$satellite %in% c("Terra", "Aqua")), all(grepl("^[0-9]{4}$", d$acq_time)))
  hh <- as.integer(substr(d$acq_time, 1, 2))
  mm <- as.integer(substr(d$acq_time, 3, 4))
  stopifnot(all(hh %in% 0:23), all(mm %in% 0:59))
  t <- as.POSIXct(paste(d$acq_date, d$acq_time), format = "%Y-%m-%d %H%M", tz = "UTC")
  stopifnot(!anyNA(t), all(t >= start & t < end_exclusive),
    identical(format(t, "%Y-%m-%d", tz = "UTC"), d$acq_date),
    identical(format(t, "%H%M", tz = "UTC"), d$acq_time))
  spatial <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), 32749)
  xy <- st_coordinates(spatial) / 1000
  stopifnot(all(inside.owin(xy[, 1], xy[, 2], window)))
  x <- ppp(xy[, 1], xy[, 2], window = window)
  near <- closepairs(x, rmax = max(spatial_limits), twice = FALSE, what = "ijd")
  pairs <- data.frame(i = near$i, j = near$j, distance_km = near$d)
  stopifnot(all(pairs$i != pairs$j), all(pairs$distance_km <= 2 + 1e-10))
  pair_key <- paste(pmin(pairs$i, pairs$j), pmax(pairs$i, pairs$j), sep = ":")
  stopifnot(!anyDuplicated(pair_key))
  seconds <- as.numeric(t)
  delta <- seconds[pairs$j] - seconds[pairs$i]
  pairs$gap_hours <- abs(delta) / 3600
  pairs$same_satellite <- d$satellite[pairs$i] == d$satellite[pairs$j]
  # Simultaneous records have no strictly earlier/later endpoint.
  pairs$later <- ifelse(delta > 0, pairs$j, ifelse(delta < 0, pairs$i, NA_integer_))
  pairs$earlier <- ifelse(delta > 0, pairs$i, ifelse(delta < 0, pairs$j, NA_integer_))
  exact_key <- paste(d$latitude, d$longitude, d$acq_date, d$acq_time, d$satellite, sep = "|")
  # Date/minute/satellite is only a grouping proxy, NOT an orbit/granule identifier.
  acquisition_key <- paste(d$acq_date, d$acq_time, d$satellite, sep = "|")
  acquisition_counts <- table(acquisition_key)
  for (r in spatial_limits) {
    within <- pairs$distance_km <= r
    acquisition_frames[[length(acquisition_frames) + 1L]] <- data.frame(
      subset = labels[s], distance_km = r,
      all_time_pairs = sum(within),
      same_minute_same_satellite_pairs = sum(within & pairs$gap_hours == 0 & pairs$same_satellite),
      same_minute_different_satellite_pairs = sum(within & pairs$gap_hours == 0 & !pairs$same_satellite),
      positive_gap_under_1h_same_satellite_pairs = sum(within & pairs$gap_hours > 0 & pairs$gap_hours < 1 & pairs$same_satellite),
      positive_gap_under_1h_different_satellite_pairs = sum(within & pairs$gap_hours > 0 & pairs$gap_hours < 1 & !pairs$same_satellite))
  }
  flags <- list()
  for (rule in rules) for (r in spatial_limits) for (h in time_limits) {
    lower_ok <- if (rule == rules[1]) pairs$gap_hours > 0 else pairs$gap_hours >= 1
    keep <- pairs$distance_km <= r & lower_ok & pairs$gap_hours <= h
    q <- pairs[keep, ]
    has_prior <- seq_len(nrow(d)) %in% q$later
    involved <- unique(c(q$i, q$j))
    flag_key <- paste(rule, r, h, sep = "|")
    flags[[flag_key]] <- has_prior
    sensitivity_frames[[length(sensitivity_frames) + 1L]] <- data.frame(
      subset = labels[s], rule = rule, distance_km = r, max_gap_hours = h,
      records = nrow(d), pair_count = nrow(q), records_with_prior = sum(has_prior),
      percent_with_prior = 100 * mean(has_prior), records_in_any_pair = length(involved),
      same_satellite_pairs = sum(q$same_satellite), different_satellite_pairs = sum(!q$same_satellite),
      different_utc_date_pairs = sum(d$acq_date[q$i] != d$acq_date[q$j]))
    # Direct-distance/time verification for deterministic sample records,
    # independent of the closepairs enumeration and its endpoint ordering.
    check_ids <- unique(round(seq(1, nrow(d), length.out = 25)))
    for (id in check_ids) {
      distances <- sqrt((xy[, 1] - xy[id, 1])^2 + (xy[, 2] - xy[id, 2])^2)
      age <- (seconds[id] - seconds) / 3600
      lower_direct <- if (rule == rules[1]) age > 0 else age >= 1
      expected <- any(distances <= r & lower_direct & age <= h)
      stopifnot(identical(has_prior[id], expected))
    }
  }
  for (sat in c("Terra", "Aqua")) {
    counts <- as.integer(table(factor(d$acq_date[d$satellite == sat], levels = as.character(dates))))
    daily_frames[[length(daily_frames) + 1L]] <- data.frame(date = dates, satellite = sat,
      detections = counts, subset = labels[s])
  }
  date_counts <- table(factor(d$acq_date, levels = as.character(dates)))
  top <- which.max(date_counts)
  summary_frames[[length(summary_frames) + 1L]] <- data.frame(
    subset = labels[s], records = nrow(d),
    duplicate_source_identifiers = sum(duplicated(paste(d$source_file, d$source_row))),
    duplicate_detection_keys = sum(duplicated(exact_key)),
    duplicate_numeric_coordinates = sum(duplicated(xy)),
    dates_with_records = sum(date_counts > 0), dates_without_records = sum(date_counts == 0),
    acquisition_group_proxies = length(acquisition_counts),
    maximum_records_per_acquisition_proxy = max(acquisition_counts),
    peak_date_utc = names(date_counts)[top], peak_day_records = as.integer(date_counts[top]),
    peak_day_percent = 100 * max(date_counts) / nrow(d),
    august_records = sum(substr(d$acq_date, 1, 7) == "2026-08"),
    august_percent = 100 * mean(substr(d$acq_date, 1, 7) == "2026-08"),
    minimum_scan_km = min(as.numeric(d$scan)), maximum_scan_km = max(as.numeric(d$scan)),
    minimum_track_km = min(as.numeric(d$track)), maximum_track_km = max(as.numeric(d$track)))
  provenance <- d[, c("source_file", "source_row", "latitude", "longitude", "acq_date",
    "acq_time", "satellite", "confidence", "scan", "track", "data_status")]
  provenance$timestamp_utc <- format(t, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  results[[labels[s]]] <- list(records = provenance, xy_km = xy, pairs_within_2km = pairs,
    prior_flags = flags, acquisition_counts = acquisition_counts)
  cat("Completed", labels[s], "with", nrow(pairs), "spatial pairs within 2 km.\n")
}
overview <- do.call(rbind, summary_frames)
sensitivity <- do.call(rbind, sensitivity_frames)
daily <- do.call(rbind, daily_frames)
acquisitions <- do.call(rbind, acquisition_frames)
stopifnot(nrow(sensitivity) == 36L, nrow(daily) == length(dates) * 4L,
  identical(before, tools::md5sum(inputs)), identical(previous_hash, tools::md5sum(previous_file)))
for (lab in labels) {
  stopifnot(sum(daily$detections[daily$subset == lab]) == overview$records[overview$subset == lab])
  for (rule in rules) for (r in spatial_limits) {
    # Use explicit indexing to avoid data-mask ambiguity for the rule variable.
    z <- sensitivity[sensitivity$subset == lab & sensitivity$rule == rule & sensitivity$distance_km == r, ]
    stopifnot(all(diff(z$records_with_prior[order(z$max_gap_hours)]) >= 0))
  }
  for (rule in rules) for (h in time_limits) {
    z <- sensitivity[sensitivity$subset == lab & sensitivity$rule == rule & sensitivity$max_gap_hours == h, ]
    stopifnot(all(diff(z$records_with_prior[order(z$distance_km)]) >= 0))
  }
}
for (i in seq_len(18L)) {
  # Higher-confidence pairs form a subset when identifiers are mapped back.
  hi <- sensitivity[19:36, ][i, ]
  lo <- sensitivity[1:18, ][i, ]
  stopifnot(hi$pair_count <= lo$pair_count, hi$records_with_prior <= lo$records_with_prior)
}
dir.create(out)
saveRDS(list(overview = overview, sensitivity = sensitivity, daily = daily,
  acquisition_proximity = acquisitions, subsets = results, source_md5 = before,
  prior_md5 = previous_hash, settings = list(projection = "EPSG:32749", units = "km",
    timezone = "UTC", start = start, end_exclusive = end_exclusive,
    distance_km = spatial_limits, max_gap_hours = time_limits, rules = rules,
    maximum_pair_distance_km = 2, no_event_merging = TRUE)),
  file.path(out, "temporal_diagnostics.rds"))
p_daily <- ggplot(daily, aes(date, detections, fill = satellite)) +
  geom_col(width = 1) + facet_wrap(~subset, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = as.Date("2026-05-01"), linetype = "dashed", colour = "#64748B") +
  scale_fill_manual(values = c(Aqua = "#0072B2", Terra = "#D55E00")) +
  scale_x_date(breaks = as.Date(paste0("2026-", sprintf("%02d", 1:8), "-01")), date_labels = "%b",
    expand = expansion(add = c(1, 1))) +
  labs(title = "Recorded MODIS detections over time", subtitle = "Ketapang | January-August 2026 | UTC dates; different vertical scales",
    x = "Acquisition date (UTC)", y = "Recorded detections per day", fill = "Satellite",
    caption = paste("Dashed line: supplied files change from standard processing to NRT on 1 May; not an estimated change point.",
      "Zero means no retained record in the supplied data, not no fire or a confirmed clear-sky observation.",
      "Nested confidence subsets; counts are detections, not independent fires. Source: supplied NASA FIRMS MODIS files.", sep = "\n")) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9), plot.margin = margin(15, 15, 15, 15))
ggsave(file.path(out, "daily_detection_timing.png"), p_daily, width = 12, height = 8, dpi = 160, bg = "white")
sensitivity$rule <- factor(sensitivity$rule, levels = rules)
p_heat <- ggplot(subset(sensitivity, rule == "At least 1 hour earlier"),
  aes(factor(distance_km), factor(max_gap_hours), fill = percent_with_prior)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%\n%d / %d", percent_with_prior, records_with_prior, records)), size = 3.5) +
  facet_wrap(~subset, nrow = 1) + scale_fill_gradient(low = "#F1F5F9", high = "#5DCAB5", limits = c(0, 100)) +
  labs(title = "Records with at least one earlier nearby detection",
    subtitle = "At least 1 hour earlier | Exploratory proximity sensitivity, not a duplicate-fire classifier",
    x = "Maximum centre-to-centre distance (km)", y = "Maximum elapsed time (hours)", fill = "Records (%)",
    caption = paste("Denominator = all retained records in that confidence subset; each focal record is counted at most once per cell.",
      "1 hour <= gap <= time limit. Same-minute pairs are excluded; the alternative >0-hour rule is tabulated in checks.txt.",
      "Thresholds are illustrative; pixel footprints vary. The 1-hour gap does not prove separate satellite passes.",
      "No records are merged or removed. Prior detections outside Ketapang or before the supplied period are unavailable.", sep = "\n")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"),
    panel.grid = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(20, 15, 15, 15))
ggsave(file.path(out, "space_time_proximity_sensitivity.png"), p_heat, width = 12, height = 6.8, dpi = 160, bg = "white")
capture.output({
  cat("EXPLORATORY TEMPORAL-REPEAT CHECKS: proximity is not fire-event identity.\n\n")
  print(overview, row.names = FALSE, digits = 7)
  cat("\nACQUISITION-TIME PROXIMITY (unordered distinct record pairs):\n")
  print(acquisitions, row.names = FALSE)
  cat("\nPRIOR-RECORD SENSITIVITY:\n")
  print(sensitivity, row.names = FALSE, digits = 6)
  cat("\nUTC HHMM acquisition timestamps round-trip checked; coordinates projected to EPSG:32749, km.\n")
  cat("Exact detection key: latitude, longitude, acquisition date/minute and satellite. Not a fire-event ID.\n")
  cat("Date/minute/satellite groups are acquisition proxies, not verified orbit or granule identifiers.\n")
  cat("Distance limits 0.5/1/2 km and time limits 6/24/72 hours are sensitivity choices, not calibrated event rules.\n")
  cat("Both satellites are included. Same-minute pairs and positive gaps below 1 hour are reported separately.\n")
  cat("A one-hour gap excludes near-simultaneous detections but cannot guarantee separate passes.\n")
  cat("MODIS coordinates represent pixel centres; scan/track dimensions vary, and distances are not footprint-overlap tests.\n")
  cat("Same fire can persist or spread; separate fires can occur nearby. Neither identity is established here.\n")
  cat("Counts are detection records, not burned area, unique fires, or measures corrected for observation opportunity.\n")
  cat("No-record dates do not establish absence of fire or availability of valid clear-sky observations.\n")
  cat("Temporal and geographic truncation may hide prior records; no prior match is not evidence of independence.\n")
  cat("Higher-confidence threshold can remove earlier nearby records. Subsets are nested, not independent.\n")
  cat("No null simulations, p-values, event merging, thinning, new bandwidth choice or significance claims.\n")
  cat("Checks: unique unordered pair keys; direct checks on 25 focal rows x 18 settings per subset;\n")
  cat("monotone counts as spatial/time limits increase; nested absolute counts; reconciled daily totals.\n")
  cat("\nSource MD5 unchanged:\n"); print(before)
  cat("\nPrevious bandwidth results unchanged:\n"); print(previous_hash)
  cat("\nReference: https://firms.modaps.eosdis.nasa.gov/descriptions/FIRMS_MODIS_Firehotspots.html\n")
  cat("\nRuntime:\n"); print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(overview, row.names = FALSE)
print(sensitivity[sensitivity$rule == rules[2] & sensitivity$distance_km == 1 & sensitivity$max_gap_hours == 24, ], row.names = FALSE)
cat("Saved temporal diagnostics to", out, "\n")
