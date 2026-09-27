# Exploratory time-window comparison, not a new primary study period or test.
# Usage: Rscript time_window_comparison.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "time_window_checks")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
prior_file <- file.path(root, "temporal_checks", "temporal_diagnostics.rds")
prior_hash <- tools::md5sum(prior_file)
prior <- readRDS(prior_file)
stopifnot(identical(unname(before), unname(prior$source_md5)))
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
sigma_km <- 10
eps_km <- 0.5
periods <- c("January-August", "January-July", "August")
days <- c(243L, 212L, 31L)
labels <- c("Confidence >=30", "Confidence >=80")
results <- list()
count_frames <- comparison_frames <- map_frames <- list()
panel_levels <- character()
for (s in 1:2) {
  d <- read.csv(inputs[s + 1L], colClasses = "character")
  date <- as.Date(d$acq_date)
  stopifnot(nrow(d) == c(4472L, 1726L)[s], !anyNA(date),
    all(date >= as.Date("2026-01-01") & date <= as.Date("2026-08-31")))
  pts <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), 32749)
  xy <- st_coordinates(pts) / 1000
  stopifnot(all(inside.owin(xy[, 1], xy[, 2], window)))
  in_august <- date >= as.Date("2026-08-01")
  masks <- list(rep(TRUE, nrow(d)), !in_august, in_august)
  n <- vapply(masks, sum, integer(1))
  stopifnot(n[1] == n[2] + n[3], all(n > 0),
    n[3] == prior$overview$august_records[prior$overview$subset == labels[s]])
  estimates <- list()
  for (k in seq_along(periods)) {
    keep <- masks[[k]]
    x <- ppp(xy[keep, 1], xy[keep, 2], window = window)
    cat("Computing", labels[s], periods[k], "with", npoints(x), "records.\n")
    flush.console()
    # Same fixed smoothing and full geographic window across all periods.
    # Pixel intensities include every retained point; leave-one-out is inapplicable.
    z <- density(x, sigma = sigma_km, kernel = "gaussian", edge = TRUE,
      diggle = TRUE, at = "pixels", eps = eps_km)
    valid <- !is.na(z$v)
    stopifnot(all(is.finite(z$v[valid])), min(z$v[valid]) >= -1e-8)
    area_pixel <- z$xstep * z$ystep
    raw_integral <- integral.im(z)
    stopifnot(abs(raw_integral - n[k]) / n[k] < 0.005)
    # Preserve raw estimates. Clip only tiny negative FFT round-off before
    # normalising a non-negative descriptive spatial mass distribution.
    cleaned <- z
    cleaned$v[valid] <- pmax(cleaned$v[valid], 0)
    clipped_mass <- sum(pmax(-z$v[valid], 0)) * area_pixel
    cleaned_integral <- integral.im(cleaned)
    relative <- cleaned / cleaned_integral
    mass <- relative$v[valid] * area_pixel
    stopifnot(abs(sum(mass) - 1) < 1e-10)
    estimates[[periods[k]]] <- list(raw_intensity = z, relative_density = relative,
      pixel_mass = mass, source_file = d$source_file[keep], source_row = d$source_row[keep])
    count_frames[[length(count_frames) + 1L]] <- data.frame(
      subset = labels[s], period = periods[k], calendar_days = days[k], records = n[k],
      percent_of_full_period = 100 * n[k] / n[1], raw_integral = raw_integral,
      integral_relative_error = abs(raw_integral - n[k]) / n[k],
      negative_roundoff_cells = sum(z$v[valid] < 0), clipped_roundoff_mass = clipped_mass,
      normalized_integral = integral.im(relative))
    grid <- as.data.frame(relative)
    names(grid)[1:3] <- c("x", "y", "relative_density")
    grid <- grid[!is.na(grid$relative_density), ]
    grid$pixel_width <- relative$xstep
    grid$pixel_height <- relative$ystep
    grid$relative_density_scaled <- grid$relative_density * 1e4
    grid$panel <- paste0(labels[s], " | ", periods[k], "\nn = ", format(n[k], big.mark = ",", trim = TRUE))
    panel_levels <- c(panel_levels, grid$panel[1])
    map_frames[[length(map_frames) + 1L]] <- grid
  }
  full <- estimates[[1]]
  early <- estimates[[2]]
  aug <- estimates[[3]]
  for (v in list(early$raw_intensity, aug$raw_intensity)) {
    stopifnot(identical(full$raw_intensity$xcol, v$xcol),
      identical(full$raw_intensity$yrow, v$yrow),
      identical(is.na(full$raw_intensity$v), is.na(v$v)))
  }
  # Direct full-period KDE must agree with the sum of the disjoint-period KDEs.
  linear_error <- max(abs(full$raw_intensity$v - early$raw_intensity$v - aug$raw_intensity$v), na.rm = TRUE)
  stopifnot(linear_error < 1e-8)
  w <- n[2] / n[1]
  mixture_error <- max(abs(full$pixel_mass - (w * early$pixel_mass + (1 - w) * aug$pixel_mass)))
  stopifnot(mixture_error < 1e-9)
  tv_full_aug <- 0.5 * sum(abs(full$pixel_mass - aug$pixel_mass))
  tv_early_aug <- 0.5 * sum(abs(early$pixel_mass - aug$pixel_mass))
  # Algebraic nesting check: full-period similarity to August is partly automatic.
  stopifnot(abs(tv_full_aug - w * tv_early_aug) < 1e-7,
    tv_full_aug <= w + 1e-7, tv_early_aug >= 0, tv_early_aug <= 1)
  comparison_frames[[length(comparison_frames) + 1L]] <- data.frame(
    subset = labels[s], early_records = n[2], august_records = n[3],
    early_fraction = w, full_vs_august_TV = tv_full_aug,
    early_vs_august_TV = tv_early_aug,
    full_vs_august_overlap = sum(pmin(full$pixel_mass, aug$pixel_mass)),
    early_vs_august_overlap = sum(pmin(early$pixel_mass, aug$pixel_mass)),
    full_vs_august_TV_upper_bound = w,
    raw_linearity_max_error = linear_error, normalized_mixture_max_error = mixture_error)
  results[[labels[s]]] <- estimates
}
counts <- do.call(rbind, count_frames)
comparison <- do.call(rbind, comparison_frames)
map_data <- do.call(rbind, map_frames)
map_data$panel <- factor(map_data$panel, levels = panel_levels)
coords <- st_coordinates(boundary)
outline <- data.frame(x = coords[, "X"] / 1000, y = coords[, "Y"] / 1000,
  ring = apply(coords[, grep("^L", colnames(coords)), drop = FALSE], 1, paste, collapse = "_"))
p <- ggplot(map_data, aes(x, y, fill = relative_density_scaled)) +
  geom_tile(aes(width = pixel_width, height = pixel_height)) +
  geom_path(data = outline, aes(x, y, group = ring), inherit.aes = FALSE, colour = "#94A3B8", linewidth = 0.2) +
  facet_wrap(~panel, ncol = 3) + coord_equal(expand = FALSE) +
  scale_fill_viridis_c(option = "magma", limits = c(0, max(map_data$relative_density_scaled)),
    name = "Relative density\n(10^-4 per km2)") +
  labs(title = "Time-window sensitivity of relative spatial distribution",
    subtitle = "Ketapang MODIS | Same illustrative 10 km bandwidth and shared colour scale | UTC dates",
    x = "UTM zone 49S easting (km)", y = "UTM zone 49S northing (km)",
    caption = paste("Each surface is normalised to integrate to one: compare spatial shape, NOT counts, rates, fire risk or burned area.",
      "January-August contains August. January-July and August are disjoint but are not assumed statistically independent.",
      "Gaussian KDE; Jones-Diggle edge correction; 0.5 km computational grid. No final bandwidth or new primary time window selected.",
      "January-July includes both standard and NRT records; August is NRT. Repeated detections and observation opportunity remain unadjusted.", sep = "\n")) +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10), strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
    plot.caption.position = "plot", plot.margin = margin(18, 18, 18, 18))
stopifnot(identical(before, tools::md5sum(inputs)), identical(prior_hash, tools::md5sum(prior_file)))
dir.create(out)
ggsave(file.path(out, "relative_spatial_distribution.png"), p, width = 14, height = 11, dpi = 160, bg = "white")
saveRDS(list(counts = counts, comparison = comparison, estimates = results,
  source_md5 = before, prior_md5 = prior_hash,
  settings = list(sigma_km = sigma_km, eps_km = eps_km, projection = "EPSG:32749",
    timezone = "UTC", periods = periods, coordinate_units = "km",
    density_edge = "Jones-Diggle", normalisation = "unit integral of nonnegative grid density",
    negative_roundoff_tolerance = 1e-8, final_parameter_selection = FALSE)),
  file.path(out, "time_window_diagnostics.rds"))
capture.output({
  cat("EXPLORATORY TIME-WINDOW COMPARISON: no new primary window or significance test.\n\n")
  print(counts, row.names = FALSE, digits = 7)
  cat("\nNORMALISED GRID COMPARISONS\n"); print(comparison, row.names = FALSE, digits = 7)
  cat("\nAll six estimates use the full polygon window, Gaussian sigma 10 km, Jones-Diggle correction and eps 0.5 km.\n")
  cat("Sigma 10 km is a shared illustrative setting from previous diagnostics, not an automatic or final choice.\n")
  cat("Raster estimates include all retained points; leave-one-out applies only to point evaluation, not these maps.\n")
  cat("Raw intensities are preserved. Tiny negative FFT round-off is clipped only before nonnegative mass normalization.\n")
  cat("Relative surfaces integrate to one. They remove record-count differences; they are not fire-risk maps.\n")
  cat("Total variation (TV) = 0.5 * sum(abs(p - q)) for unit-sum pixel masses.\n")
  cat("Overlap = sum(pmin(p,q)) = 1 - TV. Neither is a p-value or a percentage of matching records.\n")
  cat("Spatial shape similarity depends on this bandwidth and grid; this check does not establish robustness at other scales.\n")
  cat("The full-period surface is a count-weighted mixture of January-July and August at the same settings.\n")
  cat("TV(full, August) = January-July record fraction * TV(January-July, August), up to numerical error.\n")
  cat("Thus similarity of full-period and August surfaces is partly algebraic, not independent validation.\n")
  cat("January-July versus August avoids shared records, but temporal dependence and unequal sample sizes remain.\n")
  cat("No temporal thinning, record removal, coordinate jitter, bandwidth refitting, K-function fitting or simulations.\n")
  cat("Calendar days are duration descriptors, not valid-observation or cloud-free exposure denominators.\n")
  cat("All dates use UTC; August was inspected after observing temporal concentration and is a sensitivity window.\n")
  cat("Sources contain a standard/NRT transition. Repeated detections and observation opportunity remain unadjusted.\n")
  cat("Checks: geographic inclusion, record reconciliation, mass integrals, grid alignment, linearity and mixture identity.\n")
  cat("\nSource MD5 unchanged:\n"); print(before)
  cat("\nPrevious temporal results unchanged:\n"); print(prior_hash)
  cat("\nMethod reference: https://search.r-project.org/CRAN/refmans/spatstat.explore/html/density.ppp.html\n")
  cat("\nRuntime:\n"); print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(counts[, 1:5], row.names = FALSE)
print(comparison, row.names = FALSE)
cat("Saved time-window outputs to", out, "\n")
