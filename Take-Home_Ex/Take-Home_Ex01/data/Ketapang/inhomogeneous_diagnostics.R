# Exploratory intensity-reweighted summaries; not a fitted final null model.
# Usage: Rscript inhomogeneous_diagnostics.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "inhomogeneous_checks")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
baseline_file <- file.path(root, "second_order", "observed_summaries.rds")
baseline_hash <- tools::md5sum(baseline_file)
baseline <- readRDS(baseline_file)
stopifnot(identical(unname(before), unname(baseline$source_md5)),
          baseline$projection == "EPSG:32749", baseline$units == "km")
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
area_km2 <- area.owin(window)
stopifnot(abs(area_km2 - sum(as.numeric(st_area(boundary))) / 1e6) < 1e-4)
sigmas <- c(5, 10, 20)
expected <- c(4472L, 1726L)
panel_labels <- c("Confidence >=30 | n = 4,472", "Confidence >=80 | n = 1,726")
all_results <- list()
curve_frames <- list()
audit_frames <- list()
for (i in 1:2) {
  records <- read.csv(inputs[i + 1L], colClasses = "character")
  point_sf <- st_transform(st_as_sf(records, coords = c("longitude", "latitude"),
                                  crs = 4326), 32749)
  xy <- st_coordinates(point_sf) / 1000
  stopifnot(nrow(xy) == expected[i], all(inside.owin(xy[, 1], xy[, 2], window)))
  x <- ppp(xy[, 1], xy[, 2], window = window, checkdup = TRUE)
  stopifnot(npoints(x) == expected[i])
  previous <- baseline$results[[i]]
  stopifnot(previous$records == npoints(x))
  r <- previous$curve$r_km
  stopifnot(min(r) == 0, max(r) <= 10 + 1e-8)
  baseline_curve <- previous$curve
  baseline_curve$panel <- panel_labels[i]
  baseline_curve$method <- "Ordinary L (unadjusted)"
  curve_frames[[length(curve_frames) + 1L]] <-
    baseline_curve[, c("r_km", "L_minus_r_km", "panel", "method")]
  for (sigma in sigmas) {
    cat("Computing", panel_labels[i], "with sigma =", sigma, "km\n")
    flush.console()
    # Evaluate intensity at the original points. Exclude each point's own
    # kernel contribution; do not look up the earlier full-data raster KDE.
    lambda <- density(x, sigma = sigma, at = "points", leaveoneout = TRUE,
      kernel = "gaussian", edge = TRUE, diggle = TRUE, eps = 0.5)
    stopifnot(length(lambda) == npoints(x), all(is.finite(lambda)), all(lambda > 0))
    # No arbitrary floor, clipping or deletion is applied to intensity weights.
    # The existing 0-10 km grid and border correction are held constant.
    ki <- Kinhom(x, lambda = as.numeric(lambda), r = r, correction = "border",
                 renormalise = TRUE, normpower = 2)
    kd <- as.data.frame(ki)
    cat("Raw K range:", range(kd$border), "\n")
    # Cumulative numerical sums can produce tiny negative round-off near zero.
    # Keep raw K in the saved object; bound only the square-root display input.
    roundoff_tolerance_km2 <- 1e-8
    stopifnot(identical(kd$r, r), all(is.finite(kd$border)),
              all(kd$border >= -roundoff_tolerance_km2),
              max(abs(kd$theo - pi * r^2)) < 1e-8)
    method_label <- paste0("Reweighted L: sigma = ", sigma, " km")
    curve_frames[[length(curve_frames) + 1L]] <- data.frame(
      r_km = r, L_minus_r_km = sqrt(pmax(kd$border, 0) / pi) - r,
      panel = panel_labels[i], method = method_label)
    reciprocal <- 1 / as.numeric(lambda)
    shares <- reciprocal / sum(reciprocal)
    normalization <- area_km2 / sum(reciprocal)
    audit_frames[[length(audit_frames) + 1L]] <- data.frame(
      subset = panel_labels[i], sigma_km = sigma, records = npoints(x),
      minimum_intensity = min(lambda), median_intensity = median(lambda),
      maximum_intensity = max(lambda),
      maximum_reciprocal_share = max(shares),
      reciprocal_effective_count = 1 / sum(shares^2),
      normalization_c = normalization,
      correction_multiplier_c_squared = normalization^2,
      minimum_raw_K = min(kd$border),
      negative_roundoff_values = sum(kd$border < 0))
    key <- paste0("conf", c(30, 80)[i], "_sigma", sigma)
    all_results[[key]] <- list(subset = panel_labels[i], sigma_km = sigma,
      source_row = records$source_row, source_file = records$source_file,
      lambda_at_points = lambda, Kinhom = ki)
  }
}
curves <- do.call(rbind, curve_frames)
audit <- do.call(rbind, audit_frames)
method_levels <- c("Ordinary L (unadjusted)", paste0("Reweighted L: sigma = ", sigmas, " km"))
curves$method <- factor(curves$method, levels = method_levels)
curves$panel <- factor(curves$panel, levels = panel_labels)
colours <- setNames(c("#6B7280", "#0072B2", "#D55E00", "#009E73"), method_levels)
line_types <- setNames(c("dashed", "solid", "solid", "solid"), method_levels)
plot_curves <- function(d, title, subtitle) {
  ggplot(d, aes(r_km, L_minus_r_km, colour = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "#334155", linewidth = 0.4) +
    geom_line(linewidth = 0.8) + facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = colours, drop = TRUE) +
    scale_linetype_manual(values = line_types, drop = TRUE) +
    scale_x_continuous(breaks = seq(0, 10, 2)) +
    labs(title = title, subtitle = subtitle, x = "Distance r (km)",
      y = "L(r) - r (km)", colour = NULL, linetype = NULL,
      caption = paste(
        "Gaussian leave-one-out intensity at points; Jones-Diggle density correction; border-corrected Kinhom.",
        "Kinhom renormalise = TRUE, normpower = 2. Zero is a Poisson theoretical reference, not a significance threshold.",
        "Bandwidths and 0-10 km range are illustrative. Same data estimate intensity and the summary; no simulations or p-values.",
        "Nested confidence subsets; no adjustment for repeated detections or observation opportunity. Sources: NASA FIRMS; supplied boundary.",
        sep = "\n")) +
    theme_minimal(base_size = 12) + theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom", legend.text = element_text(size = 10),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.caption.position = "plot", plot.margin = margin(16, 16, 16, 16)) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
}
p_compare <- plot_curves(curves, "Exploratory L summaries: intensity reweighting",
  "Ketapang MODIS | January-August 2026 | Shared axes across confidence subsets")
p_sensitivity <- plot_curves(subset(curves, method != method_levels[1]),
  "Bandwidth sensitivity of intensity-reweighted L summaries",
  "Same curves as the comparison plot, with a closer vertical scale | Exploratory only")
stopifnot(identical(before, tools::md5sum(inputs)),
          identical(baseline_hash, tools::md5sum(baseline_file)))
dir.create(out, recursive = TRUE)
ggsave(file.path(out, "ordinary_vs_reweighted_L.png"), p_compare,
       width = 13, height = 8, dpi = 160, bg = "white")
ggsave(file.path(out, "reweighted_bandwidth_sensitivity.png"), p_sensitivity,
       width = 13, height = 8, dpi = 160, bg = "white")
saveRDS(list(estimates = all_results, curves = curves, intensity_audit = audit,
  source_md5 = before, baseline_md5 = baseline_hash, projection = "EPSG:32749",
  coordinate_units = "km", window_area_km2 = area_km2,
  settings = list(leaveoneout = TRUE, density_edge = "Jones-Diggle",
    density_eps_km = 0.5, K_edge = "border", renormalise = TRUE, normpower = 2)),
  file.path(out, "intensity_reweighted_summaries.rds"))
capture.output({
  cat("EXPLORATORY SENSITIVITY CHECKS: no final null model or significance test.\n")
  cat("Window area (km2):", area_km2, "\nProjection: EPSG:32749. Coordinates: km.\n")
  cat("Intensity: Gaussian, sigma 5/10/20 km, unit-weight detections, at points, leave-one-out.\n")
  cat("Density edge correction: Jones-Diggle; computational edge grid eps = 0.5 km.\n")
  cat("Kinhom: supplied positive intensity vector; border; renormalise TRUE; normpower 2.\n")
  cat("Distance grid matches previous 0-10 km diagnostic. K transformed as sqrt(K/pi).\n")
  cat("All original points retained, finite positive intensity estimates; K nonnegative to 1e-8 km2 tolerance.\n")
  cat("Raw K preserved; only negative round-off within tolerance is set to zero for the square-root plot.\n")
  cat("No intensity flooring, clipping, thinning or event merging applied.\n\n")
  print(audit, row.names = FALSE, digits = 6)
  cat("\nReciprocal effective count = 1/sum(normalised inverse-intensity weights squared).\n")
  cat("This is a weight-concentration diagnostic, NOT the number of independent observations.\n")
  cat("Small estimated intensities can cause high leverage; finite values alone do not establish reliability.\n")
  cat("Leave-one-out removes only the focal contribution; it is not independent validation.\n")
  cat("Intensity and summaries use the same observations. Bandwidth can absorb or retain spatial structure.\n")
  cat("Interpretation additionally requires an appropriate intensity-reweighted second-order model.\n")
  cat("Fitted intensity does not separately identify ecological causes and observation effects.\n")
  cat("Nested subsets are not independent. Full-period spatial summaries do not model time.\n")
  cat("No envelopes, p-values, significant-cluster claims or final parameter choices.\n")
  cat("A future fitted-null simulation must account for intensity estimation, not hold it fixed without justification.\n")
  cat("\nInput MD5 unchanged:\n"); print(before)
  cat("\nPrevious summary MD5 unchanged:\n"); print(baseline_hash)
  cat("\nReferences:\nhttps://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kinhom.html\n")
  cat("https://search.r-project.org/CRAN/refmans/spatstat.explore/html/density.ppp.html\n")
  print(sessionInfo())
}, file = file.path(out, "checks.txt"))
stopifnot(identical(before, tools::md5sum(inputs)),
          identical(baseline_hash, tools::md5sum(baseline_file)))
print(audit, row.names = FALSE, digits = 4)
cat("Completed. Source and baseline checksums unchanged. Outputs:", out, "\n")
