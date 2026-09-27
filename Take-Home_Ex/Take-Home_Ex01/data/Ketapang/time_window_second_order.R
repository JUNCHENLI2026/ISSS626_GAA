# Exploratory spatial summaries across time windows; no fitted-null test.
# Usage: Rscript time_window_second_order.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "time_window_second_order")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg", "ketapang_modis_2026_conf30.csv",
  "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
previous_files <- file.path(root, c("second_order/observed_summaries.rds",
  "inhomogeneous_checks/intensity_reweighted_summaries.rds",
  "time_window_checks/time_window_diagnostics.rds"))
previous_hash <- tools::md5sum(previous_files)
ordinary_prior <- readRDS(previous_files[1])
weighted_prior <- readRDS(previous_files[2])
window_prior <- readRDS(previous_files[3])
for (obj in list(ordinary_prior, weighted_prior, window_prior)) {
  stopifnot(identical(unname(before), unname(obj$source_md5)))
}
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
area_km2 <- area.owin(window)
periods <- c("January-August", "January-July", "August")
labels <- c("Confidence >=30", "Confidence >=80")
sigmas <- c(10, 20)
r <- ordinary_prior$results[[1]]$curve$r_km
stopifnot(min(r) == 0, max(r) <= 10 + 1e-8, all(diff(r) > 0))
estimates <- list()
ordinary_frames <- weighted_frames <- audit_frames <- support_frames <- list()
to_L <- function(k) {
  # Preserve raw K. Only tiny negative round-off may be bounded for sqrt.
  valid <- is.finite(k) & k >= -1e-8
  result <- rep(NA_real_, length(k))
  result[valid] <- sqrt(pmax(k[valid], 0) / pi) - r[valid]
  result
}
for (s in 1:2) {
  d <- read.csv(inputs[s + 1L], colClasses = "character")
  pts <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), 32749)
  xy <- st_coordinates(pts) / 1000
  stopifnot(all(inside.owin(xy[, 1], xy[, 2], window)))
  aug <- as.Date(d$acq_date) >= as.Date("2026-08-01")
  masks <- list(rep(TRUE, nrow(d)), !aug, aug)
  for (p in 1:3) {
    keep <- masks[[p]]
    x <- ppp(xy[keep, 1], xy[keep, 2], window = window)
    expected <- window_prior$estimates[[labels[s]]][[periods[p]]]
    stopifnot(identical(d$source_file[keep], expected$source_file),
      identical(d$source_row[keep], expected$source_row))
    bd <- bdist.points(x)
    support <- data.frame(subset = labels[s], period = periods[p], r_km = 0:10,
      eligible_centres = vapply(0:10, function(a) sum(bd > a), integer(1)), records = npoints(x))
    support$fraction <- support$eligible_centres / support$records
    stopifnot(all(support$eligible_centres > 0), all(diff(support$eligible_centres) <= 0))
    support_frames[[length(support_frames) + 1L]] <- support
    # Reuse the verified full-period ordinary estimate without refitting.
    ko <- if (p == 1) ordinary_prior$results[[s]]$K else Kest(x, r = r, correction = "border")
    kd <- as.data.frame(ko)
    stopifnot(identical(kd$r, r), all(is.finite(kd$border)), all(kd$border >= -1e-8))
    ordinary_frames[[length(ordinary_frames) + 1L]] <- data.frame(subset = labels[s], period = periods[p],
      r_km = r, L_minus_r_km = to_L(kd$border), records = npoints(x))
    result <- list(source_file = d$source_file[keep], source_row = d$source_row[keep],
      ordinary_K = ko, border_support = support, reweighted = list())
    for (sigma in sigmas) {
      cat(labels[s], periods[p], "sigma", sigma, "km; n =", npoints(x), "\n")
      flush.console()
      if (p == 1) {
        old <- weighted_prior$estimates[[paste0("conf", c(30, 80)[s], "_sigma", sigma)]]
        stopifnot(identical(old$source_file, d$source_file), identical(old$source_row, d$source_row))
        lambda <- as.numeric(old$lambda_at_points)
        ki <- old$Kinhom
      } else {
        lambda <- as.numeric(density(x, sigma = sigma, at = "points", leaveoneout = TRUE,
          kernel = "gaussian", edge = TRUE, diggle = TRUE, eps = 0.5))
        stopifnot(all(is.finite(lambda)), all(lambda > 0))
        ki <- Kinhom(x, lambda = lambda, r = r, correction = "border", renormalise = TRUE, normpower = 2)
      }
      ik <- as.data.frame(ki)
      stopifnot(identical(ik$r, r), length(lambda) == npoints(x), all(is.finite(lambda) & lambda > 0))
      inverse <- 1 / lambda
      shares <- inverse / sum(inverse)
      # Check both all-record concentration and the border-eligible denominator at 10 km.
      eligible <- bd > 10
      interior_shares <- inverse[eligible] / sum(inverse[eligible])
      top <- which.max(shares)
      audit_frames[[length(audit_frames) + 1L]] <- data.frame(subset = labels[s], period = periods[p],
        sigma_km = sigma, records = npoints(x), minimum_intensity = min(lambda),
        maximum_reciprocal_share_percent = 100 * max(shares), reciprocal_effective_count = 1 / sum(shares^2),
        normalization_c = area_km2 / sum(inverse), eligible_centres_10km = sum(eligible),
        maximum_eligible_reciprocal_share_percent_10km = 100 * max(interior_shares),
        eligible_reciprocal_effective_count_10km = 1 / sum(interior_shares^2),
        top_weight_source_file = d$source_file[keep][top], top_weight_source_row = d$source_row[keep][top],
        raw_negative_roundoff_values = sum(is.finite(ik$border) & ik$border < 0 & ik$border >= -1e-8),
        invalid_K_values = sum(!is.finite(ik$border) | ik$border < -1e-8))
      weighted_frames[[length(weighted_frames) + 1L]] <- data.frame(subset = labels[s], period = periods[p],
        sigma_km = sigma, r_km = r, L_minus_r_km = to_L(ik$border))
      result$reweighted[[as.character(sigma)]] <- list(lambda = lambda, Kinhom = ki)
    }
    estimates[[paste(labels[s], periods[p], sep = "|")]] <- result
  }
}
ordinary <- do.call(rbind, ordinary_frames)
weighted <- do.call(rbind, weighted_frames)
audit <- do.call(rbind, audit_frames)
support <- do.call(rbind, support_frames)
ordinary$period <- factor(ordinary$period, levels = periods)
weighted$period <- factor(weighted$period, levels = periods)
weighted$bandwidth <- factor(paste0("Sigma = ", weighted$sigma_km, " km"), levels = paste0("Sigma = ", sigmas, " km"))
support$period <- factor(support$period, levels = periods)
colours <- setNames(c("#334155", "#D55E00", "#0072B2"), periods)
linetypes <- setNames(c("dashed", "solid", "solid"), periods)
base_theme <- theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"),
  legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9),
  plot.caption.position = "plot", plot.margin = margin(18, 18, 18, 18))
p1 <- ggplot(ordinary, aes(r_km, L_minus_r_km, colour = period, linetype = period)) +
  geom_hline(yintercept = 0, colour = "#64748B") + geom_line(linewidth = 0.8) +
  facet_wrap(~subset, nrow = 1) + scale_colour_manual(values = colours) + scale_linetype_manual(values = linetypes) +
  labs(title = "Ordinary spatial L summaries across time windows", subtitle = "Ketapang MODIS | Shared axes | Same border correction and 0-10 km grid",
    x = "Distance r (km)", y = "L(r) - r (km)", colour = NULL, linetype = NULL,
    caption = paste("Zero is a homogeneous Poisson theoretical reference, not a significance threshold. Spatial intensity variation remains unadjusted.",
      "Full period and August overlap. January-July has only 436 / 83 records (confidence >=30 / >=80), compared with 4,036 / 1,643 in August.",
      "Descriptive summaries only: no null simulations, p-values or independent-fire interpretation. Repeated detections remain in the data.", sep = "\n")) + base_theme
p2 <- ggplot(weighted, aes(r_km, L_minus_r_km, colour = period, linetype = period)) +
  geom_hline(yintercept = 0, colour = "#64748B") + geom_line(linewidth = 0.8, na.rm = TRUE) +
  facet_grid(subset ~ bandwidth) + scale_colour_manual(values = colours) + scale_linetype_manual(values = linetypes) +
  labs(title = "Intensity-reweighted spatial L summaries across time windows",
    subtitle = "Shared axes | Same 10 and 20 km illustrative bandwidths | Check weight stability before interpretation",
    x = "Distance r (km)", y = "L(r) - r (km)", colour = NULL, linetype = NULL,
    caption = paste("Intensity: period-specific Gaussian leave-one-out at points; Jones-Diggle correction; eps 0.5 km. Full-period estimates reused.",
      "Kinhom: border correction; renormalise TRUE; normpower 2. Zero is a Poisson reference, not a significance threshold.",
      "Sparse-period intensity weights may dominate. A curve near zero or below zero does not establish randomness or inhibition.",
      "Nested periods and confidence subsets; no simulations or final model. Raw K and all weight diagnostics are retained.", sep = "\n")) + base_theme
p3 <- ggplot(support, aes(r_km, eligible_centres, colour = period, linetype = period)) +
  geom_line(linewidth = 0.8) + facet_wrap(~subset, nrow = 1) +
  scale_colour_manual(values = colours) + scale_linetype_manual(values = linetypes) +
  labs(title = "Available interior centres for border correction", subtitle = "Absolute record counts; shared axes | Not a measure of independent sample size",
    x = "Border exclusion distance (km)", y = "Eligible observed centres", colour = NULL, linetype = NULL,
    caption = "An eligible centre is farther than r from the study-window boundary. Weighted effective counts are recorded separately in checks.txt.") + base_theme
stopifnot(identical(before, tools::md5sum(inputs)), identical(previous_hash, tools::md5sum(previous_files)))
dir.create(out)
ggsave(file.path(out, "ordinary_L_by_time_window.png"), p1, width = 12, height = 7, dpi = 160, bg = "white")
ggsave(file.path(out, "reweighted_L_by_time_window.png"), p2, width = 12, height = 9, dpi = 160, bg = "white")
ggsave(file.path(out, "border_support_by_time_window.png"), p3, width = 12, height = 6, dpi = 160, bg = "white")
saveRDS(list(estimates = estimates, ordinary = ordinary, weighted = weighted, audit = audit, support = support,
  source_md5 = before, prior_md5 = previous_hash, settings = list(projection = "EPSG:32749", units = "km",
    r = r, sigma_km = sigmas, density_edge = "Jones-Diggle", density_eps_km = 0.5,
    K_edge = "border", leaveoneout = TRUE, renormalise = TRUE, normpower = 2, final_test = FALSE)),
  file.path(out, "time_window_second_order.rds"))
capture.output({
  cat("EXPLORATORY TIME-WINDOW SPATIAL SUMMARIES, NOT INFERENCE.\n\n")
  print(audit, row.names = FALSE, digits = 7)
  cat("\nBORDER SUPPORT AT 10 KM\n"); print(support[support$r_km == 10, ], row.names = FALSE)
  cat("\nThe January-August estimates are reused unchanged after source identifier and checksum checks.\n")
  cat("New period-specific intensities: Gaussian leave-one-out, Jones-Diggle, eps 0.5 km, sigma 10/20 km.\n")
  cat("Both bandwidths and the shared 0-10 km distance range are provisional sensitivity settings.\n")
  cat("All spatial points retained; full polygon unchanged; no thinning, event merging, clipping of weights or intensity flooring.\n")
  cat("Reciprocal effective counts measure weight concentration, not independent events or effective statistical sample sizes.\n")
  cat("The largest reciprocal share is not a share of K pair contributions. Border-eligible shares are separately reported.\n")
  cat("Small early-period samples, different processing status and observation opportunity affect comparability.\n")
  cat("A lack of visible departure or a negative transformed curve is not evidence of randomness or inhibition.\n")
  cat("No fitted null model, simulations, confidence envelopes, significance tests, final bandwidths or fire-risk inference.\n")
  cat("Raw K estimates are preserved. Values in [-1e-8, 0) are bounded only for sqrt display; invalid values become display NA.\n")
  cat("Invalid K values across all weighted curves:", sum(audit$invalid_K_values), "\n")
  cat("\nSource MD5 unchanged:\n"); print(before)
  cat("\nPrior result MD5 unchanged:\n"); print(previous_hash)
  cat("\nReferences:\nhttps://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kinhom.html\n")
  cat("https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kest.html\n")
  print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(audit[, c("subset", "period", "sigma_km", "records", "maximum_reciprocal_share_percent",
  "reciprocal_effective_count", "eligible_centres_10km", "invalid_K_values")], row.names = FALSE)
cat("Saved exploratory summaries to", out, "\n")
