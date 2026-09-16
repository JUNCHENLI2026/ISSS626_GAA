# Exploratory second-order summaries, not a significance test.
# Usage: Rscript second_order_diagnostics.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "second_order")
if (dir.exists(out)) stop("Output directory already exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
stopifnot(abs(area.owin(window) - sum(as.numeric(st_area(boundary))) / 1e6) < 1e-4)
labels <- c("Confidence >=30 (n = 4,472)", "Confidence >=80 (n = 1,726)")
expected <- c(4472L, 1726L)
results <- vector("list", 2)
# A bounded initial diagnostic range, not a final inferential distance range.
rmax_km <- 10
for (i in seq_along(results)) {
  records <- read.csv(inputs[i + 1L], colClasses = "character")
  points <- st_transform(st_as_sf(records, coords = c("longitude", "latitude"), crs = 4326), 32749)
  xy <- st_coordinates(points) / 1000
  stopifnot(nrow(xy) == expected[i], all(inside.owin(xy[, 1], xy[, 2], window)))
  pattern <- ppp(xy[, 1], xy[, 2], window = window, checkdup = TRUE)
  stopifnot(npoints(pattern) == expected[i])
  cat("Computing border-corrected K summary:", labels[i], "\n")
  flush.console()
  # Border correction is explicit; it is computationally tractable for the
  # detailed polygon and large point pattern. Its loss of centres is audited.
  k <- Kest(pattern, rmax = rmax_km, correction = "border")
  kd <- as.data.frame(k)
  stopifnot(all(c("r", "border", "theo") %in% names(kd)),
    all(is.finite(kd$border)), all(kd$border >= 0),
    max(abs(kd$theo - pi * kd$r^2)) < 1e-8)
  distances <- bdist.points(pattern)
  support_r <- 0:10
  support_n <- vapply(support_r, function(r) sum(distances > r), integer(1))
  stopifnot(all(diff(support_n) <= 0), all(support_n > 0))
  results[[i]] <- list(label = labels[i], records = nrow(records),
    distinct_dates = length(unique(records$acq_date)),
    duplicate_coordinates = sum(duplicated(as.data.frame(xy))),
    K = k,
    curve = data.frame(r_km = kd$r, L_minus_r_km = sqrt(kd$border / pi) - kd$r,
      subset = labels[i]),
    support = data.frame(r_km = support_r, eligible_centres = support_n,
      fraction = support_n / nrow(records), subset = labels[i]))
}
curves <- do.call(rbind, lapply(results, `[[`, "curve"))
support <- do.call(rbind, lapply(results, `[[`, "support"))
palette <- setNames(c("#0072B2", "#D55E00"), labels)
base_theme <- theme_minimal(base_size = 12) + theme(
  plot.title = element_text(face = "bold"), legend.position = "bottom",
  legend.title = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
  plot.caption.position = "plot", plot.margin = margin(15, 15, 15, 15))
p1 <- ggplot(curves, aes(r_km, L_minus_r_km, colour = subset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.85) + scale_colour_manual(values = palette) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(title = "Exploratory L(r) - r summaries",
    subtitle = "Ketapang MODIS | January-August 2026 | Border correction",
    x = "Distance r (km)", y = "L(r) - r (km)",
    caption = paste("Dashed zero line: homogeneous Poisson theoretical reference, NOT a significance threshold.",
      "No simulations or p-values. Spatial intensity variation and repeated detections are not adjusted for.",
      "The 0-10 km range is provisional. Nested confidence subsets are not independent samples.",
      "Sources: NASA FIRMS; supplied Indonesia Geospatial boundary.", sep = "\n")) + base_theme
p2 <- ggplot(support, aes(r_km, fraction * 100, colour = subset)) +
  geom_line(linewidth = 0.85) + geom_point(size = 2) +
  scale_colour_manual(values = palette) + scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(title = "Border-correction support by distance",
    subtitle = "Percentage of observed points more than r km from the window boundary",
    x = "Distance r (km)", y = "Eligible centres (%)",
    caption = paste("Border correction uses interior points as centres; neighbouring points are not filtered in the same way.",
      "This is a support diagnostic, not an uncertainty interval or proof that a distance range is valid.",
      "The supplied polygon retains holes and disconnected components; no points were removed from input files.", sep = "\n")) + base_theme
stopifnot(identical(before, tools::md5sum(inputs)))
dir.create(out, recursive = TRUE)
ggsave(file.path(out, "observed_L_summary.png"), p1, width = 10, height = 7, dpi = 160, bg = "white")
ggsave(file.path(out, "border_support.png"), p2, width = 10, height = 7, dpi = 160, bg = "white")
saveRDS(list(results = results, projection = "EPSG:32749", units = "km",
  provisional_rmax_km = rmax_km, source_md5 = before), file.path(out, "observed_summaries.rds"))
capture.output({
  cat("EXPLORATORY DIAGNOSTICS: not a significance test or final model.\n")
  cat("Projection: EPSG:32749. Coordinate and distance units: km.\n")
  cat("Window area (km2):", area.owin(window), "\n")
  cat("Estimator: homogeneous Kest, border correction; L = sqrt(K/pi).\n")
  cat("Initial diagnostic range: 0-10 km, not a final justified test range.\n")
  cat("All input points retained inside the full polygon; area conversion verified.\n")
  for (z in results) {
    cat("\n", z$label, "\nDistinct dates:", z$distinct_dates,
      "\nExact duplicate coordinates:", z$duplicate_coordinates, "\n")
    print(z$support[, c("r_km", "eligible_centres", "fraction")], row.names = FALSE)
  }
  cat("\nNo temporal interaction model, intensity adjustment, simulations or p-values.\n")
  cat("A homogeneous reference is not established as appropriate for these data.\n")
  cat("Boundary support does not address spatially varying observation opportunity.\n")
  cat("Input files unchanged; MD5:\n"); print(before)
  cat("\nReferences:\nhttps://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kest.html\n")
  cat("https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Lest.html\n")
  print(sessionInfo())
}, file = file.path(out, "checks.txt"))
stopifnot(identical(before, tools::md5sum(inputs)))
cat("Completed. Input checksums unchanged. Outputs:", out, "\n")
