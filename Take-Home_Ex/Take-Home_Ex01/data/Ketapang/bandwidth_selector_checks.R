# Exploratory bandwidth criteria and downstream weight stability, not inference.
# Usage: Rscript bandwidth_selector_checks.R <Ketapang data directory>
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "bandwidth_checks")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
inputs <- file.path(root, c("ketapang_boundary.gpkg",
  "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
before <- tools::md5sum(inputs)
prior_file <- file.path(root, "inhomogeneous_checks", "intensity_reweighted_summaries.rds")
prior_hash <- tools::md5sum(prior_file)
prior <- readRDS(prior_file)
stopifnot(identical(unname(before), unname(prior$source_md5)))
boundary <- st_transform(st_read(inputs[1], layer = "study_window", quiet = TRUE), 32749)
stopifnot(all(st_is_valid(boundary)), !any(st_is_empty(boundary)))
window <- rescale(as.owin(boundary), 1000, "km")
area_km2 <- area.owin(window)
stopifnot(abs(area_km2 - prior$window_area_km2) < 1e-6)
saved <- list()
objective_frames <- selection_frames <- audit_frames <- list()
labels <- c("Confidence >=30", "Confidence >=80")

search_selector <- function(x, method, label) {
  # Log-spaced computational search bounds, not a scientific admissibility rule.
  limits <- c(1, 40)
  history <- list()
  evaluate <- function(grid) {
    cat(label, method, "evaluating", length(grid), "bandwidths from",
        min(grid), "to", max(grid), "km\n")
    flush.console()
    if (method == "Likelihood CV") {
      # Include the integrated intensity term explicitly, without shortcut.
      bw.ppl(x, sigma = grid, shortcut = FALSE, warn = FALSE,
        kernel = "gaussian", edge = TRUE, diggle = TRUE, eps = 0.5)
    } else {
      # Installed bw.CvL ignores ...: inclusive point intensity, no edge correction.
      # Do not claim that its criterion uses the PPL/weight-audit settings.
      bw.CvL(x, sigma = grid, warn = FALSE)
    }
  }
  for (attempt in 1:3) {
    grid <- exp(seq(log(limits[1]), log(limits[2]), length.out = 20))
    fit <- evaluate(grid)
    history[[length(history) + 1L]] <- fit
    j <- attr(fit, "iopt")
    if (j > 1L && j < length(grid)) break
    # Expand once per endpoint encounter, up to three coarse searches total.
    if (j == 1L) limits[1] <- limits[1] / 4 else limits[2] <- limits[2] * 4
  }
  if (j > 1L && j < length(grid)) {
    fine <- exp(seq(log(grid[j - 1L]), log(grid[j + 1L]), length.out = 17))
    history[[length(history) + 1L]] <- evaluate(fine)
  }
  tab <- do.call(rbind, lapply(seq_along(history), function(k) {
    f <- history[[k]]
    data.frame(sigma_km = attr(f, "h"), objective = attr(f, "cv"), pass = k)
  }))
  tab <- tab[!duplicated(signif(tab$sigma_km, 12)), ]
  tab <- tab[order(tab$sigma_km), ]
  stopifnot(any(is.finite(tab$objective)), !anyNA(tab$objective))
  idx <- if (method == "Likelihood CV") which.max(tab$objective) else which.min(tab$objective)
  sigma <- tab$sigma_km[idx]
  best <- tab$objective[idx]
  # A non-negative distance from the best evaluated objective, within this panel.
  tab$objective_gap <- if (method == "Likelihood CV") best - tab$objective else tab$objective - best
  tab$subset <- label
  tab$method <- method
  summary <- data.frame(subset = label, method = method, sigma_km = sigma,
    objective = best, search_min_km = min(tab$sigma_km), search_max_km = max(tab$sigma_km),
    evaluated_candidates = nrow(tab), endpoint = idx %in% c(1L, nrow(tab)),
    lower_neighbour_km = if (idx > 1L) tab$sigma_km[idx - 1L] else NA_real_,
    upper_neighbour_km = if (idx < nrow(tab)) tab$sigma_km[idx + 1L] else NA_real_)
  list(summary = summary, curve = tab, fits = history)
}

for (i in 1:2) {
  records <- read.csv(inputs[i + 1L], colClasses = "character")
  pts <- st_transform(st_as_sf(records, coords = c("longitude", "latitude"), crs = 4326), 32749)
  xy <- st_coordinates(pts) / 1000
  stopifnot(nrow(xy) == c(4472L, 1726L)[i], all(inside.owin(xy[, 1], xy[, 2], window)))
  x <- ppp(xy[, 1], xy[, 2], window = window)
  for (method in c("Likelihood CV", "Cronie-van Lieshout")) {
    result <- search_selector(x, method, labels[i])
    selection_frames[[length(selection_frames) + 1L]] <- result$summary
    objective_frames[[length(objective_frames) + 1L]] <- result$curve
    # Audit selected value and +/-20% neighbours; these are sensitivities, not CIs.
    for (factor in c(0.8, 1, 1.2)) {
      sigma <- result$summary$sigma_km * factor
      lam <- as.numeric(density(x, sigma = sigma, at = "points", leaveoneout = TRUE,
        kernel = "gaussian", edge = TRUE, diggle = TRUE, eps = 0.5))
      positive <- all(is.finite(lam) & lam > 0)
      shares <- if (positive) (1 / lam) / sum(1 / lam) else rep(NA_real_, npoints(x))
      audit_frames[[length(audit_frames) + 1L]] <- data.frame(
        subset = labels[i], method = method, bandwidth_factor = factor, sigma_km = sigma,
        all_finite_positive = positive, minimum_intensity = min(lam),
        maximum_reciprocal_share_percent = if (positive) 100 * max(shares) else NA_real_,
        reciprocal_effective_count = if (positive) 1 / sum(shares^2) else NA_real_,
        normalization_c = if (positive) area_km2 / sum(1 / lam) else NA_real_)
      result[[paste0("lambda_factor_", factor)]] <- lam
    }
    saved[[paste0("conf", c(30, 80)[i], "_", method)]] <- result
  }
}
selection <- do.call(rbind, selection_frames)
objectives <- do.call(rbind, objective_frames)
audit <- do.call(rbind, audit_frames)
stopifnot(identical(before, tools::md5sum(inputs)),
          identical(prior_hash, tools::md5sum(prior_file)),
          all(objectives$objective_gap >= 0))
dir.create(out)
saveRDS(list(searches = saved, selection = selection, objectives = objectives,
  audit = audit, source_md5 = before, prior_md5 = prior_hash,
  settings = list(projection = "EPSG:32749", units = "km", eps_km = 0.5,
    ppl_shortcut = FALSE, ppl_edge = "Jones-Diggle", ppl_leaveoneout = TRUE,
    CvL_edge = FALSE, CvL_leaveoneout = FALSE,
    audit_edge = "Jones-Diggle", audit_leaveoneout = TRUE,
    factors = c(0.8, 1, 1.2))), file.path(out, "bandwidth_diagnostics.rds"))
p <- ggplot(subset(objectives, is.finite(objective_gap)), aes(sigma_km, objective_gap)) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_point(colour = "#0072B2", size = 1.3) +
  geom_vline(data = selection, aes(xintercept = sigma_km), colour = "#D55E00", linetype = "dashed") +
  facet_wrap(~subset + method, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 40, 80, 160)) +
  scale_y_continuous(transform = "log1p", breaks = c(0, 10, 100, 1e3, 1e4, 1e6, 1e8),
    labels = c("0", "10", "100", "1,000", "10,000", "1e6", "1e8")) +
  labs(title = "Exploratory bandwidth search", subtitle = "Lower is better within each panel; dashed line = best evaluated bandwidth",
    x = "Gaussian bandwidth sigma (km; log scale)", y = "Objective gap from best evaluated value (log1p scale)",
    caption = paste("Objective magnitudes are not comparable across panels. Grid minima are candidates, not final choices.",
      "Likelihood CV: Jones-Diggle correction, leave-one-out, full objective. CvL: no edge correction, includes focal point.",
      "Full-period MODIS detections; nested confidence subsets. No time model, simulations or significance test.", sep = "\n")) +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9), plot.margin = margin(14, 14, 14, 14))
ggsave(file.path(out, "bandwidth_objective_curves.png"), p, width = 12, height = 8, dpi = 160, bg = "white")
p2 <- ggplot(audit, aes(sigma_km, maximum_reciprocal_share_percent, colour = method)) +
  geom_line(linewidth = 0.7) + geom_point(aes(shape = bandwidth_factor == 1), size = 3) +
  facet_wrap(~subset, nrow = 1) + scale_colour_manual(values = c("#009E73", "#D55E00")) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("80% / 120% sensitivity", "Selected candidate")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Weight concentration after applying a common intensity estimator",
    subtitle = "Candidate bandwidths and +/-20% sensitivity; all original detections retained",
    x = "Gaussian bandwidth sigma (km)", y = "Largest share of reciprocal-intensity sum (%)", colour = NULL, shape = NULL,
    caption = paste("Every point uses Gaussian leave-one-out intensity, Jones-Diggle edge correction and a 0.5 km computational grid.",
      "For CvL this is a downstream transfer check, not the estimator used inside its criterion.",
      "Reciprocal shares are not K-function pair-contribution shares or independent-event percentages. No weight caps or deletions.", sep = "\n")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom",
    plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(14, 14, 14, 14)) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
ggsave(file.path(out, "candidate_weight_sensitivity.png"), p2, width = 12, height = 7, dpi = 160, bg = "white")
capture.output({
  cat("EXPLORATORY BANDWIDTH CRITERIA: no final bandwidth, null model or test.\n")
  cat("EPSG:32749; kilometres; window area:", area_km2, "km2. All observations retained.\n")
  cat("Initial search: 20 log-spaced values on 1-40 km. Expand an endpoint by factor 4, at most three coarse searches.\n")
  cat("Interior optimum: 17-point log-grid refinement between its adjacent coarse candidates.\n")
  cat("Search bounds and +/-20% audit are exploratory computational choices, not confidence intervals.\n")
  cat("Likelihood CV: sum(log(leave-one-out intensity)) minus integrated full-data intensity.\n")
  cat("Explicit Gaussian, Jones-Diggle edge correction, eps 0.5 km, shortcut FALSE.\n")
  cat("CvL: squared difference between sum(1/intensity) and window area.\n")
  cat("Installed CvL uses inclusive intensity at points with edge FALSE; extra arguments are ignored.\n")
  cat("Common downstream audit: Gaussian leave-one-out, Jones-Diggle, eps 0.5 km.\n")
  cat("Transferring CvL bandwidth to this different estimator is a sensitivity diagnostic, not an equivalent optimisation.\n\n")
  print(selection, row.names = FALSE, digits = 7)
  cat("\nCOMMON-ESTIMATOR WEIGHT AUDIT\n")
  print(audit, row.names = FALSE, digits = 7)
  cat("\nEffective count = 1/sum(normalised reciprocal weights squared), NOT independent events.\n")
  cat("No floors, caps, deletions, temporal thinning, K-function refits, simulation envelopes or p-values.\n")
  cat("Criteria can respond to repeated detections and small-scale structure; they do not identify first-order causes.\n")
  cat("An interior grid optimum does not establish scientific suitability, numerical-grid convergence or a unique continuous optimum.\n")
  cat("Intensity estimation, bandwidth choice and observation effects still require consideration in any fitted-null test.\n")
  cat("\nInput MD5 checksums unchanged:\n"); print(before)
  cat("\nPrevious result unchanged:\n"); print(prior_hash)
  cat("\nInstalled selector implementations (for reproducibility):\n"); print(bw.ppl); print(bw.CvL)
  cat("\nRuntime:\n"); print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(selection, row.names = FALSE)
print(audit, row.names = FALSE)
cat("Saved diagnostics to", out, "\n")
