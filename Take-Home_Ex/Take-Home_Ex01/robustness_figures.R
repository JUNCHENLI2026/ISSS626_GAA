# Figures use saved, independently checked results without rerunning simulations.
reference_assignment_sensitivity <- function(recurrence) {
  base <- recurrence$cells[["Confidence >=30 | 5 | 0"]]
  offsets <- expand.grid(dx = c(-.001, .001), dy = c(-.001, .001))
  counts <- matrix(0L, nrow(base), nrow(offsets))
  for (q in seq_len(nrow(offsets))) for (tab in recurrence$cells) {
    width <- tab$cell_km[1]; shift <- width * tab$shift_fraction[1]
    ids <- paste(floor(((base$gx + .5) * 5 + offsets$dx[q] - shift) / width),
      floor(((base$gy + .5) * 5 + offsets$dy[q] - shift) / width), sep = "_")
    hit <- match(ids, tab$cell)
    cutoff <- rank(-tab$recorded_dates, ties.method = "max") / nrow(tab) <= .1
    counts[, q] <- counts[, q] + as.integer(!is.na(hit) & cutoff[hit])
  }
  # Independent rectangle-inequality and direct-count check of all 32 combinations.
  direct <- matrix(0L, nrow(base), nrow(offsets))
  for (q in seq_len(nrow(offsets))) for (tab in recurrence$cells) {
    width <- tab$cell_km[1]; shift <- width * tab$shift_fraction[1]
    left <- tab$gx * width + shift; bottom <- tab$gy * width + shift
    for (i in seq_len(nrow(base))) {
      x <- (base$gx[i] + .5) * 5 + offsets$dx[q]
      y <- (base$gy[i] + .5) * 5 + offsets$dy[q]
      j <- which(left <= x & x < left + width & bottom <= y & y < bottom + width)
      stopifnot(length(j) <= 1L)
      if (length(j)) direct[i, q] <- direct[i, q] + as.integer(mean(tab$recorded_dates >= tab$recorded_dates[j]) <= .1)
    }
  }
  stopifnot(identical(counts, direct))
  data.frame(cell = base$cell, min_settings = apply(counts, 1, min),
    max_settings = apply(counts, 1, max), all_32 = rowSums(counts) == 32L,
    district_context = base$district)
}

plot_pair_denominators <- function(audit) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- subset(audit$summary, subset == "Confidence >=30" & period == "January-August" & lower_hours == 1)
  labels <- c("Daily patterns", "3-day blocks: start day 1", "3-day blocks: start day 2")
  d$model <- factor(labels, levels = rev(labels))
  long <- rbind(transform(d, measure = "All qualifying pairs", value = total_ratio),
    transform(d, measure = "Reassignable component", value = variable_ratio))
  ggplot2::ggplot(long, ggplot2::aes(value, model)) +
    ggplot2::geom_vline(xintercept = 1, colour = "#94a3b8") +
    ggplot2::geom_point(size = 3, colour = "#0072b2") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), nudge_y = .19, size = 4) +
    ggplot2::facet_wrap(~measure, nrow = 1) +
    ggplot2::scale_x_continuous(limits = c(.8, 5.2)) +
    ggplot2::labs(x = "Observed / permutation mean", y = NULL,
      title = "Fixed pairs change the enrichment denominator",
      subtitle = "Main subset, January-August, 1 km and 1-24 hours",
      caption = "Each model holds a different set of pairs fixed. Ratios are model-specific, not comparable fire-risk effects.") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 11))
}

plot_grid_stability <- function(audit, recurrence) {
  stopifnot(requireNamespace("sf", quietly = TRUE), requireNamespace("ggplot2", quietly = TRUE))
  cells <- recurrence$main_grid
  cells$agreements <- audit$support$settings_in_top_decile[match(cells$cell, audit$support$cell)]
  stopifnot(!anyNA(cells$agreements), sum(cells$agreements == 8) == 11L)
  # Plot the reference centres, not polygons, because membership is evaluated at a point.
  points <- sf::st_as_sf(transform(audit$support, x = x_km * 1000, y = y_km * 1000), coords = c("x", "y"), crs = 32749)
  points$agreement <- factor(points$settings_in_top_decile, levels = 0:8)
  stable <- points[points$settings_in_top_decile == 8, ]
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = recurrence$districts, fill = "#f1f5f9", colour = "#94a3b8", linewidth = .2) +
    ggplot2::geom_sf(data = points, ggplot2::aes(colour = agreement), size = 1.5, alpha = .85) +
    ggplot2::geom_sf(data = stable, shape = 21, size = 3, fill = NA, colour = "#111827", stroke = .7) +
    ggplot2::scale_colour_viridis_d(option = "C", drop = FALSE, name = "Settings in\ntop decile") +
    ggplot2::labs(title = "Recurrence ranking depends on spatial support",
      subtitle = "621 primary occupied-cell centres, assessed against eight grid/confidence settings",
      caption = "Outlines: 8/8 settings under the stated half-open assignment convention only.\nNone qualifies in all 32 direction/setting combinations. These are not significant hotspots or area fractions.",
      x = "Longitude", y = "Latitude") +
    ggplot2::coord_sf(datum = sf::st_crs(4326)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 10))
}
