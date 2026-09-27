# Exploratory robustness audit fixed before this extension's results are inspected.
# Rscript robustness_extension.R <Ketapang directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "robustness_extension")
if (dir.exists(out)) stop("Output exists; use a fresh rebuild directory.")
prior <- readRDS(file.path(root, "temporal_checks/temporal_diagnostics.rds"))
old <- readRDS(file.path(root, "permutation_tests/conditional_permutation_results.rds"))
day <- readRDS(file.path(root, "day_block_tests/day_block_results.rds"))
spatial <- readRDS(file.path(root, "spatial_recurrence/spatial_recurrence.rds"))
inputs <- file.path(root, c("ketapang_boundary.gpkg", "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
hashes <- tools::md5sum(inputs)
stopifnot(identical(unname(hashes), unname(day$source_md5)),
  identical(unname(hashes), unname(spatial$source_md5)))
B <- 9999L
lower <- c(0, 1, 3, 6)
schemes <- c("Daily within week", "Three-day blocks offset 0", "Three-day blocks offset 1")
results <- cases <- list()
case_id <- 0L
for (lab in c("Confidence >=30", "Confidence >=80")) for (period in c("January-August", "August")) {
  full <- prior$subsets[[lab]]
  keep <- if (period == "August") substr(full$records$acq_date, 1, 7) == "2026-08" else rep(TRUE, nrow(full$records))
  d <- full$records[keep, ]
  remap <- rep(NA_integer_, length(keep)); remap[keep] <- seq_len(nrow(d))
  p <- full$pairs_within_2km
  p <- p[keep[p$i] & keep[p$j] & p$distance_km <= 1, ]
  p$i <- remap[p$i]; p$j <- remap[p$j]
  dates <- seq(as.Date(if (period == "August") "2026-08-01" else "2026-01-01"), as.Date("2026-08-31"), by = "day")
  calendar <- as.numeric(dates); month <- format(dates, "%Y-%m")
  id <- match(as.Date(d$acq_date), dates)
  tod <- as.integer(substr(d$acq_time, 1, 2)) + as.integer(substr(d$acq_time, 3, 4)) / 60
  time <- calendar[id] * 24 + tod
  gaps <- abs(time[p$i] - time[p$j])
  count <- function(g) vapply(lower, function(lo) as.integer(sum(g >= lo & g <= 24)), integer(1))
  observed <- count(gaps)
  for (scheme in schemes) {
    case_id <- case_id + 1L
    groups <- list(); units <- rep(-1L, length(dates)); uid <- 0L
    if (scheme == schemes[1]) {
      week <- (as.integer(format(dates, "%d")) - 1L) %/% 7L
      groups <- lapply(split(seq_along(dates), interaction(month, week, drop = TRUE)), function(g) matrix(g, ncol = 1))
    } else {
      offset <- if (scheme == schemes[2]) 0L else 1L
      for (g in split(seq_along(dates), month)) {
        usable <- g[seq.int(1L + offset, length(g))]
        n <- length(usable) %/% 3L
        groups[[length(groups) + 1L]] <- matrix(usable[seq_len(n * 3L)], ncol = 3, byrow = TRUE)
      }
    }
    for (g in groups) for (j in seq_len(nrow(g))) {
      uid <- uid + 1L; units[g[j, ]] <- uid
    }
    fixed <- units[id[p$i]] == units[id[p$j]]
    fixed_counts <- count(gaps[fixed])
    varying <- p[!fixed, ]
    simulations <- matrix(NA_integer_, B, length(lower))
    checked <- list()
    seed <- 2026092300L + case_id
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    key <- paste(lab, period, scheme, sep = " | ")
    cat("START", key, "\n"); flush.console()
    for (b in seq_len(B)) {
      assigned <- calendar
      for (g in groups) assigned[as.vector(t(g))] <- calendar[as.vector(t(g[sample.int(nrow(g)), , drop = FALSE]))]
      newtime <- assigned[id] * 24 + tod
      simulations[b, ] <- fixed_counts + count(abs(newtime[varying$i] - newtime[varying$j]))
      if (b <= 3L) {
        stopifnot(identical(simulations[b, ], count(abs(newtime[p$i] - newtime[p$j]))),
          identical(sort(assigned), calendar), all(format(as.Date(assigned, origin = "1970-01-01"), "%Y-%m") == month),
          all(assigned[units == -1L] == calendar[units == -1L]))
        for (g in groups) for (j in seq_len(nrow(g))) stopifnot(all(diff(assigned[g[j, ]]) == diff(calendar[g[j, ]])))
        checked[[b]] <- assigned
      }
    }
    exceed <- colSums(sweep(simulations, 2, observed, ">="))
    means <- colMeans(simulations)
    ci <- t(vapply(exceed, function(k) binom.test(k, B)$conf.int, numeric(2)))
    tab <- data.frame(subset = lab, period = period, scheme = scheme, lower_hours = lower,
      observed = observed, guaranteed_fixed = fixed_counts, variable_observed = observed - fixed_counts,
      null_mean = means, variable_null_mean = means - fixed_counts,
      total_ratio = observed / means, variable_ratio = (observed - fixed_counts) / (means - fixed_counts),
      fixed_percent = 100 * fixed_counts / observed, p_mc = (exceed + 1) / (B + 1),
      exceedances = exceed, tail_ci_low = ci[, 1], tail_ci_high = ci[, 2],
      null_q025 = apply(simulations, 2, quantile, .025), null_q975 = apply(simulations, 2, quantile, .975), seed = seed)
    results[[key]] <- tab
    cases[[key]] <- list(summary = tab, simulations = simulations, checked_assignments = checked,
      dates = dates, units = units, groups = groups, day_id = id, hour_of_day = tod,
      source_file = d$source_file, source_row = d$source_row)
    cat("DONE", key, "raw p", paste(tab$p_mc, collapse = ", "), "\n"); flush.console()
  }
}
summary <- do.call(rbind, results)
stopifnot(nrow(summary) == 48L)
# Four repeated primary settings remain in the family: deliberate conservative duplication.
joint <- p.adjust(c(old$summary$p_mc, day$summary$p_mc, summary$p_mc), "holm")
summary$p_holm_192 <- tail(joint, 48)

# Evaluate the same 621 occupied reference-cell centres across all eight grids.
# This measures consistency on a fixed reference support, not the fraction of regency area.
base <- spatial$cells[["Confidence >=30 | 5 | 0"]]
centres <- data.frame(cell = base$cell, x_km = (base$gx + .5) * 5, y_km = (base$gy + .5) * 5)
stability <- list()
for (key in names(spatial$cells)) {
  tab <- spatial$cells[[key]]
  size <- tab$cell_km[1]; shift <- size * tab$shift_fraction[1]
  ids <- paste(floor((centres$x_km - shift) / size), floor((centres$y_km - shift) / size), sep = "_")
  hit <- match(ids, tab$cell)
  # Maximum tie ranks avoid arbitrarily selecting a fraction of an equal-date group.
  percent <- 100 * rank(-tab$recorded_dates, ties.method = "max") / nrow(tab)
  top <- !is.na(hit) & percent[hit] <= 10
  stability[[key]] <- data.frame(reference_cell = centres$cell, setting = key,
    subset = tab$subset[1], cell_km = size, shift_fraction = tab$shift_fraction[1],
    containing_cell = ids, dates = ifelse(is.na(hit), 0, tab$recorded_dates[hit]),
    upper_rank_percent = percent[hit], top_decile = top)
}
stability <- do.call(rbind, stability)
support <- aggregate(top_decile ~ reference_cell, stability, sum)
names(support)[2] <- "settings_in_top_decile"
support <- merge(centres, support, by.x = "cell", by.y = "reference_cell", sort = FALSE)
anchors <- stability[stability$reference_cell %in% c("80_1969", "81_1957"), ]
stopifnot(nrow(anchors) == 16L, nrow(support) == 621L,
  identical(hashes, tools::md5sum(inputs)))
dir.create(out)
saveRDS(list(summary = summary, cases = cases, previous_p_holm_192 = head(joint, 144),
  stability = stability, support = support, anchors = anchors,
  settings = list(B = B, lower_hours = lower, distance_km = 1, max_gap_hours = 24,
    family_size = 192L, grid_top_percent = 10, reference_support = "621 occupied primary-grid centres"),
  source_md5 = hashes, script_md5 = tools::md5sum(file.path(root, "robustness_extension.R"))),
  file.path(out, "robustness_results.rds"))
write.csv(summary, file.path(out, "permutation_summary.csv"), row.names = FALSE)
write.csv(anchors, file.path(out, "anchor_stability.csv"), row.names = FALSE)
write.csv(support, file.path(out, "reference_support.csv"), row.names = FALSE)
capture.output({print(summary); print(anchors); print(table(support$settings_in_top_decile)); print(sessionInfo())},
  file = file.path(out, "checks.txt"))
cat("COMPLETE. New comparisons", nrow(summary), "significant after Holm192", sum(summary$p_holm_192 < .05), "\n")
