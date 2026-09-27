# Independent all-pairs checks for the robustness extension. Does not alter results.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], mustWork = TRUE)
x <- readRDS(file.path(root, "robustness_extension/robustness_results.rds"))
prior <- readRDS(file.path(root, "temporal_checks/temporal_diagnostics.rds"))
old <- readRDS(file.path(root, "permutation_tests/conditional_permutation_results.rds"))
day <- readRDS(file.path(root, "day_block_tests/day_block_results.rds"))
for (key in names(x$cases)) {
  z <- x$cases[[key]]
  full <- prior$subsets[[z$summary$subset[1]]]
  id <- match(paste(z$source_file, z$source_row), paste(full$records$source_file, full$records$source_row))
  stopifnot(!anyNA(id))
  distmat <- as.matrix(dist(full$xy_km[id, , drop = FALSE]))
  pairs <- which(lower.tri(distmat) & distmat <= 1, arr.ind = TRUE)
  evaluate <- function(assigned) {
    time <- assigned[z$day_id] * 24 + z$hour_of_day
    gap <- abs(time[pairs[, 1]] - time[pairs[, 2]])
    vapply(c(0, 1, 3, 6), function(lo) as.integer(sum(gap >= lo & gap <= 24)), integer(1))
  }
  stopifnot(identical(evaluate(as.numeric(z$dates)), z$summary$observed))
  for (b in 1:3) {
    assigned <- z$checked_assignments[[b]]
    stopifnot(identical(evaluate(assigned), z$simulations[b, ]),
      identical(sort(assigned), as.numeric(z$dates)),
      identical(format(as.Date(assigned, origin = "1970-01-01"), "%Y-%m"), format(z$dates, "%Y-%m")))
    for (g in z$groups) for (j in seq_len(nrow(g))) {
      stopifnot(all(diff(assigned[g[j, ]]) == diff(as.numeric(z$dates[g[j, ]]))))
    }
  }
  p <- (1 + colSums(sweep(z$simulations, 2, z$summary$observed, ">="))) / (nrow(z$simulations) + 1)
  stopifnot(isTRUE(all.equal(unname(p), z$summary$p_mc)))
  cat("Verified", key, "\n")
}
p <- c(old$summary$p_mc, day$summary$p_mc, x$summary$p_mc)
o <- order(p); adjusted <- numeric(length(p))
adjusted[o] <- pmin(1, cummax((length(p):1) * p[o]))
stopifnot(length(p) == 192L,
  isTRUE(all.equal(tail(adjusted, 48), x$summary$p_holm_192)),
  isTRUE(all.equal(head(adjusted, 144), x$previous_p_holm_192)))
spatial <- readRDS(file.path(root, "spatial_recurrence/spatial_recurrence.rds"))
base <- spatial$cells[["Confidence >=30 | 5 | 0"]]
for (key in unique(x$stability$setting)) {
  tab <- spatial$cells[[key]]
  result <- x$stability[x$stability$setting == key, ]
  ref <- match(result$reference_cell, base$cell)
  width <- tab$cell_km[1]; shift <- width * tab$shift_fraction[1]
  expected <- paste(floor(((base$gx[ref] + .5) * 5 - shift) / width),
    floor(((base$gy[ref] + .5) * 5 - shift) / width), sep = "_")
  hit <- match(expected, tab$cell)
  # Direct counting verifies the maximum-tie rank independently of rank().
  fractions <- vapply(tab$recorded_dates, function(v) 100 * mean(tab$recorded_dates >= v), numeric(1))
  expected_top <- !is.na(hit) & fractions[hit] <= 10
  stopifnot(identical(expected, result$containing_cell), identical(expected_top, result$top_decile))
}
stopifnot(identical(unname(x$source_md5),
  unname(tools::md5sum(file.path(root, basename(names(x$source_md5)))))),
  identical(unname(x$script_md5), unname(tools::md5sum(file.path(root, "robustness_extension.R")))))
cat("PASS: 192-test Holm family, 12 cases, sampled block permutations and all eight grid projections.\n")
