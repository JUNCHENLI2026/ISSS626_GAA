# Independent arithmetic, permutation invariants and grid-assignment checks.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], mustWork = TRUE)
x <- readRDS(file.path(root, "day_block_tests/day_block_results.rds"))
prior <- readRDS(file.path(root, "temporal_checks/temporal_diagnostics.rds"))
spatial <- readRDS(file.path(root, "spatial_recurrence/spatial_recurrence.rds"))
stopifnot(identical(unname(x$source_md5), unname(spatial$source_md5)))
for (key in names(x$cases)) {
  z <- x$cases[[key]]; lab <- z$summary$subset[1]
  full <- prior$subsets[[lab]]
  id <- match(paste(z$source_file, z$source_row), paste(full$records$source_file, full$records$source_row))
  stopifnot(!anyNA(id))
  distances <- as.matrix(dist(full$xy_km[id, , drop = FALSE]))
  ij <- which(lower.tri(distances) & distances <= 2, arr.ind = TRUE)
  ds <- distances[ij]
  count <- function(assignment) {
    times <- assignment[z$day_id] * 24 + z$hour_of_day
    gap <- abs(times[ij[,1]] - times[ij[,2]])
    vapply(seq_len(9), function(k) sum(ds <= z$summary$distance_km[k] & gap >= 1 & gap <= z$summary$max_gap_hours[k]), integer(1))
  }
  stopifnot(identical(count(as.numeric(z$dates)), as.integer(z$summary$observed_pairs)))
  for (b in 1:3) {
    stopifnot(identical(count(z$checked_assignments[[b]]), z$simulations[b, ]))
    for (g in split(seq_along(z$dates), z$grouping)) stopifnot(identical(sort(z$checked_assignments[[b]][g]), sort(as.numeric(z$dates[g]))))
  }
  p <- (1 + colSums(sweep(z$simulations, 2, z$summary$observed_pairs, ">="))) / 10000
  stopifnot(max(abs(p - z$summary$p_mc)) < 1e-12)
  cat("Verified direct all-pairs counts and saved permutations:", key, "\n")
  rm(distances); gc(verbose = FALSE)
}
p <- c(x$earlier_tests$p_mc, x$summary$p_mc)
ord <- order(p); corrected <- pmin(1, cummax((length(p):1) * p[ord])); restored <- numeric(length(p)); restored[ord] <- corrected
stopifnot(max(abs(restored - c(x$earlier_tests$p_holm_144, x$summary$p_holm_144))) < 1e-12)
for (key in names(spatial$cells)) {
  tab <- spatial$cells[[key]]; lab <- tab$subset[1]
  d <- prior$subsets[[lab]]$records; xy <- prior$subsets[[lab]]$xy_km
  size <- tab$cell_km[1]; shift <- size * tab$shift_fraction[1]
  ids <- paste(floor((xy[,1] - shift)/size), floor((xy[,2] - shift)/size), sep = "_")
  stopifnot(identical(ids, spatial$assignments[[key]]$cell))
  counts <- table(ids)
  days <- tapply(d$acq_date, ids, function(v) length(unique(v)))
  stopifnot(all(as.integer(counts[tab$cell]) == tab$records), all(as.integer(days[tab$cell]) == tab$recorded_dates))
}
stopifnot(identical(unname(x$source_md5),
  unname(tools::md5sum(file.path(root, basename(names(x$source_md5)))))))
cat("All 144 adjusted p-values and all eight spatial grid assignments verified. Source hashes unchanged.\n")
