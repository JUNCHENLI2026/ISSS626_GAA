# Whole-calendar-day permutations preserve every observed pattern within a UTC day.
# Run: Rscript day_block_permutation.R <Ketapang data directory>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "day_block_tests")
if (dir.exists(out)) stop("Output directory exists; preserve it before a deliberate rerun.")
prior <- readRDS(file.path(root, "temporal_checks/temporal_diagnostics.rds"))
old <- readRDS(file.path(root, "permutation_tests/conditional_permutation_results.rds"))
inputs <- file.path(root, c("ketapang_boundary.gpkg", "ketapang_modis_2026_conf30.csv", "ketapang_modis_2026_conf80.csv"))
hashes <- tools::md5sum(inputs)
stopifnot(identical(unname(hashes), unname(prior$source_md5)), identical(unname(hashes), unname(old$source_md5)))
B <- 9999L
thresholds <- expand.grid(max_gap_hours = c(6, 24, 72), distance_km = c(0.5, 1, 2))
schemes <- c("Whole days within month", "Whole days within month blocks")
periods <- c("January-August", "August")
labels <- c("Confidence >=30", "Confidence >=80")
all_results <- daily_checks <- cases <- list()
counter <- function(gap, spatial_bin) {
  valid <- gap >= 1 & gap <= 72
  temporal_bin <- 1L + (gap[valid] > 6) + (gap[valid] > 24)
  z <- matrix(tabulate((spatial_bin[valid] - 1L) * 3L + temporal_bin, nbins = 9), nrow = 3)
  z[2, ] <- z[2, ] + z[1, ]; z[3, ] <- z[3, ] + z[2, ]
  z[, 2] <- z[, 2] + z[, 1]; z[, 3] <- z[, 3] + z[, 2]
  as.integer(z)
}
case_index <- 0L
for (lab in labels) for (period in periods) {
  full <- prior$subsets[[lab]]
  d <- full$records
  keep <- if (period == periods[1]) rep(TRUE, nrow(d)) else substr(d$acq_date, 1, 7) == "2026-08"
  index <- which(keep); d <- d[keep, ]
  map <- rep(NA_integer_, nrow(full$records)); map[index] <- seq_along(index)
  p <- full$pairs_within_2km
  p <- p[keep[p$i] & keep[p$j], c("i", "j", "distance_km")]
  p$i <- map[p$i]; p$j <- map[p$j]
  dates <- seq(as.Date(if (period == periods[1]) "2026-01-01" else "2026-08-01"), as.Date("2026-08-31"), by = "day")
  calendar <- as.numeric(dates)
  day_id <- match(as.Date(d$acq_date), dates)
  tod <- as.integer(substr(d$acq_time, 1, 2)) + as.integer(substr(d$acq_time, 3, 4)) / 60
  observed_time <- calendar[day_id] * 24 + tod
  spatial_bin <- 1L + (p$distance_km > 0.5) + (p$distance_km > 1)
  same_day <- day_id[p$i] == day_id[p$j]
  fixed_counts <- counter(abs(tod[p$i[same_day]] - tod[p$j[same_day]]), spatial_bin[same_day])
  cross <- p[!same_day, ]
  cross_bin <- spatial_bin[!same_day]
  observed_gap <- abs(observed_time[p$i] - observed_time[p$j])
  observed <- counter(observed_gap, spatial_bin)
  primary <- p$distance_km <= 1 & observed_gap >= 1 & observed_gap <= 24
  day_a <- day_id[p$i[primary]]; day_b <- day_id[p$j[primary]]
  influence <- vapply(seq_along(dates), function(j) sum(day_a == j | day_b == j), integer(1))
  daily_checks[[length(daily_checks) + 1L]] <- data.frame(subset = lab, period = period,
    date = dates, records = tabulate(day_id, nbins = length(dates)),
    pairs_touching_date = influence, pairs_remaining_without_date = observed[5] - influence)
  for (scheme in schemes) {
    case_index <- case_index + 1L
    key <- paste(lab, period, scheme, sep = " | ")
    month <- format(dates, "%Y-%m")
    block <- (as.integer(format(dates, "%d")) - 1L) %/% 7L + 1L
    strata <- if (scheme == schemes[1]) factor(month) else interaction(month, block, drop = TRUE)
    groups <- split(seq_along(dates), strata)
    seed <- 2026091900L + case_index
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    simulations <- matrix(NA_integer_, B, 9)
    checked <- list()
    cat("Starting", key, "\n"); flush.console()
    for (b in seq_len(B)) {
      assigned <- calendar
      for (g in groups) assigned[g] <- calendar[g[sample.int(length(g))]]
      times <- assigned[day_id] * 24 + tod
      gap <- abs(times[cross$i] - times[cross$j])
      simulations[b, ] <- fixed_counts + counter(gap, cross_bin)
      if (b <= 3L) {
        all_gap <- abs(times[p$i] - times[p$j])
        direct <- vapply(seq_len(9), function(k) sum(p$distance_km <= thresholds$distance_km[k] &
          all_gap >= 1 & all_gap <= thresholds$max_gap_hours[k]), integer(1))
        stopifnot(identical(direct, simulations[b, ]), all(abs((times - observed_time) %% 24) < 1e-8))
        for (g in groups) stopifnot(identical(sort(assigned[g]), sort(calendar[g])))
        checked[[b]] <- assigned
      }
    }
    exceed <- colSums(sweep(simulations, 2, observed, ">="))
    result <- data.frame(subset = lab, period = period, scheme = scheme, thresholds,
      records = nrow(d), observed_pairs = observed, fixed_within_day_pairs = fixed_counts,
      observed_cross_day_pairs = observed - fixed_counts, null_mean = colMeans(simulations),
      null_q025 = apply(simulations, 2, quantile, .025, names = FALSE),
      null_q975 = apply(simulations, 2, quantile, .975, names = FALSE),
      null_max = apply(simulations, 2, max), observed_to_mean = observed / colMeans(simulations),
      exceedances = exceed, p_mc = (exceed + 1) / (B + 1), seed = seed)
    all_results[[length(all_results) + 1L]] <- result
    cases[[key]] <- list(summary = result, simulations = simulations, dates = dates,
      grouping = strata, checked_assignments = checked, source_file = d$source_file,
      source_row = d$source_row, day_id = day_id, hour_of_day = tod)
    cat("Completed", key, "\n"); flush.console()
  }
}
summary <- do.call(rbind, all_results)
stopifnot(nrow(summary) == 72L)
joint_p <- p.adjust(c(old$summary$p_mc, summary$p_mc), method = "holm")
summary$p_holm_144 <- tail(joint_p, 72)
summary$p_holm_day72 <- p.adjust(summary$p_mc, "holm")
old_adjusted <- old$summary
old_adjusted$p_holm_144 <- head(joint_p, 72)
daily <- do.call(rbind, daily_checks)
stopifnot(identical(hashes, tools::md5sum(inputs)))
dir.create(out)
saveRDS(list(summary = summary, earlier_tests = old_adjusted, cases = cases, daily = daily,
  settings = list(B = B, thresholds = thresholds, minimum_gap_hours = 1, timezone = "UTC",
    calendar_empty_days_included = TRUE, family_size = 144L, seed_base = 2026091900L,
    null = "Whole daily observed spatial patterns are exchangeable among calendar dates within strata"),
  source_md5 = hashes), file.path(out, "day_block_results.rds"))
write.csv(summary, file.path(out, "day_block_results.csv"), row.names = FALSE)
write.csv(daily, file.path(out, "daily_pair_influence.csv"), row.names = FALSE)
primary <- subset(summary, distance_km == 1 & max_gap_hours == 24)
display <- rbind(data.frame(model = "Individual records", subset = old_adjusted$subset,
    period = old_adjusted$period, scheme = old_adjusted$scheme, distance_km = old_adjusted$distance_km,
    max_gap_hours = old_adjusted$max_gap_hours, ratio = old_adjusted$observed_to_null_mean, p = old_adjusted$p_holm_144),
  data.frame(model = "Whole daily patterns", subset = summary$subset, period = summary$period,
    scheme = summary$scheme, distance_km = summary$distance_km, max_gap_hours = summary$max_gap_hours,
    ratio = summary$observed_to_mean, p = summary$p_holm_144))
display <- subset(display, distance_km == 1 & max_gap_hours == 24 & period == "January-August")
display$scheme <- factor(display$scheme, levels = rev(c("Monthly strata", "Within-month weekly strata", schemes)))
fig <- ggplot(display, aes(ratio, scheme, colour = model)) + geom_vline(xintercept = 1, colour = "#94A3B8") +
  geom_point(size = 3) + geom_text(aes(label = sprintf("Holm p = %.4f", p)), nudge_y = .2, size = 3.5) +
  facet_wrap(~subset) + scale_colour_manual(values = c("#D55E00", "#0072B2")) +
  scale_x_continuous(limits = c(0.65, 5.25), breaks = 1:5) +
  labs(x = "Observed pairs / permutation mean", y = NULL, colour = NULL,
    title = "Results depend on what the null preserves",
    subtitle = "January-August 2026, pairs within 1 km and 1-24 hours",
    caption = "Holm adjustment covers all 144 old and new comparisons. Ratios are pair enrichment, not fire-risk ratios.") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "null_model_comparison.png"), fig, width = 13, height = 6, dpi = 160, bg = "white")
capture.output({ print(primary, row.names = FALSE); print(sessionInfo()) }, file = file.path(out, "checks.txt"))
print(primary[, c("subset", "period", "scheme", "observed_pairs", "fixed_within_day_pairs", "null_mean", "p_mc", "p_holm_144")], row.names = FALSE)
