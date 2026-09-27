# Conditional space-time association of recorded MODIS detections.
# Usage: Rscript conditional_permutation.R <Ketapang data directory>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the Ketapang data directory.")
root <- normalizePath(args[1], mustWork = TRUE)
out <- file.path(root, "permutation_tests")
if (dir.exists(out)) stop("Output directory exists; review before rerunning.")
input <- file.path(root, "temporal_checks", "temporal_diagnostics.rds")
prior_hash <- tools::md5sum(input)
prior <- readRDS(input)
before <- tools::md5sum(file.path(root, basename(names(prior$source_md5))))
stopifnot(identical(unname(before), unname(prior$source_md5)))
# Freeze settings before the randomisation results are inspected.
B <- 9999L
alpha <- 0.05
base_seed <- 2026091700L
thresholds <- expand.grid(max_gap_hours = c(6, 24, 72), distance_km = c(0.5, 1, 2))
minimum_gap_hours <- 1
schemes <- c("Monthly strata", "Within-month weekly strata")
periods <- c("January-August", "August")
subset_labels <- c("Confidence >=30", "Confidence >=80")
all_results <- list()
summary_frames <- group_frames <- plot_frames <- list()
case_number <- 0L

make_counter <- function(pairs) {
  # Disjoint distance bins are accumulated to the nine nested thresholds.
  sb <- 1L + (pairs$distance_km > 0.5) + (pairs$distance_km > 1)
  force(sb); force(pairs)
  function(times) {
    gap <- abs(times[pairs$i] - times[pairs$j]) / 3600
    valid <- gap >= minimum_gap_hours & gap <= 72
    tb <- 1L + (gap[valid] > 6) + (gap[valid] > 24)
    cells <- matrix(tabulate((sb[valid] - 1L) * 3L + tb, nbins = 9L), nrow = 3L)
    cells[2, ] <- cells[2, ] + cells[1, ]
    cells[3, ] <- cells[3, ] + cells[2, ]
    cells[, 2] <- cells[, 2] + cells[, 1]
    cells[, 3] <- cells[, 3] + cells[, 2]
    as.integer(cells)
  }
}
direct_counter <- function(times, pairs) {
  gap <- abs(times[pairs$i] - times[pairs$j]) / 3600
  vapply(seq_len(nrow(thresholds)), function(k) sum(pairs$distance_km <= thresholds$distance_km[k] &
    gap >= minimum_gap_hours & gap <= thresholds$max_gap_hours[k]), integer(1))
}
for (lab in subset_labels) {
  full <- prior$subsets[[lab]]
  full_file <- file.path(root, if (lab == subset_labels[1]) "ketapang_modis_2026_conf30.csv" else "ketapang_modis_2026_conf80.csv")
  original <- read.csv(full_file, colClasses = "character")
  stopifnot(identical(original$source_file, full$records$source_file),
    identical(original$source_row, full$records$source_row), all(original$daynight %in% c("D", "N")))
  for (period in periods) {
    keep <- if (period == periods[1]) rep(TRUE, nrow(original)) else substr(original$acq_date, 1, 7) == "2026-08"
    d <- original[keep, ]
    mapping <- rep(NA_integer_, nrow(original)); mapping[keep] <- seq_len(sum(keep))
    p <- full$pairs_within_2km
    p <- p[keep[p$i] & keep[p$j], c("i", "j", "distance_km")]
    p$i <- mapping[p$i]; p$j <- mapping[p$j]
    stopifnot(!anyNA(p), all(p$i != p$j))
    times <- as.numeric(as.POSIXct(paste(d$acq_date, d$acq_time), format = "%Y-%m-%d %H%M", tz = "UTC"))
    stopifnot(!anyNA(times))
    counter <- make_counter(p)
    observed <- counter(times)
    stopifnot(identical(observed, direct_counter(times, p)))
    if (period == periods[1]) {
      previous_counts <- prior$sensitivity[prior$sensitivity$subset == lab &
        prior$sensitivity$rule == "At least 1 hour earlier", ]
      for (k in seq_len(9L)) {
        hit <- previous_counts$distance_km == thresholds$distance_km[k] &
          previous_counts$max_gap_hours == thresholds$max_gap_hours[k]
        stopifnot(sum(hit) == 1L, observed[k] == previous_counts$pair_count[hit])
      }
    }
    for (scheme in schemes) {
      case_number <- case_number + 1L
      key <- paste(lab, period, scheme, sep = " | ")
      month <- substr(d$acq_date, 1, 7)
      week <- (as.integer(substr(d$acq_date, 9, 10)) - 1L) %/% 7L + 1L
      grouping <- if (scheme == schemes[1]) {
        interaction(month, d$satellite, d$daynight, d$data_status, drop = TRUE)
      } else {
        interaction(month, week, d$satellite, d$daynight, d$data_status, drop = TRUE)
      }
      groups <- split(seq_len(nrow(d)), grouping)
      active <- groups[vapply(groups, function(g) length(unique(times[g])) > 1L, logical(1))]
      details <- do.call(rbind, lapply(names(groups), function(gname) {
        g <- groups[[gname]]
        data.frame(case = key, group = gname, records = length(g), distinct_times = length(unique(times[g])))
      }))
      group_frames[[length(group_frames) + 1L]] <- details
      seed <- base_seed + case_number
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      simulations <- matrix(NA_integer_, nrow = B, ncol = 9L)
      checked <- list()
      cat("Starting", key, "with", length(active), "movable strata and", B, "permutations.\n")
      flush.console()
      for (b in seq_len(B)) {
        permuted <- times
        for (g in active) permuted[g] <- times[g[sample.int(length(g))]]
        simulations[b, ] <- counter(permuted)
        if (b <= 3L) {
          stopifnot(identical(simulations[b, ], direct_counter(permuted, p)))
          for (g in groups) stopifnot(identical(sort(permuted[g]), sort(times[g])))
          checked[[b]] <- permuted
        }
        if (b %% 2500L == 0L) {
          cat(" ", key, ":", b, "permutations completed.\n"); flush.console()
        }
      }
      exceed <- colSums(sweep(simulations, 2, observed, ">="))
      p_mc <- (exceed + 1) / (B + 1)
      lo_hi <- t(vapply(exceed, function(k) stats::binom.test(k, B)$conf.int, numeric(2)))
      result <- data.frame(subset = lab, period = period, scheme = scheme, thresholds,
        records = nrow(d), observed_pairs = observed, permutation_mean = colMeans(simulations),
        permutation_sd = apply(simulations, 2, sd),
        null_q025 = apply(simulations, 2, quantile, probs = 0.025, names = FALSE),
        null_q975 = apply(simulations, 2, quantile, probs = 0.975, names = FALSE),
        null_maximum = apply(simulations, 2, max),
        observed_to_null_mean = observed / colMeans(simulations),
        excess_pairs_vs_null_mean = observed - colMeans(simulations),
        exceedances = exceed, permutations = B, p_mc = p_mc,
        tail_probability_lower_mc95 = lo_hi[, 1], tail_probability_upper_mc95 = lo_hi[, 2],
        seed = seed, strata = length(groups), movable_strata = length(active),
        records_in_movable_strata = sum(lengths(active)))
      summary_frames[[length(summary_frames) + 1L]] <- result
      all_results[[key]] <- list(summary = result, simulations = simulations,
        observed_pairs = observed, grouping = grouping,
        source_file = d$source_file, source_row = d$source_row,
        times_utc_seconds = times, checked_permutations = checked)
      primary <- which(thresholds$distance_km == 1 & thresholds$max_gap_hours == 24)
      plot_frames[[length(plot_frames) + 1L]] <- data.frame(subset = lab, period = period, scheme = scheme,
        close_pairs = simulations[, primary])
    }
  }
}
summary <- do.call(rbind, summary_frames)
stopifnot(nrow(summary) == 72L, all(summary$p_mc > 0), all(summary$p_mc <= 1))
# One family across both null schemes, both time windows, both confidence subsets,
# and all nine threshold settings. No significance-driven choice of settings.
summary$p_holm_72 <- p.adjust(summary$p_mc, method = "holm")
summary$reject_at_0_05 <- summary$p_holm_72 <= alpha
groups <- do.call(rbind, group_frames)
sim_plot <- do.call(rbind, plot_frames)
primary <- summary[summary$distance_km == 1 & summary$max_gap_hours == 24, ]
primary$panel <- paste(primary$subset, primary$period, sep = "\n")
sim_plot$panel <- paste(sim_plot$subset, sim_plot$period, sep = "\n")
p1 <- ggplot(sim_plot, aes(close_pairs)) + geom_histogram(bins = 35, fill = "#0072B2", colour = "white") +
  geom_vline(data = primary, aes(xintercept = observed_pairs), colour = "#D55E00", linewidth = 0.8) +
  facet_grid(scheme ~ panel, scales = "free_x") +
  labs(title = "Conditional time-permutation reference distributions",
    subtitle = "Pairs within 1 km and 1-24 hours | Orange line: observed count | 9,999 permutations per panel",
    x = "Number of qualifying unordered detection pairs", y = "Permutation frequency",
    caption = paste("Fixed locations and within-stratum timestamps. Strata preserve satellite, day/night and processing status, plus month or month/day block.",
      "These are model-conditional reference distributions, not confidence intervals for independent fires. Full-period and August datasets overlap.",
      "Monthly blocks allow within-month time reassignment. Weekly blocks use days 1-7, 8-14, 15-21, 22-28 and 29-end within each month.", sep = "\n")) +
  theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 9), plot.caption = element_text(hjust = 0, size = 8), plot.margin = margin(16, 16, 16, 16))
summary$panel <- paste(summary$subset, summary$period, sep = "\n")
p2 <- ggplot(summary, aes(factor(distance_km), factor(max_gap_hours), fill = p_holm_72)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.4f", p_holm_72)), size = 2.7) +
  facet_grid(scheme ~ panel) + scale_fill_gradient(low = "#5DCAB5", high = "#F1F5F9", limits = c(0, 1), name = "Holm p") +
  labs(title = "Multiplicity-adjusted conditional permutation results",
    subtitle = "One-sided excess-pair tests | One Holm family of 72 comparisons | All time gaps are at least 1 hour",
    x = "Maximum spatial separation (km)", y = "Maximum time separation (hours)",
    caption = paste("Significance level 0.05. Adjustments cover the declared comparisons, not all earlier exploratory choices or model misspecification.",
      "Per-record time exchangeability within strata is required. Shared overpasses, cloud cover and repeated fires can violate ecological interpretations.", sep = "\n")) +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 9), plot.caption = element_text(hjust = 0, size = 8), plot.margin = margin(16, 16, 16, 16))
stopifnot(identical(before, tools::md5sum(names(before))), identical(prior_hash, tools::md5sum(input)))
dir.create(out)
saveRDS(list(summary = summary, cases = all_results, strata = groups,
  source_md5 = before, prior_md5 = prior_hash, settings = list(B = B, alpha = alpha,
    base_seed = base_seed, thresholds = thresholds, minimum_gap_hours = minimum_gap_hours,
    schemes = schemes, periods = periods, alternatives = "greater", family_size = 72L,
    adjustment = "Holm", event_unit = "MODIS detection record", timezone = "UTC",
    null = "Conditional exchangeability of timestamps among fixed record locations within each stratum")),
  file.path(out, "conditional_permutation_results.rds"))
ggsave(file.path(out, "primary_permutation_distributions.png"), p1, width = 15, height = 8, dpi = 160, bg = "white")
ggsave(file.path(out, "adjusted_p_value_grid.png"), p2, width = 15, height = 8, dpi = 160, bg = "white")
capture.output({
  cat("CONDITIONAL SPACE-TIME PERMUTATION ANALYSIS\n\n")
  print(summary, row.names = FALSE, digits = 7)
  cat("\nPRIMARY DISPLAY: 1 KM AND 1-24 HOURS\n")
  print(primary, row.names = FALSE, digits = 7)
  cat("\nConditional null: within each stratum, timestamp labels are exchangeable over the observed locations.\n")
  cat("Statistic: unordered pairs within the spatial threshold with elapsed time from 1 hour through the temporal threshold, inclusive.\n")
  cat("The lower time cutoff is an explicit adaptation of a Knox-type close-pair statistic. The surveillance::knox default is not used.\n")
  cat("All source detections are retained in their specified period. Only time labels change in simulations.\n")
  cat("Monthly strata: month x satellite x day/night x processing status.\n")
  cat("Weekly sensitivity adds within-month blocks 1-7,8-14,15-21,22-28,29-end. Short final blocks remain short.\n")
  cat("One-sided Monte Carlo p = (1 + number of simulated statistics >= observed)/(9999 + 1). Ties count as exceedances.\n")
  cat("Minimum attainable p is 0.0001, never zero. Binomial intervals describe Monte Carlo tail uncertainty, not effect uncertainty.\n")
  cat("Holm adjustment uses all 72 declared comparisons together, regardless of result, at alpha 0.05.\n")
  cat("Thresholds and periods follow earlier exploratory inspection. The results are exploratory, not prospectively confirmatory.\n")
  cat("The observed/null-mean ratio describes pair-count enrichment, not relative fire risk or an independent-event ratio.\n")
  cat("Conditioning on locations preserves spatial inhomogeneity and the study-window geometry. No new K-function or edge-unbiased intensity estimate is claimed.\n")
  cat("This tests time-location association, not a pure spatial inhomogeneous Poisson null, causal contagion, or independent fires.\n")
  cat("Spatially changing cloud cover, observation effort, overpass groups and persistence can invalidate mechanistic interpretations.\n")
  cat("Stratification addresses selected marginal variation but does not demonstrate exchangeability. Pair observations are not treated as independent Bernoulli trials.\n")
  cat("Checks: prior source hashes, observed counts versus earlier diagnostics, direct versus binned counters, and within-stratum timestamp multisets for first three permutations.\n")
  cat("\nStratum sizes:\n"); print(groups, row.names = FALSE)
  cat("\nSource MD5 unchanged:\n"); print(before)
  cat("\nPrior temporal diagnostic MD5 unchanged:\n"); print(prior_hash)
  cat("\nReferences:\nhttps://search.r-project.org/CRAN/refmans/surveillance/html/knox.html\n")
  cat("https://arxiv.org/abs/1603.05766\nhttps://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html\n")
  print(sessionInfo())
}, file = file.path(out, "checks.txt"))
print(primary[, c("subset", "period", "scheme", "observed_pairs", "permutation_mean", "observed_to_null_mean", "p_mc", "p_holm_72")], row.names = FALSE)
cat("Total rejected among 72 declared comparisons:", sum(summary$reject_at_0_05), "\n")
cat("Saved results to", out, "\n")
