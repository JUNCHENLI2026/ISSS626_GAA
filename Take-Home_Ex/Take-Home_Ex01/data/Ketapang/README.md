# Ketapang geographic extraction

Prepared on 2026-09-16. This is a geographic subset, not a final cleaned analytical dataset or a count of independent forest fires. Original downloads were not modified.

## Navigation and current status

The extraction is followed by exploratory spatial diagnostics, individual-record time-label permutations (17 September 2026), whole-day permutations and spatial recurrence (19 September 2026), and lower-lag/three-day-block robustness checks with common-reference geographic screening (23 September 2026). All original detection subsets remain available; diagnostic flags are not exclusion rules. Current inference uses one Holm family across 192 comparisons, with adjusted p = 0.0192. The historical 144-family p = 0.0144 remains in earlier saved objects. Inference concerns recorded-detection time-location association under specified exchangeability models, not independent physical-fire interaction.

The final stage is `robustness_extension.R`, writing `robustness_extension/robustness_results.rds` and complete CSV tables. It runs 9,999 permutations for each of 12 cases and evaluates lower lags 0/1/3/6 hours at 1 km and an upper lag of 24 hours. Three-day patterns retain internal order; incomplete fragments remain fixed. All 48 comparisons enter the correction, including four repeated settings. `verify_robustness.R` independently checks counts, sampled assignments, all final adjustments and reference-point grid screening. The original data and historical analytical outputs remain unchanged.

| Question | Outputs | Reproduction script |
|---|---|---|
| Where are the records, and when were they observed? | `figures/` | `plot_data_checks.R` |
| How does fixed-bandwidth smoothing change the maps? | `learning/` (earlier technical outputs) | `learning/kde_learning.R` |
| What do the ordinary spatial summaries show descriptively? | `second_order/` | `second_order_diagnostics.R` |
| How sensitive are reweighted summaries to illustrative bandwidths? | `inhomogeneous_checks/` | `inhomogeneous_diagnostics.R` |
| Are high reciprocal weights traceable to the original records? | `weight_audit/` | `weight_leverage_audit.R` |
| What do two bandwidth criteria suggest, and are weights stable? | `bandwidth_checks/` | `bandwidth_selector_checks.R` |
| Are there earlier nearby detections, or strong temporal concentration? | `temporal_checks/` | `temporal_repeat_checks.R` |
| Does the descriptive spatial shape change with the time window? | `time_window_checks/` | `time_window_comparison.R` |
| Do spatial L summaries and weight diagnostics change with the time window? | `time_window_second_order/` | `time_window_second_order.R` |
| Are close space-time pairs enriched under conditional time reassignment? | `permutation_tests/` | `conditional_permutation.R` |
| Does association persist when complete daily patterns remain intact? | `day_block_tests/` | `day_block_permutation.R` |
| Where do records recur across dates, and how does this depend on spatial support? | `spatial_recurrence/` | `spatial_recurrence.R` |
| Do direct pair counts and grid assignments reproduce the new calculations? | Console verification | `verify_revision.R` |

The later scripts accept this Ketapang directory as their first argument and protect existing output directories. Where a previous result is required, scripts verify its source checksums before continuing. Individual `checks.txt` files record settings and numerical verification; RDS files preserve structured calculations. New diagnostics do not replace the source data or earlier results. The following sections document each stage and its limitations.

## Revision dated 19 September 2026

`day_block_permutation.R` reads the geographic subsets, temporal diagnostics and original individual-record permutation object. It applies one common date mapping to every detection on a UTC date. The alternatives exchange complete daily patterns within months or within days 1-7, 8-14, 15-21, 22-28 and 29-end. Calendar dates without retained records are included as empty recorded patterns, not verified observation opportunities. The full-period calendar has 243 dates; the August calendar has 31. All within-date arrangements remain fixed. Date exchangeability, UTC midnight boundaries and cross-date dependence remain assumptions or limitations.

The resulting `day_block_tests/day_block_results.rds` contains 72 new summaries, 9,999 simulated statistics per case, the first three audited date mappings per case, source identifiers and checksums, and the original 72 summaries with a recalculated `p_holm_144` field. The original permutation RDS remains unchanged. All 144 adjusted values are 0.0144. Among the new comparisons, 69 raw p-values are 0.0001, one is 0.0003 and two are 0.0004. The primary main-subset count is 3,597, including 1,946 fixed within-date pairs; the stricter whole-day null mean is 2,937.63. `daily_pair_influence.csv` is an arithmetic influence diagnostic, not a set of deletion-based significance tests.

`spatial_recurrence.R` assigns detections to 5 and 10 km UTM grids, each with an unshifted and half-cell-shifted origin, for both confidence subsets. It records cell assignments, detection totals, distinct UTC dates and months, and busiest-date shares. The saved object contains the primary-grid ranked cell list and all eight sensitivity settings. `subdistrict_context.csv` summarizes assignments to 20 unions derived from the supplied village attributes. Point-level modal subdistrict names label cells; they do not imply exclusive cell ownership. These summaries are descriptive and have neither a valid-observation denominator nor local significance tests.

`verify_revision.R` recomputes all observed threshold counts from an independent full distance matrix, checks stored within-stratum daily permutations, recalculates Monte Carlo p-values and the expanded Holm adjustment, and reconstructs grid assignments and date counts. It checks source hashes and stops on any mismatch. It verifies computation, not the scientific validity of exchangeability or the identification of independent fires.

For a complete build beginning with the original national boundary and NASA downloads, use `../../reproduce.R` as documented in the exercise-level README. An optional second argument to `filter_ketapang.R` selects a new extraction destination; the original one-argument usage is unchanged. The full rebuild writes into a new directory and does not overwrite original project outputs.

## Files

- `ketapang_modis_2026_jan_aug.csv`: 4,738 MODIS detection records inside the supplied Ketapang boundary. Original field text, including leading zeros in `acq_time`, is preserved. `source_file`, `source_row` (one-based data row, excluding the header) and `data_status` identify provenance. The absent NRT `type` field remains blank, not zero.
- `ketapang_boundary.gpkg`: `study_window` is the union of the selected village polygons; `villages` contains all 262 selected source features. Both layers are stored in WGS 84 (EPSG:4326).
- `monthly_record_counts.csv`: record counts by month and processing status, including zero combinations.
- `filter_checks.csv`: extraction diagnostics.
- `source_checksums.csv`: source-file MD5 checksums for reproducibility, not a security guarantee.
- `filter_ketapang.R`: extraction helper; requires R, sf and dplyr. Run with the parent data directory as its first argument. It refuses to overwrite an existing Ketapang output folder.
- `session-info.txt`: R and package versions used.

## Sources

The administrative data are `Batas_Wilayah_KelurahanDesa_10K_AR.*` under `../Indonesia-geospasial/`, downloaded by the user from https://www.indonesia-geospasial.com/2023/05/download-shapefile-batas-administrasi.html . Its XML metadata creation date is 2023-05-28, despite the download page's 2024 title; this is not evidence of an authoritative 2024 boundary reference date.

Fire detections are the user's NASA FIRMS downloads under `../NASA FIRMS/`: `fire_archive_M-C61_803220.csv` and `fire_nrt_M-C61_803220.csv`. Source: https://firms.modaps.eosdis.nasa.gov/download/ . Product explanation: https://firms.modaps.eosdis.nasa.gov/download/Readme.txt . Supplied NASA confirmation emails establish request 803220 and download availability on 10 September 2026, for MODIS C6.1, Indonesia, CSV format and 1 January-31 August 2026. The request body gives 2026-09-10 13:49:56 without a timezone; no timezone is inferred. The request identifier matches both CSV filenames. These are request and availability dates, not independent proof of local file-transfer time. Personal email details and screenshots are excluded from the website.

## Selection and verification

1. Selected village records with `WADMPR = Kalimantan Barat` and `WADMKK = Ketapang`. Their populated `KDPKAB` values are `61.04` (the code convention in this source, not assumed to be a BPS code).
2. Removed the Z dimension for two-dimensional topology operations and transformed a working copy to EPSG:32749. No invalid village geometries were detected. `st_make_valid` was applied defensively, and all resulting geometries were valid and non-empty.
3. Unioned all 262 village geometries without simplifying, buffering, filling holes, or discarding disconnected components. This is a boundary derived from the supplied village data, not an independently verified official regency polygon. Geometry validity alone does not prove that coverage is complete or current.
4. Matched all 46,449 source detections against that union using intersection, including exact boundary points. Selected 4,738 and excluded 41,711 from this output only. There were zero selected points lying only on the boundary.
5. Preserved all confidence levels and satellite/day-night categories. No type-based filtering, temporal thinning, deduplication, or event-merging was applied. No duplicate latitude/longitude/date/time/satellite/instrument keys were found in the selection. Repeated observations of one physical fire may still occur.
6. The unnamed West Kalimantan polygon has zero area overlap with the derived Ketapang boundary and contains zero input detections. It was not assigned to Ketapang.
7. Reopened the exported CSV and GeoPackage to verify record count, time strings and selected-point inclusion.

## Temporal coverage and limitations

The requested observation window is 2026-01-01 through 2026-08-31. Selected detections actually run from 2026-01-14 through 2026-08-31. Absence of a record on a date does not establish absence of fire or complete satellite coverage.

| Month | Standard | NRT | Total detections |
|---|---:|---:|---:|
| 2026-01 | 26 | 0 | 26 |
| 2026-02 | 22 | 0 | 22 |
| 2026-03 | 82 | 0 | 82 |
| 2026-04 | 11 | 0 | 11 |
| 2026-05 | 0 | 5 | 5 |
| 2026-06 | 0 | 50 | 50 |
| 2026-07 | 0 | 255 | 255 |
| 2026-08 | 0 | 4287 | 4287 |
| Total | 141 | 4597 | 4738 |

Standard and NRT data differ in processing status. All 141 selected standard records have `type = 0`; the 4,597 NRT records have no supplied type classification. Do not impute that missing classification as zero or assume all selected records are confirmed forest fires. Monthly counts are detection counts, not area burned, fire risk or independent event counts.

## Additional quality checks

Read-only checks of the exported 4,738 records found zero repeated exact coordinate pairs and observations on 98 distinct dates. This does not rule out repeated detection of the same physical fire at slightly different coordinates or times.

| Confidence range (diagnostic bins only) | Standard | NRT | Total |
|---|---:|---:|---:|
| Below 30 | 7 | 259 | 266 |
| 30 to below 80 | 99 | 2647 | 2746 |
| 80 and above | 35 | 1691 | 1726 |

The geographic subset retains all 4,738 records. The following separately saved subsets apply the confidence thresholds subsequently confirmed by the user; no original records were overwritten.

## Confidence subsets

- `ketapang_modis_2026_conf30.csv`: main-analysis input, `confidence >= 30`, 4,472 detections (134 standard and 4,338 NRT).
- `ketapang_modis_2026_conf80.csv`: sensitivity-check input, `confidence >= 80`, 1,726 detections (35 standard and 1,691 NRT).
- `monthly_confidence_comparison.csv`: eight monthly rows comparing the unchanged geographic subset with both confidence subsets, including excluded counts and standard/NRT breakdowns. Counts refer to satellite detections, not independent fire incidents.

Filtering was performed numerically on `confidence` only. All original columns, values, order, source-row identifiers and four-character HHMM time strings are retained. Missing NRT `type` values remain blank. No thinning, type filtering, event merging, or further spatial changes were performed.

| Month | All detections | Confidence >=30 | Confidence >=80 |
|---|---:|---:|---:|
| 2026-01 | 26 | 23 | 8 |
| 2026-02 | 22 | 22 | 4 |
| 2026-03 | 82 | 79 | 23 |
| 2026-04 | 11 | 10 | 0 |
| 2026-05 | 5 | 5 | 1 |
| 2026-06 | 50 | 49 | 5 |
| 2026-07 | 255 | 248 | 42 |
| 2026-08 | 4287 | 4036 | 1643 |
| Total | 4738 | 4472 | 1726 |

Both subsets have detections dated 2026-01-16 through 2026-08-31, within the unchanged study window of January 1 to August 31. They contain records on 95 and 67 distinct dates, respectively. April has no retained detections at the >=80 threshold; this does not imply there were no fires in April.

Verification reconciled all monthly counts, confirmed the >=80 data are a subset of the >=30 data, round-tripped every exported cell, and checked that the original geographic CSV's SHA256 was unchanged (`7afc972022d77fc153e44beb995f63b1ade98282232fe4c9a3ed4bdb791f4ec6`). CSV is plain text: when importing, treat `acq_time` as text to retain leading zeros.

## Descriptive figures

- `figures/ketapang_detection_locations.png`: side-by-side maps of the >=30 and >=80 confidence subsets, with identical geographic extent, point size and transparency. The saved study window is displayed without simplifying its geometry.
- `figures/ketapang_monthly_detections.png`: grouped monthly counts for all detections and the two nested subsets. The vertical axis starts at zero; every bar has a count label, including the zero for April in the >=80 subset. The series must not be summed.
- `plot_data_checks.R`: reproduces these figures using R, sf and ggplot2. Pass this directory as its first argument. Existing figures are protected from accidental replacement.

Every plotted point was checked against the saved boundary, and plotted monthly counts were reconciled to the three input CSV files. All five data inputs (the three detection CSVs, monthly comparison CSV and boundary GeoPackage) retained unchanged MD5 checksums during plotting. Both exported images were visually checked for legible labels and complete content.

These figures describe recorded detections only. They are not density estimates, tests of clustering, estimates of burned area, or evidence of causal effects. No report findings or statistical interpretations were generated in this step.

## Kernel intensity outputs

Separate exploratory outputs are saved under `learning/`:

- `kde_bandwidth_demo.png`: confidence >=30 at Gaussian bandwidths of 5, 10 and 20 km, using a common colour scale.
- `kde_confidence_demo.png`: the two confidence subsets at the same 10 km bandwidth and colour scale. Intensities are absolute, not divided by the number of records.
- `kde_learning.R`: reproducible calculation and plotting script.
- `kde_checks.txt`: numerical checks, source checksums and runtime versions.

These bandwidths are illustrative, not an optimised final selection. Calculations use EPSG:32749 with coordinates rescaled to kilometres, the full derived polygon window, a Gaussian kernel, Jones-Diggle edge correction and a 0.5 km computational grid. Legend units are detections per square kilometre over the January-August observation period, not probabilities, daily rates or independent-fire counts. Pixel size does not describe the satellite's measurement resolution. Display-only negative FFT round-off values, if any, are clipped to zero.

Checks confirm all input points are inside the window, polygon area is preserved when converting spatial objects, and each surface's grid integral matches its corresponding point count within numerical tolerance. The three data inputs retain unchanged checksums. No null model, significance envelope or clustering test has been fitted. See the [spatstat authors' intensity-estimation notes](https://spatstat.org/ECAS2019/notes/notes02.html) for method background.

## Exploratory second-order diagnostics

`second_order_diagnostics.R` produces the following files under `second_order/`:

- `observed_L_summary.png`: observed border-corrected L(r) - r summaries for the two nested confidence subsets.
- `border_support.png`: the proportion of observed points eligible as interior centres for border correction at each whole-kilometre distance from 0 to 10 km.
- `observed_summaries.rds`: K estimates, transformed curves, support counts and calculation metadata.
- `checks.txt`: input checks, method settings, checksums and runtime versions.

Run the script with this data directory as its sole argument. It refuses to overwrite an existing `second_order` output directory.

Calculations retain the full polygon window and all input points, using EPSG:32749 with coordinates in kilometres. The initial 0-10 km distance range is provisional, not an independently justified final testing range. Border correction is specified explicitly for computational tractability with the large point patterns and detailed boundary; its loss of eligible centres is reported separately. No coordinate simplification, temporal thinning or event merging is applied.

The zero line in the L(r) - r figure is the homogeneous Poisson theoretical reference, not a significance threshold. The homogeneous estimator assumes spatial stationarity; its suitability for these detections has not been established. Spatially varying intensity, observation opportunity and repeated detections of the same physical event remain unadjusted. The curves do not isolate interaction from these effects. No simulations, p-values, significance envelopes or inferential conclusions are provided.

Verification confirms that all points remain inside the window, the area conversion is consistent, the curves contain finite non-negative K estimates, and eligible-centre counts are positive and non-increasing over the diagnostic distances. Input checksums are unchanged. Method references: [Kest](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kest.html) and [Lest](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Lest.html).

## Intensity-reweighted sensitivity diagnostics

`inhomogeneous_diagnostics.R` adds six exploratory calculations: the two existing confidence subsets at Gaussian bandwidths of 5, 10 and 20 km. It reads the previous summary object to use the same distance grid and verifies that its input checksums match the current data. Existing outputs are protected from overwriting.

Files under `inhomogeneous_checks/`:

- `ordinary_vs_reweighted_L.png`: ordinary and intensity-reweighted L(r) - r curves, with shared axes across the two confidence subsets.
- `reweighted_bandwidth_sensitivity.png`: the same three reweighted curves without the ordinary curve, using a closer vertical scale. This is not an independent analysis.
- `intensity_reweighted_summaries.rds`: raw K estimates, pointwise intensity vectors, plotted values, provenance identifiers and settings.
- `checks.txt`: numerical checks, inverse-intensity weight diagnostics, source checksums and runtime versions.

Intensity is evaluated at each original point using Gaussian leave-one-out smoothing with Jones-Diggle edge correction and a 0.5 km computational edge grid. This excludes the focal point's own kernel contribution; it does not create an independent validation sample. Positive intensity values are passed directly to `Kinhom`, with border correction, `renormalise = TRUE` and `normpower = 2`. The latter is the normalization power recommended in the function documentation rather than its legacy default of 1. All settings remain exploratory, including the 0-10 km range and bandwidths. No final null model or significance test has been selected or run.

All six intensity vectors are finite and positive. No intensity floor, weight cap, data deletion or temporal thinning is applied. Raw K estimates are preserved, including near-zero negative round-off (minimum approximately -5.2e-12 km2). Only the square-root plotting input is bounded at zero, after verifying that any negative value is within the stated 1e-8 km2 numerical tolerance.

The 5 km calculations have maximum normalised inverse-intensity shares of approximately 45.1% (confidence >=30) and 42.7% (confidence >=80). These are shares of single-point reciprocal-intensity sums, not shares of all pair contributions or percentages of independent events. They flag weight concentration requiring investigation before inference. The reported reciprocal effective count is a weight-concentration diagnostic, not an effective number of independent fires.

Fitting intensity and computing the summary from the same observations makes bandwidth and estimation choices consequential. Intensity reweighting does not establish that a suitable second-order model holds, distinguish environmental variation from observation effects, or account for repeated detections of a physical event. Future fitted-null simulations would need to address intensity estimation and parameter selection. No report findings or claims of significant clustering or inhibition are supplied here. All source files and the earlier summary object retain unchanged checksums.

References: [Kinhom](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kinhom.html) and [density.ppp](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/density.ppp.html).

## Source and proximity audit of high reciprocal weights

`weight_leverage_audit.R` traces the five highest reciprocal-intensity weights in each of the six existing subset/bandwidth combinations. It does not refit the curves or exclude observations. Run it with this data directory as its sole argument; it protects an existing `weight_audit` directory from replacement.

Outputs under `weight_audit/`:

- `checks.txt`: source-row identifiers, original attributes, weight shares, distances to the study-window boundary, nearest retained records and nearest records in the all-confidence Ketapang subset.
- `record_audit.rds`: the structured 30-row audit and cumulative weight diagnostics.
- `reciprocal_weight_concentration.png`: cumulative reciprocal-intensity shares of the 20 highest-weight records, separately ranked for each setting.
- `flagged_record_locations.png`: the locations of the five highest-weight records at 5 km, overlaid on each unchanged confidence subset.

All original source columns in both retained subsets were compared against the supplied NASA downloads using the recorded file and one-based data-row identifiers. The comparisons match, and all retained points are inside the supplied window. This verifies transfer and geographic selection against these inputs, not the accuracy of the satellite classification or the administrative boundary itself.

At 5 km, the highest-weight record in the >=30 subset is NRT source row 39148, dated 2026-08-29, confidence 59. It is approximately 12.08 km from the window boundary and 17.27 km from its nearest retained record; the nearest distance remains 17.27 km when all confidence levels inside Ketapang are included. The highest-weight record in the >=80 subset is NRT source row 38442, dated 2026-08-29, confidence 86. Its nearest retained record is approximately 15.74 km away, whereas the all-confidence subset has a record approximately 1.41 km away. The closer record is excluded by the >=80 threshold.

Both leading 5 km records have no other retained record within 10 km. Their reciprocal-intensity shares therefore must not be presented as their shares of nearby K-function pair contributions. Normalization and border denominators also use intensity weights; a full influence calculation would require a separately defined diagnostic. All proximity searches here pool the observation period and are not time-constrained or searches outside Ketapang.

No mismatch established by this audit justifies record removal. No records were deleted, no weights were capped, and no final bandwidth or inferential model was selected. Input and prior-result checksums remain unchanged. The next methodological decision requires a justified intensity-estimation approach and sensitivity assessment, not selection of whichever curve gives a preferred conclusion.

## Exploratory bandwidth criteria and transfer checks

`bandwidth_selector_checks.R` compares likelihood cross-validation (`bw.ppl`) and the Cronie-van Lieshout criterion (`bw.CvL`) for the two unchanged confidence subsets. Run the script with this directory as its sole argument. It refuses to replace an existing `bandwidth_checks` output directory.

The initial search evaluates 20 log-spaced bandwidths from 1 to 40 km. If a best candidate is at an endpoint, that bound is expanded by a factor of four, with at most three coarse searches. An interior candidate is refined using 17 log-spaced values between its adjacent coarse candidates. These are explicit computational search choices, not scientifically established bandwidth limits. The output records endpoint status and the neighbouring evaluated bandwidths; an interior grid optimum is not proof of a unique continuous optimum.

The installed selectors use different estimators. Likelihood CV is explicitly run with Gaussian smoothing, Jones-Diggle edge correction, leave-one-out point intensity, a 0.5 km computational grid and `shortcut = FALSE`, retaining the integrated-intensity term. The installed CvL implementation instead uses inclusive point intensity and no edge correction, and ignores additional arguments. Their objective values therefore must not be compared as if they measure the same criterion.

For each selected candidate and its 80% and 120% bandwidths, a separate weight audit applies the same Gaussian leave-one-out, Jones-Diggle estimator used in the earlier intensity-reweighted diagnostics. Applying a CvL candidate to this different estimator is explicitly a transfer sensitivity check, not an equivalent optimisation. The 80% and 120% values are exploratory perturbations, not confidence limits. No arbitrary intensity floor, weight cap, point removal or temporal thinning is used.

Outputs under `bandwidth_checks/`:

- `bandwidth_objective_curves.png`: the gap from the best evaluated objective within each panel, with log-scaled bandwidth and a log1p vertical scale. Objective magnitudes are not comparable between panels.
- `candidate_weight_sensitivity.png`: the largest reciprocal-intensity share under the common downstream estimator.
- `bandwidth_diagnostics.rds`: complete selector objects, objective grids, candidate summaries, pointwise audit intensities, settings and input checksums.
- `checks.txt`: numerical results, explicit estimator differences, installed selector implementations and runtime versions.

These criteria do not distinguish repeated detections from separate physical fires, identify environmental causes of intensity variation, or validate the assumptions of a subsequent interaction test. Lower weight concentration alone is not a reason to select a bandwidth. No final bandwidth, null model, significance test or report conclusion is selected in this step. Numerical-grid convergence also remains a separate check before final use.

Method references: [likelihood cross-validation](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/bw.ppl.html) and [Cronie-van Lieshout criterion](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/bw.CvL.html).

The evaluated likelihood-CV candidates are 1.2143 km (confidence >=30) and 1.6647 km (confidence >=80); CvL candidates are 13.0986 km and 17.1066 km, respectively. All four are interior optima of the evaluated grids, with neighbouring grid values recorded in the checks. Under the common leave-one-out estimator, the likelihood-CV candidates concentrate approximately 100% of the reciprocal-intensity sum on one record in each subset. The CvL candidates have maximum shares of approximately 2.0973% and 1.4846%. Their 80%-120% bandwidth sensitivity checks yield maximum shares of approximately 1.79%-2.10% and 1.13%-2.07%, respectively.

These results flag severe reciprocal-weight concentration at the small likelihood-CV bandwidths. They do not establish that CvL is a final appropriate choice: the estimator-transfer caveat, time pooling, observation process and repeated detections remain unresolved. All twelve intensity vectors are finite and positive, showing that positivity alone is insufficient as a stability check. Source and previous-result checksums are unchanged.

## Temporal-repeat and acquisition-proximity checks

`temporal_repeat_checks.R` adds exploratory space-time proximity diagnostics without merging, thinning or removing records. Run it with this directory as its sole argument; it protects an existing `temporal_checks` output directory. All acquisition dates and four-digit HHMM times are parsed and round-trip checked in UTC. Coordinates use EPSG:32749, rescaled to kilometres.

Outputs under `temporal_checks/`:

- `daily_detection_timing.png`: daily recorded detections by satellite, with separate vertical scales for the two confidence subsets. The dashed line marks the supplied standard/NRT file transition on 1 May, not a fitted change point.
- `space_time_proximity_sensitivity.png`: percentages and counts of records with an earlier nearby detection, requiring an elapsed gap of at least one hour.
- `temporal_diagnostics.rds`: original-row provenance, projected coordinates, all unordered spatial pairs within 2 km, acquisition-time gaps, record-level proximity flags, daily counts, diagnostic tables, settings and checksums.
- `checks.txt`: exact-key checks, acquisition-proximity counts, all 36 subset/distance/time/rule combinations, method caveats and runtime versions.

The distance limits are 0.5, 1 and 2 km; upper elapsed-time limits are 6, 24 and 72 hours. Two lower-gap rules are compared: strictly positive elapsed time and at least one hour. These are illustrative sensitivity thresholds, not calibrated rules for identifying the same physical fire. Each focal record is counted once per setting if any earlier record qualifies; pair counts are reported separately. Same-minute pairs do not have a strictly earlier endpoint and are excluded from both prior-record rules. Records sharing date, acquisition minute and satellite are grouped only as an acquisition proxy, not a verified orbit or granule.

Neither subset contains duplicated source-row identifiers, duplicated latitude/longitude/date/minute/satellite keys, or repeated numeric coordinates. This does not establish independent fire events. Within 2 km, the >=30 and >=80 subsets contain 3,335 and 1,449 same-minute, same-satellite unordered pairs, respectively. These are distinct neighbouring pixel detections, not exact duplicate rows. No within-2-km pair has a strictly positive gap below one hour in the supplied subsets, so the two lower-gap rules give identical counts here.

For the illustrative 1-km and 1-24-hour rule, 1,632 of 4,472 records (36.49%) in the >=30 subset and 644 of 1,726 records (37.31%) in the >=80 subset have at least one earlier nearby detection. These are proportions of retained records with a proximity match, not estimated percentages of duplicated fires or observations that should be deleted. The two subsets are nested and are not independent comparisons.

August contains 4,036 records (90.25%) and 1,643 records (95.19%), respectively. The highest recorded daily totals occur on 29 August UTC: 674 and 318 detections. Across the 243 supplied calendar days, the subsets have records on 95 and 67 dates. A date with zero retained records is not proof of no fire or a valid cloud-free observation. The data do not supply an observation-opportunity denominator, and the standard/NRT transition remains an additional limitation.

MODIS coordinates are fire-pixel centres, not known ignition locations. Pixel size varies: supplied scan dimensions range from 1.00 to 4.82 km in the main subset, and track dimensions from 1.00 to 2.00 km. Centre-distance tests therefore do not establish overlapping footprints. A continuing fire may be recorded repeatedly or spread; separate fires can also be nearby. Neighbouring locations alone cannot distinguish these possibilities. The one-hour rule also cannot prove that records belong to different satellite passes.

All searches are restricted to retained records inside Ketapang and the supplied January-August period. Earlier detections outside this geographic or temporal window can be missed, so a record without a match is not proven independent. Source data and prior bandwidth results remain unchanged. Internal checks cover unique unordered pairs, sampled direct distance/time calculations, monotone counts as thresholds widen, nested absolute counts and reconciled daily totals. An independent exhaustive distance check also verified every enumerated spatial neighbour pair and every focal flag under the 1-km/1-24-hour rule; exported daily-bar totals were reconciled to the inputs, and both figures were visually checked.

No event merging, temporal thinning, final time-window selection, simulation, significance test or independent-fire estimate is performed. Before inferential use, the analysis needs an explicit distinction between recorded-detection patterns and physical-fire events, plus consideration of temporal concentration and observation opportunity. A narrower time window chosen after viewing the data would be a sensitivity analysis, not a pre-specified design.

Reference: [NASA FIRMS MODIS attribute definitions](https://firms.modaps.eosdis.nasa.gov/descriptions/FIRMS_MODIS_Firehotspots.html).

## Time-window comparison at common smoothing settings

`time_window_comparison.R` compares January-August, January-July and August 2026 in each confidence subset. August is an exploratory sensitivity window inspected after observing temporal concentration, not a replacement for the primary January-August period. January-July is included because comparing the full period with its dominant August component alone provides limited information about differences between periods.

Outputs under `time_window_checks/`:

- `relative_spatial_distribution.png`: six maps of normalised spatial shape on a common colour scale, with record counts in each panel.
- `time_window_diagnostics.rds`: raw raster intensities, unit-integral relative densities, pixel masses, source-row identifiers, count tables, comparison metrics, settings and checksums.
- `checks.txt`: count reconciliation, numerical integrals, grid alignment, KDE linearity, mixture checks and runtime versions.

All six estimates use the same full Ketapang window, EPSG:32749 with kilometres, a Gaussian kernel, Jones-Diggle edge correction and a 0.5 km computational grid. The shared 10 km bandwidth is an illustrative setting retained from earlier diagnostics to hold smoothing constant. It is not a newly selected final bandwidth, and this comparison does not establish robustness at other smoothing scales. Raster estimates include all retained points; the leave-one-out setting used for previous pointwise weight checks does not apply to these pixel maps.

The raw intensity integral is checked against each period's record count. Raw estimates are retained, including any negligible negative FFT round-off. Only negative values within the specified 1e-8 tolerance are clipped before constructing nonnegative relative densities, which are normalised to integrate to one. This removes differences in total record count for a spatial-shape comparison; the resulting maps do not show absolute detection rates, independent-fire risk, burned area or exposure-adjusted activity. Calendar durations are not valid-observation denominators.

For two unit-sum pixel-mass vectors, the total-variation distance is `0.5 * sum(abs(p - q))`, and the overlap is `sum(pmin(p, q)) = 1 - TV`. These are descriptive, bandwidth-dependent map comparisons, not p-values or fractions of matching detection records. The full-period intensity is verified against the sum of the two disjoint-period intensities. Consequently, up to numerical error, `TV(full period, August) = January-July record fraction * TV(January-July, August)`. A high overlap between the full period and August is partly an algebraic consequence of August's dominance, not independent confirmation of temporal stability.

January-July contains only 436 records in the >=30 subset and 83 in the >=80 subset, compared with 4,036 and 1,643 in August. The disjoint-period comparison avoids shared records but does not establish independence, equal precision or equal observation opportunity. January-July also combines standard and NRT files, whereas August uses NRT records. No record thinning, event merging, coordinate modification, fitted null model, significance test or causal explanation is added. Source files and previous temporal-check results remain unchanged.

Method reference: [spatstat kernel intensity estimation](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/density.ppp.html).

At these fixed smoothing settings, the descriptive map-overlap results are:

| Confidence subset | January-August vs August | January-July vs August |
|---|---:|---:|
| >=30 | 97.08% | 70.05% |
| >=80 | 98.25% | 63.65% |

These percentages describe the overlap of normalised smoothed grid mass, not the percentage of repeated fires or matching observations. The direct full-period KDE and the sum of the disjoint-period KDEs agree to a maximum absolute error below 7e-16; source identifiers partition correctly, and all six normalised surfaces integrate to one. The small January-July samples, especially the 83 high-confidence records, require caution when interpreting differences. Neither overlap comparison is a significance test. Independent checks recomputed the metrics and source partitions, and the six-panel figure was visually verified.

## Spatial-summary sensitivity across time windows

`time_window_second_order.R` adds ordinary and intensity-reweighted spatial L summaries for January-August, January-July and August in both confidence subsets. These are exploratory learning diagnostics, not inferential findings or a final analysis specification. The script retains the complete study window and every record in each stated period. It checks input hashes and source-row identifiers before reusing the earlier full-period estimates unchanged.

Outputs under `time_window_second_order/`:

- `ordinary_L_by_time_window.png`: six ordinary border-corrected L(r)-r curves, with shared axes across confidence subsets.
- `reweighted_L_by_time_window.png`: twelve intensity-reweighted curves for the three periods, two confidence subsets and two illustrative bandwidths.
- `border_support_by_time_window.png`: absolute numbers of observed centres farther than each diagnostic distance from the window boundary.
- `time_window_second_order.rds`: raw estimates, pointwise intensities, source-row identifiers, plotted curves, weight diagnostics, support counts, settings and checksums.
- `checks.txt`: method settings, numerical checks, weight concentration, border support, provenance and runtime versions.

The shared distance grid is the earlier provisional 0-10 km grid. Both 10 and 20 km Gaussian bandwidths are held constant across periods as sensitivity settings, not selected final parameters. Each new period's intensity is estimated at its own points with leave-one-out smoothing, Jones-Diggle edge correction and a 0.5 km computational grid. These are pointwise detection intensities, not the unit-integral raster surfaces from the preceding map comparison. Positive numeric intensity vectors are supplied to `Kinhom`, using border correction, `renormalise = TRUE` and `normpower = 2`.

All twelve reweighted curves contain valid finite K values under the stated round-off tolerance. Raw values are preserved. Only negative round-off within 1e-8 of zero is bounded for the square-root display; more negative or nonfinite values would be reported as invalid and omitted from that display. No observations, inverse-intensity weights or intensities are capped or removed. Independent verification reproduced the plotted transformations and weight summaries and confirmed that the four cached full-period estimates are unchanged. All three figures were visually checked.

The January-July high-confidence subset contains 83 detections. Its maximum single-record reciprocal-intensity share is 17.58% at a 10 km bandwidth, compared with 3.85% at 20 km. At r = 10 km, only 55 observed centres qualify for border correction; within that eligible-centre denominator the maximum reciprocal share is 19.77% and 4.60%, respectively. These are weight-concentration diagnostics, not shares of K-function pair contributions or percentages of independent fires. Lower concentration at 20 km does not establish that it is the correct final bandwidth.

For comparison, the August high-confidence subset retains 960 eligible centres at r = 10 km, out of 1,643 records. Its maximum all-record reciprocal shares are 3.64% and 1.17% at the same two bandwidths. The January-July main subset retains 297 of 436 centres, with maximum all-record reciprocal shares of 6.52% and 2.08%. Changes in sample composition and weight concentration must be considered alongside any visual curve differences.

The zero line is a Poisson theoretical reference, not a simulation envelope or significance threshold. Ordinary curves do not adjust for spatially varying intensity. Reweighted curves depend on the intensity estimate and its assumptions; a negative segment does not by itself demonstrate inhibition, and a near-zero segment does not establish randomness. The same data estimate intensity and the summary. Nested periods, confidence filtering, repeated detections, temporal concentration and observation opportunity remain relevant limitations. No null model, simulation design, p-value, final bandwidth or formal clustering conclusion is selected here.

References: [inhomogeneous K](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kinhom.html) and [ordinary K](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kest.html).

## Conditional time-label permutation analysis

`conditional_permutation.R` takes this directory as its sole argument and protects an existing `permutation_tests/` directory. It uses the verified projected locations and pair lists from `temporal_checks/temporal_diagnostics.rds`. Original records, coordinates and timestamps are retained; reassignment occurs only in simulation copies.

The test counts unordered pairs within 0.5, 1 or 2 km and an elapsed interval from 1 hour to an upper limit of 6, 24 or 72 hours. It holds observed locations fixed and permutes timestamps within month, satellite, day/night and processing-status strata. A more restrictive sensitivity scheme adds within-month day blocks 1-7, 8-14, 15-21, 22-28 and 29-end. These are not ISO calendar weeks. Four datasets combine the two nested confidence thresholds with the full period or August alone.

Each of eight dataset/stratum cases uses 9,999 permutations, with seeds 2026091701-2026091708. The nine distance/time combinations share permutations within a case. The one-sided Monte Carlo p-value is `(1 + exceedances) / 10000`, including ties. Holm adjustment is applied across the complete family of 72 comparisons at alpha 0.05. Thresholds and time windows followed earlier data inspection, so the analysis is exploratory, not prospectively confirmatory.

Outputs under `permutation_tests/`:

- `conditional_permutation_results.rds`: all 72 results, all 79,992 simulated count vectors, grouping assignments, settings, source hashes, provenance and three saved timestamp permutations per case for verification.
- `primary_permutation_distributions.png`: observed counts and conditional reference distributions at 1 km and 1-24 hours.
- `adjusted_p_value_grid.png`: the declared threshold family and Holm-adjusted values.
- `checks.txt`: numerical results, assumptions, Monte Carlo uncertainty and runtime information.

At 1 km and 1-24 hours, the main full-period subset has 3,597 observed pairs. Null means are 1,091.8116 under monthly strata and 2,343.9909 under the more restrictive blocks: observed/null-mean ratios of 3.2945 and 1.5346. Every declared comparison has zero simulated exceedances, giving p = 0.0001 and Holm-adjusted p = 0.0072. Equal p-values reflect finite simulation resolution, not identical effect sizes. The saved two-sided 95% binomial intervals describe uncertainty in the permutation tail probability, not ecological effect uncertainty.

These results support significant time-location association of the recorded detections relative to the declared conditional models. They do not establish independent-fire clustering, causal spread or a pure spatial inhomogeneous-Poisson result. Record-level timestamp exchangeability may be undermined by shared overpasses, persistent fires or spatially varying observation opportunity; finer strata and a one-hour cutoff do not remove every such dependence. Spatial and temporal pair windows remain fixed across permutations, rather than using an edge-corrected K estimator.

Independent computational checks re-enumerated all spatial pairs from projected coordinates, reconciled the nine observed statistics and three stored permutations per case, checked within-stratum timestamp multisets, recomputed every Monte Carlo p-value, and reconstructed Holm adjustment manually. Source and prior-result hashes remained unchanged.

The technical report and revealjs executive summary are in the exercise root as `Take-Home_Ex01.qmd` and `executive_summary.qmd`. Their sources and renderable outputs are English; the report contains the complete permutation implementation, methodological assumptions and references. Rendering reads saved permutation results rather than rerunning all simulations.
