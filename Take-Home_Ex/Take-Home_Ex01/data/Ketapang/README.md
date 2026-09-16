# Ketapang geographic extraction

Prepared on 2026-09-16. This is a geographic subset, not a final cleaned analytical dataset or a count of independent forest fires. Original downloads were not modified.

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

Fire detections are the user's NASA FIRMS downloads under `../NASA FIRMS/`: `fire_archive_M-C61_803220.csv` and `fire_nrt_M-C61_803220.csv`. Source: https://firms.modaps.eosdis.nasa.gov/download/ . Product explanation: https://firms.modaps.eosdis.nasa.gov/download/Readme.txt . Consult the user's download confirmation for the actual request and access dates.

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
