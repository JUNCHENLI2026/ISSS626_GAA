# Take-Home Exercise 1: Ketapang MODIS detections

The technical report is `Take-Home_Ex01.qmd`. The executive summary is
`executive_summary.qmd` (10 content slides plus a cover). Both use Quarto and R.

## Requirements

- R with `sf`, `dplyr`, `ggplot2`, `spatstat.geom`, `spatstat.explore`, `knitr`
  and `rmarkdown`.
- Quarto. The report uses HTML and the summary uses revealjs.
- The derived boundary, detection CSV files and saved RDS results under
  `data/Ketapang/` for a fast render.

See `data/Ketapang/README.md` for source provenance, analytical stages, data
restrictions and the distinction between detection records and independent fires.

## Render the existing results

On Windows, the recommended command from this exercise directory is:

```powershell
.\render.ps1 -Verify
```

This runs the independent checks and renders the report, summary and submission
page. It temporarily uses a Windows-supported UTF-8 locale and a local temporary
directory, so paths containing Chinese characters work. It restores the calling
session's environment when finished. R and Quarto must already be installed.

On other platforms, or when the R locale already supports the project path, run:

```powershell
quarto render Take-Home_Ex01.qmd
quarto render executive_summary.qmd
```

When rendered inside the coursework website, output goes to
`../../_site/Take-Home_Ex/Take-Home_Ex01/`. Rendering uses saved calculations;
it does not rerun the complete permutation analysis.

If Windows reports an encoding error for the account's temporary directory,
set `TEMP` and `TMP` to an existing writable directory with an ASCII-only path
for that terminal session, then rerun Quarto. This does not require changing
permanent system settings.

## Independently check the saved calculations

```powershell
Rscript data/Ketapang/verify_revision.R data/Ketapang
Rscript data/Ketapang/verify_robustness.R data/Ketapang
```

These checks independently reconstruct pair counts, selected stored permutations,
Monte Carlo p-values, Holm adjustments and grid assignments. They verify
computation, not the scientific validity of exchangeability assumptions.
Rendering the report also checks source hashes and the 32 directional/grid
assignments using an independent rectangle-based calculation.

## Rebuild from raw inputs

```powershell
Rscript reproduce.R "C:/path/to/Take-Home_Ex01" "C:/path/to/new-build" "C:/path/to/quarto.exe"
```

Replace these example paths with the local paths. Use an unused, writable output
directory, preferably with an ASCII-only path on Windows. The runner refuses to
overwrite an existing directory. It sets an output-local temporary directory for
Quarto, writes stage logs and session information, and renders both documents.

Raw inputs must be present under `data/NASA FIRMS/` and
`data/Indonesia-geospasial/` as documented in the data README. The large national
boundary is intentionally excluded from Git. Reproducing extraction requires
obtaining those original source files. The retained Ketapang subset supports
the later analytical stages without that national boundary.

## Display changes

The report's main conditional-test section presents the final 192-test family.
Detailed specifications and historical comparison families remain in its appendix.
The summary uses the same final results. Its recurrence map comes from
`recurrence_map.R`; the optional `presentation` argument only adjusts label size
and the visible map extent. It does not alter input records or calculations.

The coursework navbar links to the report and summary. Publishing and GitHub
submission are separate steps from local rendering. The temporary whole-exercise Git exclusion has been removed so that source files
can be included in a submission commit. Do not add the excluded national boundary.

## Course-aligned presentation and portable checks

The report exposes short `library`, `st_read`, `st_transform`, `as.owin`,
`rescale`, `ppp`, `filter` and `select` steps, with links to the course chapters.
The displayed coordinates are checked against the saved analytical coordinates.
This does not replace the conditional test with a different reference model.

Checksum checks now resolve input filenames within the supplied data directory.
Saved hashes still validate the file contents after the project is copied to
another directory. No simulation result or source dataset is rewritten.

`index.qmd` groups the report, summary and GitHub repository links. The companion
submission source archive has a standalone `_quarto.yml`; it can render these
three pages without the rest of the coursework website.
