# Hunan data for Hands-on Exercise 4

The exercise uses the Hunan county boundaries and `Hunan_2012.csv` specified in Chapter 8. The statistical reference year is **2012**, not the download year. The boundary's effective reference date is not independently established by its filename.

Source: [the instructor's current course repository](https://github.com/tskam/ISSS626-AY2026-27Aug/tree/966eea6df498b38f10bc209d07a922f4d87b95d9/lesson/Lesson04/data).

Pinned commit: `966eea6df498b38f10bc209d07a922f4d87b95d9`.

Retrieved: 2026-09-16. The six downloaded files are unmodified. `download-manifest.txt` records their source URLs and MD5 checksums for reproducibility.

- `aspatial/Hunan_2012.csv`: 88 counties and 29 indicator/identifier columns.
- `geospatial/Hunan.shp`: 88 polygon features in WGS 84, with seven attribute fields.
- `geospatial/Hunan.dbf`, `.shx`, `.prj`, `.qpj`: shapefile companion files. Keep them together with the `.shp` file.

From the exercise directory, `Rscript data/download-data.R` downloads the files into an empty data location. Existing destination files are protected from overwriting. Normal page rendering uses the existing local files and does not download them again.

The Quarto page validates the county-key join, missing GDPPC values, geometry validity, neighbour counts, weight row sums and lag identities. It also checks that the raw files remain unchanged. No source GDPPC values, county names or boundary coordinates are rewritten.
