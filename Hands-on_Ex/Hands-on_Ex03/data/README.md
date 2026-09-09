# Exercise 3 data provenance

Retrieved 2026-09-09 from https://github.com/xhwang0226/ISSS626/tree/main/Hands-on/Hands-on_Ex05/data . This is a third-party course data mirror, not a direct NASA download. The textbook attributes the archive to NASA FIRMS (MODIS) and boundaries to Indonesia Geospatial.

The original `fire_archive_M-C61_802038.csv` contains 2023 observations. All four original shapefile components are preserved under `rawdata/`. The mirrored province boundary has 393 features, including Belitung.

Run `source('data/prepare-data.R')` from the exercise directory to recreate the derived data. The largest connected polygon is Bangka main island; it intersects 297 original administrative units and matches the chapter's extent (105.108509 to 106.848839 E; -3.116593 to -1.501603 N). Output: `bangka_window.gpkg` (EPSG:32748) and `forestfires.csv` (744 detections, 2023-01-10 through 2023-12-18).

The textbook reports 741 observations. The reason for this difference is unverified. No artificial records or count-matching deletions are used. All plots and numerical results use the 744 local observations.
