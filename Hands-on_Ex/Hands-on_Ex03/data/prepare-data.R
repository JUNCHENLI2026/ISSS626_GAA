library(sf)
library(readr)
library(dplyr)
# Run from the exercise directory; retain the original archive unchanged.
boundary <- st_read('data/rawdata/Kepulauan_Bangka_Belitung.shp', quiet=TRUE) |>
  st_zm() |> st_make_valid() |> st_transform(32748)
parts <- st_cast(st_union(boundary), 'POLYGON')
region <- parts[which.max(st_area(parts))]
st_write(st_sf(geometry=region), 'data/bangka_window.gpkg', delete_dsn=TRUE, quiet=TRUE)
raw <- read_csv('data/rawdata/fire_archive_M-C61_802038.csv', show_col_types=FALSE)
stopifnot(all(format(as.Date(raw$acq_date), '%Y') == '2023'))
points <- raw |> st_as_sf(coords=c('longitude','latitude'),crs=4326,remove=FALSE) |> st_transform(32748)
selected <- st_filter(points, region)
write_csv(st_drop_geometry(selected), 'data/forestfires.csv')
cat('Main-island subdistricts:', sum(lengths(st_intersects(boundary,region))>0), '\nFire detections:', nrow(selected), '\n')
print(range(selected$acq_date))
