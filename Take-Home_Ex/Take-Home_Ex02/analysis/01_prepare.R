# Prepare a complete township-month panel without changing the raw files.
suppressPackageStartupMessages({library(sf);library(dplyr);library(tidyr);library(readr);library(lubridate);library(spdep)})
dir.create('outputs',showWarnings=FALSE)
set.seed(20261018)
raw_paths <- c('data/Aspatial/ACLED_Data_Myanmar_Jan2021-Sep2025.csv', list.files('data/Geospatial',full.names=TRUE))
raw_hashes <- tools::md5sum(raw_paths)
raw <- read_csv(raw_paths[1],show_col_types=FALSE)
stopifnot(!anyDuplicated(raw$event_id_cnty))
raw$event_date <- as.Date(raw$event_date)
stopifnot(!anyNA(raw$event_date), all(raw$fatalities>=0), !anyNA(raw$fatalities))
event_types <- c('Battles','Explosions/Remote violence','Violence against civilians')
eligible <- raw %>% filter(country=='Myanmar',event_date>=as.Date('2021-01-01'),event_date<=as.Date('2025-09-30'),event_type %in% event_types)
stopifnot(all(is.finite(eligible$longitude)),all(is.finite(eligible$latitude)))
townships_raw <- st_read('data/Geospatial/mmr_polbnda_adm3_250k_mimu_1.shp',quiet=TRUE)
stopifnot(nrow(townships_raw)==330,all(as.character(townships_raw$PCode_V)=='9.4'),!anyDuplicated(townships_raw$TS_PCODE))
invalid_before <- sum(!st_is_valid(townships_raw))
crs_area <- '+proj=aea +lat_1=10 +lat_2=28 +lat_0=19 +lon_0=96 +datum=WGS84 +units=m +no_defs'
townships <- townships_raw %>% st_make_valid() %>% st_transform(crs_area) %>% arrange(TS_PCODE)
stopifnot(all(st_is_valid(townships)),!any(st_is_empty(townships)))
townships$area_km2 <- as.numeric(st_area(townships))/1e6
coords <- eligible %>% distinct(longitude,latitude)
pts <- st_as_sf(coords,coords=c('longitude','latitude'),crs=4326,remove=FALSE) %>% st_transform(st_crs(townships))
hits <- st_intersects(pts,townships)
coords$match_count <- lengths(hits)
coords$township_id <- vapply(hits,function(h) if(length(h)==1)townships$TS_PCODE[h] else NA_character_,character(1))
nearest <- st_nearest_feature(pts,townships)
coords$nearest_id <- townships$TS_PCODE[nearest]
coords$nearest_m <- as.numeric(st_distance(pts,townships[nearest,],by_element=TRUE))
events <- eligible %>% left_join(coords,by=c('longitude','latitude')) %>% mutate(month=as.Date(format(event_date,'%Y-%m-01')))
assigned <- events %>% filter(!is.na(township_id))
unassigned <- events %>% filter(is.na(township_id))
months <- seq(as.Date('2021-01-01'),as.Date('2025-09-01'),by='month')
panel <- expand_grid(township_id=townships$TS_PCODE,month=months) %>%
 left_join(assigned %>% group_by(township_id,month) %>% summarise(events=n(),fatalities=sum(fatalities),.groups='drop'),by=c('township_id','month')) %>%
 mutate(across(c(events,fatalities),~replace_na(.x,0))) %>% arrange(month,township_id)
stopifnot(nrow(panel)==330*57,sum(panel$events)==nrow(assigned),sum(panel$fatalities)==sum(assigned$fatalities),!anyDuplicated(panel[c('township_id','month')]))
totals <- panel %>% group_by(township_id) %>% summarise(events=sum(events),fatalities=sum(fatalities),active_months=sum(events>0),.groups='drop')
townships <- townships %>% left_join(totals,by=c('TS_PCODE'='township_id'))
audit <- tibble(stage=c('Input rows','Outside Myanmar','Outside requested dates','Other event types within Myanmar/time','Eligible armed-conflict records','Unique spatially assigned records','Unassigned or ambiguous records','Assigned reported fatalities','Invalid polygons repaired','Townships','Months','Township-month bins'),
 value=c(nrow(raw),sum(raw$country!='Myanmar'),sum(raw$event_date<min(months)|raw$event_date>as.Date('2025-09-30')),sum(raw$country=='Myanmar' & !(raw$event_type %in% event_types)),nrow(events),nrow(assigned),nrow(unassigned),sum(assigned$fatalities),invalid_before,nrow(townships),length(months),nrow(panel)))
precision <- assigned %>% count(geo_precision,time_precision,name='records')
name_check <- assigned %>% left_join(st_drop_geometry(townships) %>% select(TS_PCODE,TS),by=c('township_id'='TS_PCODE')) %>% summarise(missing_admin3=sum(is.na(admin3)),exact_name_disagreement=sum(!is.na(admin3)&tolower(trimws(admin3))!=tolower(trimws(TS))))
quality <- list(audit=audit,precision=precision,name_check=name_check,
 unmatched=unassigned %>% count(geo_precision,match_count,name='records'),
 near_boundary=sum(unassigned$match_count==0 & unassigned$nearest_m<=1000),
 event_types=raw %>% count(event_type,name='records'), hashes=raw_hashes)
saveRDS(list(townships=townships,panel=panel,assigned=assigned,unassigned=unassigned,months=months,quality=quality), 'outputs/prepared.rds')
write_csv(audit,'outputs/data_audit.csv')
print(audit,n=Inf);print(quality$unmatched);print(name_check)
cat('Unassigned within 1 km:',quality$near_boundary,'\n')
