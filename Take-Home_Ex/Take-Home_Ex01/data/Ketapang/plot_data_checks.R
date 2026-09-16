# Descriptive data checks only; no spatial significance tests or event merging.
# Run: Rscript plot_data_checks.R "path/to/data/Ketapang"
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
out <- file.path(root, "figures")
outputs <- file.path(out, c("ketapang_detection_locations.png", "ketapang_monthly_detections.png"))
if (any(file.exists(outputs))) stop("A plot already exists; review before replacing it.")
inputs <- file.path(root, c("ketapang_modis_2026_jan_aug.csv", "ketapang_modis_2026_conf30.csv",
                           "ketapang_modis_2026_conf80.csv", "monthly_confidence_comparison.csv", "ketapang_boundary.gpkg"))
before <- tools::md5sum(inputs)
read_detections <- function(f) read.csv(f, colClasses=c(acq_time="character"), stringsAsFactors=FALSE)
all_points <- read_detections(inputs[1])
main <- read_detections(inputs[2])
high <- read_detections(inputs[3])
counts <- read.csv(inputs[4], stringsAsFactors=FALSE)
boundary <- st_read(inputs[5], layer="study_window", quiet=TRUE)
stopifnot(nrow(all_points)==4738, nrow(main)==4472, nrow(high)==1726)
stopifnot(all(main$confidence>=30), all(high$confidence>=80))
for (i in seq_len(nrow(counts))) {
  month <- counts$month[i]
  stopifnot(sum(substr(all_points$acq_date,1,7)==month)==counts$all_detections[i],
            sum(substr(main$acq_date,1,7)==month)==counts$confidence_ge30[i],
            sum(substr(high$acq_date,1,7)==month)==counts$confidence_ge80[i])
}
main$subset <- "Confidence >=30 | n = 4,472"
high$subset <- "Confidence >=80 | n = 1,726"
map_data <- rbind(main, high)
map_data$subset <- factor(map_data$subset, levels=c(main$subset[1], high$subset[1]))
points_sf <- st_as_sf(map_data, coords=c("longitude","latitude"), crs=4326)
# Validate inclusion in the saved study window using the same working projection
# as extraction. Plot labels remain geographic coordinates for readability.
inside <- lengths(st_intersects(st_transform(points_sf,32749),st_transform(boundary,32749)))>0
stopifnot(all(inside))
theme_set(theme_minimal(base_size=12, base_family="sans"))
base_theme <- theme(
  plot.background=element_rect(fill="white",colour=NA),
  panel.background=element_rect(fill="#F4F7FA",colour=NA),
  plot.title=element_text(face="bold",size=18),
  plot.subtitle=element_text(colour="#475569",size=11),
  plot.caption=element_text(colour="#475569",size=9,hjust=0),
  plot.caption.position="plot",
  plot.margin=margin(16,20,16,16),
  strip.text=element_text(face="bold",size=12),
  panel.grid.major=element_line(colour="#DFE5EB",linewidth=0.25),
  panel.grid.minor=element_blank()
)
map <- ggplot() +
  geom_sf(data=boundary,fill="#E8EEE9",colour="#52665B",linewidth=0.35) +
  geom_sf(data=points_sf,colour="#B83E32",size=0.5,alpha=0.45,stroke=0) +
  facet_wrap(~subset,nrow=1) +
  coord_sf(crs=4326,datum=st_crs(4326),expand=TRUE) +
  labs(title="MODIS detections in Ketapang",
       subtitle="Observation window: 1 January-31 August 2026 | Same extent and point styling in both panels",
       x="Longitude",y="Latitude",
       caption=paste("Sources: NASA FIRMS MODIS; supplied Indonesia Geospatial village boundaries (dissolved to Ketapang).",
                     "Panels are nested subsets, not independent samples. Detections are not counts of unique fires; darker areas may reflect overlapping symbols.",sep="\n")) + base_theme
series <- c("All detections", "Confidence >=30", "Confidence >=80")
monthly <- do.call(rbind,lapply(seq_along(series),function(i) {
  data.frame(month=factor(month.abb[1:8],levels=month.abb[1:8]),
             subset=series[i], detections=counts[[c("all_detections","confidence_ge30","confidence_ge80")[i]]])
}))
monthly$subset <- factor(monthly$subset,levels=series)
dodge <- position_dodge(width=0.84)
bars <- ggplot(monthly,aes(month,detections,fill=subset)) +
  geom_col(position=dodge,width=0.78) +
  geom_text(aes(label=format(detections,big.mark=",",trim=TRUE)),position=dodge,vjust=-0.45,size=3.1) +
  scale_fill_manual(values=c("All detections"="#94A3B8","Confidence >=30"="#176B80","Confidence >=80"="#BE632B"),name=NULL) +
  scale_y_continuous(labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE),
                     limits=c(0,NA),expand=expansion(mult=c(0,0.12))) +
  labs(title="Monthly MODIS detection counts",
       subtitle="Ketapang | January-August 2026 | Counts compared under the confirmed confidence thresholds",
       x=NULL,y="Number of detections",
       caption=paste("Source: NASA FIRMS MODIS, spatially filtered to the supplied Ketapang boundary. January-April: standard; May-August: NRT.",
                     "Series are nested and must not be added. A zero means no records retained under that filter, not proof of no fires.",sep="\n")) +
  base_theme + theme(legend.position="top",legend.justification="left",panel.grid.major.x=element_blank())
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(outputs[1],map,width=14,height=8,dpi=160,bg="white")
ggsave(outputs[2],bars,width=14,height=7,dpi=160,bg="white")
stopifnot(identical(before,tools::md5sum(inputs)))
cat("All displayed points are inside the saved boundary. Monthly counts match all three CSV inputs.\n")
cat("Input hashes unchanged.\n")
cat(paste(outputs,collapse="\n"),"\n")
