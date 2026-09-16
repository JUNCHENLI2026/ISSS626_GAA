# Learning demonstration: fixed, illustrative bandwidths, not an optimal selection.
# Run with the Ketapang data directory as the first argument.
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else ".."
out <- file.path(root,"learning")
targets <- file.path(out,c("kde_bandwidth_demo.png","kde_confidence_demo.png","kde_checks.txt"))
if(any(file.exists(targets)) && !"--replace-generated" %in% args) stop("Learning outputs exist; review before overwriting.")
inputs <- file.path(root,c("ketapang_boundary.gpkg","ketapang_modis_2026_conf30.csv","ketapang_modis_2026_conf80.csv"))
hash_before <- tools::md5sum(inputs)
boundary <- st_transform(st_read(inputs[1],layer="study_window",quiet=TRUE),32749)
stopifnot(all(st_is_valid(boundary)),!any(st_is_empty(boundary)))
window_m <- as.owin(boundary)
stopifnot(abs(area.owin(window_m)-sum(as.numeric(st_area(boundary))))/area.owin(window_m)<1e-8)
window_km <- rescale(window_m,1000,"km")
make_pattern <- function(file) {
  d <- read.csv(file,colClasses=c(acq_time="character"))
  p <- st_transform(st_as_sf(d,coords=c("longitude","latitude"),crs=4326),32749)
  xy <- st_coordinates(p)/1000
  stopifnot(all(inside.owin(xy[,1],xy[,2],window_km)))
  x <- ppp(xy[,1],xy[,2],window=window_km,checkdup=TRUE)
  stopifnot(npoints(x)==nrow(d))
  x
}
main <- make_pattern(inputs[2]); high <- make_pattern(inputs[3])
stopifnot(npoints(main)==4472,npoints(high)==1726)
sigmas <- c(5,10,20)
# Distances and bandwidths are in km; each point contributes weight 1.
# 0.5 km pixels control numerical rendering, not the satellite's resolution.
kde <- function(x,sigma) density(x,sigma=sigma,kernel="gaussian",edge=TRUE,diggle=TRUE,eps=0.5)
surfaces <- lapply(sigmas,function(s) kde(main,s))
high_surface <- kde(high,10)
checks <- data.frame(subset=c(rep("confidence_ge30",3),"confidence_ge80"),
  sigma_km=c(sigmas,10),records=c(rep(npoints(main),3),npoints(high)),
  integral=vapply(c(surfaces,list(high_surface)),integral.im,numeric(1)))
checks$relative_mass_error <- abs(checks$integral-checks$records)/checks$records
stopifnot(all(checks$relative_mass_error<0.05))
for(z in c(surfaces,list(high_surface))) {
  stopifnot(all(is.finite(z$v[!is.na(z$v)])),min(z$v,na.rm=TRUE)>-1e-8)
}
to_grid <- function(z,label) {
  d <- as.data.frame(z)
  names(d)[1:3] <- c("x","y","intensity")
  d <- d[!is.na(d$intensity),]
  # FFT round-off can produce negligible negatives; truncate only for display.
  d$intensity <- pmax(d$intensity,0)
  d$pixel_width <- z$xstep
  d$pixel_height <- z$ystep
  d$panel <- label
  d
}
bw_labels <- paste0("sigma = ",sigmas," km")
bw_grid <- do.call(rbind,Map(to_grid,surfaces,bw_labels))
bw_grid$panel <- factor(bw_grid$panel,levels=bw_labels)
conf_labels <- c("Confidence >=30 | n = 4,472","Confidence >=80 | n = 1,726")
conf_grid <- rbind(to_grid(surfaces[[2]],conf_labels[1]),to_grid(high_surface,conf_labels[2]))
conf_grid$panel <- factor(conf_grid$panel,levels=conf_labels)
coords <- st_coordinates(boundary)
outline <- data.frame(x=coords[,"X"]/1000,y=coords[,"Y"]/1000,
  ring=apply(coords[,grep("^L",colnames(coords)),drop=FALSE],1,paste,collapse="_"))
draw <- function(d,title,subtitle,caption) {
  ggplot(d,aes(x,y,fill=intensity))+
    geom_tile(aes(width=pixel_width,height=pixel_height))+
    geom_path(data=outline,aes(x,y,group=ring),inherit.aes=FALSE,colour="#64748B",linewidth=0.25)+
    facet_wrap(~panel,nrow=1)+coord_equal(expand=FALSE)+
    scale_fill_viridis_c(option="magma",limits=c(0,max(d$intensity)),
      name="Detections / km²\n(over Jan-Aug)",na.value="white")+
    labs(title=title,subtitle=subtitle,x="UTM zone 49S easting (km)",y="UTM zone 49S northing (km)",caption=caption)+
    theme_minimal(base_size=12)+theme(
      plot.title=element_text(face="bold",size=17),
      plot.subtitle=element_text(size=11,colour="#475569"),
      strip.text=element_text(face="bold",size=11),
      panel.grid=element_blank(),panel.background=element_rect(fill="white",colour=NA),
      plot.background=element_rect(fill="white",colour=NA),
      plot.caption=element_text(hjust=0,size=9,colour="#475569"),
      plot.caption.position="plot",plot.margin=margin(16,16,16,16),legend.position="right")
}
p1 <- draw(bw_grid,"Learning example: effect of KDE bandwidth",
  "Ketapang MODIS | Confidence >=30 | 1 January-31 August 2026 | Common colour scale",
  paste("Illustrative bandwidths, not an optimised or final choice. Gaussian kernel; Jones-Diggle edge correction; 0.5 km computational grid.",
        "Intensity refers to recorded detections, not independent fires or risk. Sources: NASA FIRMS and the supplied Indonesia Geospatial boundary.",sep="\n"))
p2 <- draw(conf_grid,"Learning example: confidence thresholds at the same bandwidth",
  "Ketapang MODIS | sigma = 10 km in both panels | 1 January-31 August 2026 | Common colour scale",
  paste("Absolute intensity is not normalised by record count; fewer retained records can reduce intensity. The >=80 subset is nested within >=30.",
        "These are descriptive maps, not tests of clustering. Sources: NASA FIRMS and the supplied Indonesia Geospatial boundary.",sep="\n"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(targets[1],p1,width=15,height=9,dpi=160,bg="white")
ggsave(targets[2],p2,width=12,height=9,dpi=160,bg="white")
stopifnot(identical(hash_before,tools::md5sum(inputs)))
capture.output({
  cat("Working projection: EPSG:32749; computational distance unit: km\n")
  cat("Derived study-window area (km2):",area.owin(window_km),"\n")
  cat("All records inside polygon window; no dropped points.\n")
  print(checks,row.names=FALSE)
  cat("Original inputs unchanged; MD5:\n");print(hash_before)
  cat("\nRuntime versions:\n");print(sessionInfo())
},file=targets[3])
print(checks,row.names=FALSE)
cat("Study-window area (km2):",area.owin(window_km),"\n")
cat("Outputs:\n",paste(targets,collapse="\n"),"\n")
