# Stage-specific content cache: figure edits do not repeat expensive inference.
required_packages <- c('sf','dplyr','tidyr','readr','lubridate','spdep','Kendall','ggplot2','scales','knitr','Rcpp')
missing <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages before rendering: ',paste(missing,collapse=', '))
versions <- vapply(required_packages,function(p)as.character(packageVersion(p)),character(1))
prepare_signature <- list(hashes=tools::md5sum(c(list.files('data',recursive=TRUE,full.names=TRUE),'analysis/01_prepare.R')),versions=versions)
statistics_signature <- list(preparation=prepare_signature,script=tools::md5sum(c('analysis/02_statistics.R','analysis/fast_gi.R')))
figures_signature <- list(statistics=statistics_signature,script=tools::md5sum('analysis/03_figures.R'))
extensions_signature <- list(statistics=statistics_signature,script=tools::md5sum('analysis/04_extensions.R'))
extension_figures_signature <- list(extensions=extensions_signature,script=tools::md5sum('analysis/05_extension_figures.R'))
evidence_signature <- list(extensions=extensions_signature,script=tools::md5sum('analysis/06_evidence_audit.R'))
evidence_figures_signature <- list(evidence=evidence_signature,script=tools::md5sum('analysis/07_evidence_figures.R'))
validation_signature <- list(evidence=evidence_signature,script=tools::md5sum('analysis/08_severity_validation.R'))
global_signature <- list(statistics=statistics_signature,script=tools::md5sum('analysis/09_global_moran.R'))
cache_file <- 'outputs/build_signature.rds'
cache <- if(file.exists(cache_file))readRDS(cache_file) else list()
force <- Sys.getenv('TH02_REBUILD')=='1'
prepare_changed <- force || !file.exists('outputs/prepared.rds') || !identical(prepare_signature,cache$prepare)
if(prepare_changed){
  source('analysis/01_prepare.R',local=TRUE)
  cache$prepare<-prepare_signature
  cache$statistics<-NULL;cache$figures<-NULL
  saveRDS(cache,cache_file)
}
statistics_changed <- prepare_changed || !file.exists('outputs/statistics.rds') || !identical(statistics_signature,cache$statistics)
if(statistics_changed){
  source('analysis/02_statistics.R',local=TRUE)
  cache$statistics<-statistics_signature;cache$figures<-NULL
  saveRDS(cache,cache_file)
}
figure_files <- paste0('outputs/figures/',c('burden','monthly','events_lisa','events_hotspot','fatalities_lisa','fatalities_hotspot','events_ehsa','fatalities_ehsa','gi_heatmap','multiplicity','top_townships','events_lisa_slide','fatalities_hotspot_slide','events_ehsa_slide','fatalities_ehsa_slide','multiplicity_slide'),'.png')
if(statistics_changed || !all(file.exists(c('outputs/figure_data.rds',figure_files))) || !identical(figures_signature,cache$figures)){
  source('analysis/03_figures.R',local=TRUE)
  cache$figures<-figures_signature
  saveRDS(cache,cache_file)
}
extensions_changed <- statistics_changed || !file.exists('outputs/extensions.rds') || !identical(extensions_signature,cache$extensions)
if(extensions_changed){
  source('analysis/04_extensions.R',local=TRUE)
  cache$extensions<-extensions_signature
  cache$extension_figures<-NULL
  saveRDS(cache,cache_file)
}
extension_files <- paste0('outputs/figures/',c('trend_robustness','outcome_concordance','conditional_severity',
 'conditional_severity_slide','case_locator',paste0('case_',1:4)),'.png')
if(extensions_changed || !all(file.exists(extension_files)) || !identical(extension_figures_signature,cache$extension_figures)){
  source('analysis/05_extension_figures.R',local=TRUE)
  cache$extension_figures<-extension_figures_signature
  saveRDS(cache,cache_file)
}
evidence_changed <- extensions_changed || !file.exists('outputs/evidence_audit.rds') || !identical(evidence_signature,cache$evidence)
if(evidence_changed){
  source('analysis/06_evidence_audit.R',local=TRUE)
  cache$evidence<-evidence_signature
  saveRDS(cache,cache_file)
}
evidence_files<-paste0('outputs/figures/',c('national_share_events','national_share_fatalities',
 'robust_severity_types','robust_severity_periods','robust_severity_neighbours'),'.png')
if(evidence_changed || !all(file.exists(evidence_files)) || !identical(evidence_figures_signature,cache$evidence_figures)){
  source('analysis/07_evidence_figures.R',local=TRUE)
  cache$evidence_figures<-evidence_figures_signature
  saveRDS(cache,cache_file)
}
if(evidence_changed || !file.exists('outputs/severity_validation.rds') || !identical(validation_signature,cache$validation)){
  source('analysis/08_severity_validation.R',local=TRUE)
  cache$validation<-validation_signature
  saveRDS(cache,cache_file)
}
if(statistics_changed || !all(file.exists(c('outputs/global_moran.rds','outputs/global_moran_summary.csv','outputs/global_moran_sensitivity.csv','outputs/figures/global_moran_scatter.png'))) || !identical(global_signature,cache$global)){
  source('analysis/09_global_moran.R',local=new.env())
  cache$global<-global_signature
  saveRDS(cache,cache_file)
}
