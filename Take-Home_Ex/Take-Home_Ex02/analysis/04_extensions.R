# Additional uncertainty checks and explanatory cases. Raw data are read only.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(Kendall)})
d <- readRDS('outputs/prepared.rds')
s <- readRDS('outputs/statistics.rds')
g <- d$townships
n <- nrow(g); nt <- length(d$months)
# Reuse the authoritative functions without executing the statistical pipeline.
for(expr in parse('analysis/02_statistics.R')) {
  if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
     as.character(expr[[2]]) %in% c('classify_pattern','local_results')) eval(expr)
}
source('analysis/fast_gi.R')

sen_slope <- function(y) {
  ij <- combn(seq_along(y),2)
  median((y[ij[2,]]-y[ij[1,]])/(ij[2,]-ij[1,]))
}
mk_result <- function(y) {
  if(length(unique(y))<2L) return(c(tau=0,p=1))
  a <- Kendall::MannKendall(y)
  if(!is.finite(a$sl)) return(c(tau=0,p=1))
  c(tau=as.numeric(a$tau),p=as.numeric(a$sl))
}
block_test <- function(y,width,B=999L,seed=20261018L) {
  set.seed(seed)
  blocks <- split(seq_along(y),ceiling(seq_along(y)/width))
  observed <- mk_result(y)[['tau']]
  null <- replicate(B,mk_result(y[unlist(blocks[sample.int(length(blocks))])])[['tau']])
  (1+sum(abs(null)>=abs(observed)-1e-12))/(B+1)
}
space_time_weights <- function(nb,lag,n,nt) {
  out <- vector('list',n*nt)
  for(t in seq_len(nt)) for(i in seq_len(n)) {
    spatial <- sort(unique(c(i,nb[[i]][nb[[i]]>0])))
    out[[(t-1L)*n+i]] <- sort(as.integer(unlist(lapply(seq.int(max(1,t-lag),t),function(j)(j-1L)*n+spatial))))
  }
  class(out)<-'nb';attr(out,'region.id')<-as.character(seq_along(out))
  attr(out,'type')<-'space-time';attr(out,'sym')<-FALSE;attr(out,'self.included')<-TRUE
  nb2listw(out,style='W',zero.policy=TRUE)
}
summarise_cube <- function(z,p,variable,scenario) {
  zm<-matrix(z,nrow=n);pm<-matrix(p,nrow=n)
  bind_rows(lapply(seq_len(n),function(i){
    mk<-mk_result(zm[i,])
    cl<-classify_pattern(zm[i,],pm[i,],mk[['tau']],mk[['p']])
    tibble(code=g$TS_PCODE[i],variable=variable,scenario=scenario,
           tau=mk[['tau']],mk_p=mk[['p']],classification=cl,
           map_class=ifelse(mk[['p']]<.05 && cl!='No pattern detected',cl,'Not displayed'))
  }))
}

temporal <- list(); lag_results <- list(); bin_adjusted <- list(); slopes <- list()
for(v in c('events','fatalities')) {
  message('Trend robustness: ',v)
  gi<-s$gis %>% filter(variable==v)
  stopifnot(identical(gi$township_id,d$panel$township_id),identical(gi$month,d$panel$month))
  zm<-matrix(gi$gi_z,nrow=n);xm<-matrix(d$panel[[v]],nrow=n)
  baseline<-s$ehsa %>% filter(variable==v)
  temporal[[v]]<-bind_rows(lapply(seq_len(n),function(i){
    tibble(code=g$TS_PCODE[i],variable=v,block_3=baseline$block_p[i],
      block_6=block_test(zm[i,],6,seed=20262000L+i),
      block_12=block_test(zm[i,],12,seed=20263000L+i),
      lag1_acf=as.numeric(acf(zm[i,],plot=FALSE,lag.max=1)$acf[2]))
  }))
  slopes[[v]]<-tibble(code=g$TS_PCODE,variable=v,
    gi_slope=apply(zm,1,sen_slope),count_slope=apply(xm,1,sen_slope))
  bin_adjusted[[v]]<-summarise_cube(gi$gi_z,gi$q,v,'BH-adjusted space-time bins')
  for(lag in c(0L,2L)) {
    message('Alternative time neighbourhood: ',v,', lag ',lag)
    w<-space_time_weights(s$nb,lag,n,nt)
    result<-fast_space_time_gi(d$panel[[v]],w,nsim=999L,seed=20261018L)
    lag_results[[paste(v,lag)]]<-summarise_cube(result$z,result$p,v,paste0('Previous months: ',lag))
  }
}
temporal<-bind_rows(temporal);lag_results<-bind_rows(lag_results)
bin_adjusted<-bind_rows(bin_adjusted);slopes<-bind_rows(slopes)
checks<-s$ehsa %>% left_join(temporal,by=c('code','variable')) %>% left_join(slopes,by=c('code','variable')) %>%
  mutate(all_blocks=block_3<.05 & block_6<.05 & block_12<.05,
         trend_supported=mk_p<.05 & mk_q<.05 & all_blocks,
         evidence=case_when(map_class=='Not displayed'~'No displayed trend-pattern',
           trend_supported~'BH and all block checks',TRUE~'Sensitivity-dependent pattern'))
temporal_summary<-checks %>% group_by(variable) %>% summarise(nominal=sum(mk_p<.05),
  block_3=sum(mk_p<.05 & block_3<.05),block_6=sum(mk_p<.05 & block_6<.05),
  block_12=sum(mk_p<.05 & block_12<.05),BH_and_all_blocks=sum(trend_supported),
  displayed=sum(map_class!='Not displayed'),supported_displayed=sum(map_class!='Not displayed' & trend_supported),.groups='drop')

comparison<-bind_rows(lag_results,bin_adjusted) %>%
  left_join(s$ehsa %>% select(code,variable,baseline_class=classification,baseline_map=map_class),by=c('code','variable'))
scenario_summary<-comparison %>% group_by(variable,scenario) %>% summarise(
  class_agreement=mean(classification==baseline_class),
  displayed=sum(map_class!='Not displayed'),
  retained_same_pattern=sum(map_class==baseline_map & baseline_map!='Not displayed'),
  baseline_displayed=sum(baseline_map!='Not displayed'),.groups='drop')

# Extend the original spatial sensitivity checks to the second outcome.
xy<-st_coordinates(st_point_on_surface(g));kn6<-make.sym.nb(knn2nb(knearneigh(xy,k=6)))
lw<-nb2listw(s$nb,style='W',zero.policy=TRUE);star<-nb2listw(include.self(s$nb),style='W',zero.policy=TRUE)
precise<-d$assigned %>% filter(geo_precision==1) %>% group_by(township_id) %>% summarise(fatalities=sum(fatalities),.groups='drop')
near<-d$unassigned %>% filter(match_count==0,nearest_m<=1000) %>% group_by(nearest_id) %>% summarise(fatalities=sum(fatalities),.groups='drop')
sens_fatalities<-bind_rows(
 bind_cols(tibble(scenario='Six-nearest-neighbour graph',code=g$TS_PCODE),local_results(g$fatalities,nb2listw(kn6,style='W'),nb2listw(include.self(kn6),style='W'))),
 bind_cols(tibble(scenario='Precisely located events only',code=g$TS_PCODE),local_results(replace_na(precise$fatalities[match(g$TS_PCODE,precise$township_id)],0),lw,star)),
 bind_cols(tibble(scenario='Add records within 1 km of boundary',code=g$TS_PCODE),local_results(g$fatalities+replace_na(near$fatalities[match(g$TS_PCODE,near$nearest_id)],0),lw,star)))
spatial<-bind_rows(mutate(s$sensitivity,variable='events'),mutate(sens_fatalities,variable='fatalities')) %>%
 left_join(s$static %>% select(code,variable,baseline_hotspot=hotspot),by=c('code','variable'))
spatial_summary<-spatial %>% group_by(variable,scenario) %>% summarise(
 baseline_hot=sum(baseline_hotspot=='Hot spot'),alternative_hot=sum(hotspot=='Hot spot'),
 retained=sum(hotspot=='Hot spot' & baseline_hotspot=='Hot spot'),
 gained=sum(hotspot=='Hot spot' & baseline_hotspot!='Hot spot'),
 lost=sum(hotspot!='Hot spot' & baseline_hotspot=='Hot spot'),
 Jaccard=sum(hotspot=='Hot spot' & baseline_hotspot=='Hot spot')/sum(hotspot=='Hot spot' | baseline_hotspot=='Hot spot'),.groups='drop')

# Value-added question: where do frequency and fatality concentration disagree?
# Use separately BH-adjusted Gi* families. No ranking of individual risk is implied.
ev<-s$static %>% filter(variable=='events') %>% select(code,event_z=gi_z,event_q=gi_q)
fa<-s$static %>% filter(variable=='fatalities') %>% select(code,fatality_z=gi_z,fatality_q=gi_q)
concordance<-inner_join(ev,fa,by='code') %>%
 mutate(event_hot=event_z>0 & event_q<.05,fatality_hot=fatality_z>0 & fatality_q<.05,
  category=case_when(event_hot & fatality_hot~'Both outcomes',event_hot~'Events only',fatality_hot~'Fatalities only',TRUE~'Neither outcome')) %>%
 left_join(st_drop_geometry(g) %>% select(code=TS_PCODE,township=TS,region=ST,events,fatalities),by='code')
concordance_summary<-concordance %>% group_by(category) %>% summarise(townships=n(),events=sum(events),fatalities=sum(fatalities),.groups='drop')
joint_stability<-spatial %>% mutate(hot=gi_z>0 & gi_q<.05) %>% select(code,variable,scenario,hot) %>%
 pivot_wider(names_from=variable,values_from=hot) %>%
 left_join(concordance %>% select(code,baseline_both=category),by='code') %>% group_by(scenario) %>%
 summarise(both=sum(events & fatalities),retained_baseline_both=sum(events & fatalities & baseline_both=='Both outcomes'),.groups='drop')

# Conditional spatial random-labelling test of incident severity.
# Hold incident locations, months and types fixed. Shuffle fatality marks only
# within month x event-type strata, preserving each stratum's fatality distribution.
# The self-including queen graph aggregates marks into local neighbourhoods.
severity_test <- function(records,scenario,B=999L,seed=20264000L) {
  idx<-match(records$township_id,g$TS_PCODE)
  aggregate_town<-function(x) as.numeric(rowsum(c(x,rep(0,n)),c(idx,seq_len(n)),reorder=TRUE))
  strata<-split(seq_len(nrow(records)),interaction(records$month,records$event_type,drop=TRUE))
  marks<-records$fatalities
  expected_marks<-marks
  for(ids in strata) expected_marks[ids]<-mean(marks[ids])
  adjacency<-nb2mat(include.self(s$nb),style='B',zero.policy=TRUE)
  observed<-as.numeric(adjacency %*% aggregate_town(marks))
  expected<-as.numeric(adjacency %*% aggregate_town(expected_marks))
  upper<-integer(n)
  set.seed(seed)
  for(b in seq_len(B)) {
    shuffled<-marks
    for(ids in strata) shuffled[ids]<-marks[ids[sample.int(length(ids))]]
    simulated<-as.numeric(adjacency %*% aggregate_town(shuffled))
    upper<-upper+as.integer(simulated>=observed)
  }
  p<-(upper+1)/(B+1)
  tibble(code=g$TS_PCODE,scenario=scenario,neighbourhood_events=as.numeric(adjacency %*% aggregate_town(rep(1,length(idx)))),
    observed=observed,expected=expected,ratio=ifelse(expected>0,observed/expected,NA_real_),p=p,q=p.adjust(p,'BH'))
}
message('Conditional severity permutation, all assigned events')
severity<-severity_test(d$assigned,'All assigned events')
message('Conditional severity permutation, precision code 1')
severity_precise<-severity_test(filter(d$assigned,geo_precision==1),'Precision code 1')
severity<-severity %>% left_join(severity_precise %>% select(code,precise_ratio=ratio,precise_q=q),by='code') %>%
 left_join(concordance %>% select(code,township,region,category),by='code') %>%
 mutate(excess=q<.05 & observed>expected,precise_excess=precise_q<.05 & precise_ratio>1)
write.csv(severity,'outputs/conditional_severity.csv',row.names=FALSE)

# Contrastive, reproducible selection. These four examples are not a random sample.
case_pool<-concordance %>% left_join(checks %>% filter(variable=='events') %>%
 select(code,classification,map_class,trend_supported),by='code')
used<-character();cases<-list()
choose_case<-function(pool,role,rankvar='events') {
  a<-pool %>% filter(!code %in% used) %>% arrange(desc(.data[[rankvar]]),code) %>% slice_head(n=1)
  stopifnot(nrow(a)==1)
  used<<-c(used,a$code);a$selection_reason<-role;a
}
cases[[1]]<-choose_case(filter(case_pool,map_class=='Oscillating hot spot'),'Largest event total among displayed oscillating event hot spots')
cases[[2]]<-choose_case(filter(case_pool,map_class=='Sporadic hot spot'),'Largest event total among displayed sporadic event hot spots')
cases[[3]]<-choose_case(filter(case_pool,map_class=='Not displayed'),'Largest event total without a displayed event trend-pattern')
cases[[4]]<-choose_case(filter(case_pool,category=='Fatalities only'),'Largest fatality total among fatality-only BH hot spots','fatalities')
cases<-bind_rows(cases) %>% mutate(case_id=paste0('Case ',row_number()))
case_stats<-checks %>% filter(code %in% cases$code) %>%
 left_join(cases %>% select(code,case_id,township,region,selection_reason),by='code')
case_trends<-d$panel %>% filter(township_id %in% cases$code) %>%
 pivot_longer(c(events,fatalities),names_to='variable',values_to='value') %>%
 left_join(s$gis,by=c('township_id','month','variable')) %>%
 left_join(cases %>% select(township_id=code,case_id,township,region),by='township_id')
matched_period<-case_trends %>% filter(format(month,'%Y') %in% c('2024','2025'),as.integer(format(month,'%m'))<=9) %>%
 group_by(township_id,variable,year=format(month,'%Y')) %>% summarise(total=sum(value),.groups='drop') %>%
 pivot_wider(names_from=year,values_from=total,names_prefix='Jan_Sep_') %>%
 mutate(absolute_change=Jan_Sep_2025-Jan_Sep_2024,
        percentage_change=ifelse(Jan_Sep_2024>0,100*absolute_change/Jan_Sep_2024,NA_real_))
case_stats<-case_stats %>% left_join(matched_period,by=c('code'='township_id','variable'))
saveRDS(list(checks=checks,temporal_summary=temporal_summary,comparison=comparison,scenario_summary=scenario_summary,
 spatial=spatial,spatial_summary=spatial_summary,concordance=concordance,concordance_summary=concordance_summary,
 joint_stability=joint_stability,severity=severity,cases=cases,case_stats=case_stats,case_trends=case_trends),'outputs/extensions.rds')
write.csv(checks,'outputs/trend_robustness.csv',row.names=FALSE)
write.csv(comparison,'outputs/ehsa_sensitivity.csv',row.names=FALSE)
write.csv(spatial_summary,'outputs/hotspot_stability.csv',row.names=FALSE)
write.csv(concordance,'outputs/outcome_concordance.csv',row.names=FALSE)
write.csv(case_stats,'outputs/township_cases.csv',row.names=FALSE)
print(temporal_summary);print(scenario_summary);print(spatial_summary);print(concordance_summary)
print(cases %>% select(case_id,township,region,events,fatalities,classification,category))
