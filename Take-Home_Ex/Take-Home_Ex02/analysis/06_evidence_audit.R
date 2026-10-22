# Diagnose naming flags, separate national change, and profile stable severity cases.
# Only aggregated outputs are saved. Raw incident records are never modified.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep)})
d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds');e<-readRDS('outputs/extensions.rds')
g<-d$townships;n<-nrow(g)
for(expr in parse('analysis/02_statistics.R')) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
    as.character(expr[[2]])=='local_results')eval(expr)
}
normalise_name<-function(x)gsub('[^a-z0-9]','',tolower(x))
lookup<-st_drop_geometry(g) %>% transmute(code=TS_PCODE,name=TS,region=ST,
 key=normalise_name(TS),base_key=normalise_name(sub(' *\\(.*$','',TS)),region_key=normalise_name(ST))
a<-d$assigned %>% select(-region) %>% left_join(lookup %>% select(code,name,region),by=c('township_id'='code')) %>%
 mutate(name_missing=is.na(admin3)|trimws(admin3)=='',
 exact_match=!name_missing & tolower(trimws(admin3))==tolower(trimws(name)),
 name_flag=!exact_match,key=normalise_name(admin3),region_key=normalise_name(admin1))
candidate_code<-function(key,region_key) {
 if(is.na(key)||key=='')return(NA_character_)
 hit<-which(lookup$key==key & lookup$region_key==region_key)
 if(length(hit)==0L)hit<-which(lookup$base_key==key & lookup$region_key==region_key)
 if(length(hit)==1L)lookup$code[hit] else NA_character_
}
a$candidate<-mapply(candidate_code,a$key,a$region_key,USE.NAMES=FALSE)
a<-a %>% mutate(audit_class=case_when(
 name_missing~'Missing township name',exact_match~'Exact name match',
 key==normalise_name(name)~'Formatting-only difference',
 candidate==township_id~'Boundary primary-name match',
 !is.na(candidate) & candidate!=township_id~'Conflicting unique township label',
 TRUE~'Unresolved name'))
audit_summary<-a %>% group_by(audit_class) %>% summarise(records=n(),fatalities=sum(fatalities),.groups='drop')
name_pairs<-a %>% filter(name_flag) %>% count(admin1,admin3,region,name,audit_class,name='records',sort=TRUE)
# Geometry evidence for conflicting labels, without exposing incident coordinates.
conflicting<-a %>% filter(audit_class=='Conflicting unique township label')
if(nrow(conflicting)) {
 pts<-st_as_sf(conflicting,coords=c('longitude','latitude'),crs=4326) %>% st_transform(st_crs(g))
 conflict_distance<-as.numeric(st_distance(pts,g[match(conflicting$candidate,g$TS_PCODE),],by_element=TRUE))/1000
 conflicts<-conflicting %>% mutate(label_distance_km=conflict_distance,
   labelled_township=lookup$name[match(candidate,lookup$code)],
   adjacent=mapply(function(from,to)match(to,g$TS_PCODE)%in%s$nb[[match(from,g$TS_PCODE)]],township_id,candidate)) %>%
   group_by(admin3,name,region,labelled_township) %>% summarise(records=n(),fatalities=sum(fatalities),
    minimum_km=min(label_distance_km),maximum_km=max(label_distance_km),
    all_adjacent=all(adjacent),precise_records=sum(geo_precision==1),.groups='drop')
}else conflicts<-tibble()
unresolved<-a$audit_class %in% c('Missing township name','Conflicting unique township label','Unresolved name')
audit_sets<-list('Exclude unresolved or missing labels'=a[!unresolved,],
 'Exclude every initial name flag'=a[!a$name_flag,],
 'Use unique conflicting labels instead of coordinates'=a %>% mutate(township_id=if_else(
   audit_class=='Conflicting unique township label',candidate,township_id)))
lw<-nb2listw(s$nb,style='W',zero.policy=TRUE)
star<-nb2listw(include.self(s$nb),style='W',zero.policy=TRUE)
audit_results<-list();audit_counts<-list()
for(scenario in names(audit_sets)) {
 z<-audit_sets[[scenario]]
 totals<-z %>% group_by(township_id) %>% summarise(events=n(),fatalities=sum(fatalities),.groups='drop')
 audit_counts[[scenario]]<-tibble(scenario=scenario,retained=nrow(z),excluded=nrow(a)-nrow(z),
  relabelled=if(scenario=='Use unique conflicting labels instead of coordinates')nrow(conflicting) else 0L)
 for(v in c('events','fatalities')){
  x<-replace_na(totals[[v]][match(g$TS_PCODE,totals$township_id)],0)
  out<-local_results(x,lw,star)
  base<-s$static %>% filter(variable==v)
  audit_results[[paste(scenario,v)]]<-bind_cols(tibble(code=g$TS_PCODE,scenario=scenario,variable=v),out) %>%
   mutate(baseline_hot=base$gi_z>0 & base$gi_q<.05,hot=gi_z>0 & gi_q<.05,
    baseline_lisa=ifelse(base$moran_q<.05,base$lisa,'Not significant'),
    adjusted_lisa=ifelse(moran_q<.05,lisa,'Not significant'))
 }
}
audit_results<-bind_rows(audit_results)
audit_stability<-audit_results %>% group_by(scenario,variable) %>% summarise(
 baseline_n=sum(baseline_hot),alternative_hot=sum(hot),retained=sum(hot & baseline_hot),
 gained=sum(hot & !baseline_hot),lost=sum(!hot & baseline_hot),
 Jaccard=sum(hot & baseline_hot)/sum(hot|baseline_hot),
 lisa_label_changes=sum(adjusted_lisa!=baseline_lisa),.groups='drop') %>% rename(baseline_hot=baseline_n)

# Descriptive national comparison: constant-share benchmark, not a causal model.
national<-d$panel %>% group_by(month) %>% summarise(across(c(events,fatalities),sum),.groups='drop') %>%
 pivot_longer(c(events,fatalities),names_to='variable',values_to='national')
relative_series<-d$panel %>% pivot_longer(c(events,fatalities),names_to='variable',values_to='local') %>%
 left_join(national,by=c('month','variable')) %>% mutate(share=ifelse(national>0,local/national,NA_real_))
national_windows<-national %>% filter(as.integer(format(month,'%m'))<=9,format(month,'%Y')%in%c('2024','2025')) %>%
 group_by(variable,year=format(month,'%Y')) %>% summarise(total=sum(national),.groups='drop') %>%
 pivot_wider(names_from=year,values_from=total,names_prefix='national_') %>%
 mutate(national_growth_pct=100*(national_2025/national_2024-1))
relative_comparison<-relative_series %>% filter(as.integer(format(month,'%m'))<=9,format(month,'%Y')%in%c('2024','2025')) %>%
 group_by(township_id,variable,year=format(month,'%Y')) %>% summarise(total=sum(local),.groups='drop') %>%
 pivot_wider(names_from=year,values_from=total,names_prefix='local_') %>% left_join(national_windows,by='variable') %>%
 mutate(share_2024=100*local_2024/national_2024,share_2025=100*local_2025/national_2025,
 share_change_pp=share_2025-share_2024,expected_2025=local_2024/national_2024*national_2025,
 local_growth_pct=ifelse(local_2024>0,100*(local_2025/local_2024-1),NA_real_),
 relative_growth=ifelse(local_2024>0,(local_2025/local_2024)/(national_2025/national_2024),NA_real_))
case_comparison<-relative_comparison %>% inner_join(e$cases %>% select(code,case_id,township,region),by=c('township_id'='code')) %>% arrange(case_id,variable)
case_shares<-relative_series %>% inner_join(e$cases %>% select(code,case_id,township,region),by=c('township_id'='code')) %>%
 mutate(quarter=as.Date(paste0(format(month,'%Y'),'-',sprintf('%02d',3*((as.integer(format(month,'%m'))-1)%/%3)+1),'-01'))) %>%
 group_by(case_id,township,region,quarter,variable) %>% summarise(local=sum(local),national=sum(national),.groups='drop') %>%
 mutate(share_pct=100*local/national)

# Select the largest lower-of-two O/E ratios from the already supported set.
# A Jaccard ceiling prevents selecting nearly identical neighbourhoods.
pool<-e$severity %>% filter(excess,precise_excess) %>% mutate(minimum_ratio=pmin(ratio,precise_ratio)) %>%
 arrange(desc(minimum_ratio),code)
neighbours<-function(code)sort(unique(c(match(code,g$TS_PCODE),s$nb[[match(code,g$TS_PCODE)]])))
selected<-character()
for(code in pool$code){
 ids<-neighbours(code)
 overlaps<-vapply(selected,function(other){b<-neighbours(other);length(intersect(ids,b))/length(union(ids,b))},numeric(1))
 if(!length(overlaps)||all(overlaps<=.5))selected<-c(selected,code)
 if(length(selected)==3L)break
}
stopifnot(length(selected)==3L)
robust_cases<-pool %>% filter(code%in%selected) %>% arrange(match(code,selected)) %>% mutate(case_label=paste0('S',row_number()))
marked<-d$assigned %>% group_by(month,event_type) %>% mutate(expected_mark=mean(fatalities)) %>% ungroup()
decomposition<-list();annual<-list();concentration<-list();overlap<-list()
for(code in selected){
 z<-marked %>% filter(township_id%in%g$TS_PCODE[neighbours(code)])
 label<-robust_cases$township[match(code,robust_cases$code)]
 decomposition[[code]]<-z %>% group_by(event_type) %>% summarise(events=n(),observed=sum(fatalities),expected=sum(expected_mark),.groups='drop') %>%
  mutate(code=code,township=label,excess=observed-expected)
 annual[[code]]<-z %>% group_by(year=as.integer(format(month,'%Y'))) %>% summarise(events=n(),observed=sum(fatalities),expected=sum(expected_mark),.groups='drop') %>%
  mutate(code=code,township=label,excess=observed-expected)
 concentration[[code]]<-tibble(code=code,neighbour_townships=length(neighbours(code)),
  largest_incident_pct=100*max(z$fatalities)/sum(z$fatalities),
  largest_month_pct=100*max(tapply(z$fatalities,z$month,sum))/sum(z$fatalities),
  precise_incident_pct=100*mean(z$geo_precision==1),
  name_unresolved_records=sum(z$event_id_cnty%in%a$event_id_cnty[unresolved]))
}
robust_cases<-robust_cases %>% left_join(bind_rows(concentration),by='code')
for(i in seq_len(length(selected)-1))for(j in seq.int(i+1,length(selected))){
 ids<-neighbours(selected[i]);other<-neighbours(selected[j])
 overlap[[paste(i,j)]]<-tibble(first=robust_cases$township[i],second=robust_cases$township[j],
  shared_townships=length(intersect(ids,other)),Jaccard=length(intersect(ids,other))/length(union(ids,other)))
}
decomposition<-bind_rows(decomposition);annual<-bind_rows(annual)
drivers<-decomposition %>% group_by(code) %>% slice_max(excess,n=1,with_ties=FALSE) %>% ungroup() %>%
 left_join(robust_cases %>% select(code,total_excess=observed,expected_total=expected),by='code') %>%
 mutate(net_excess_share_pct=100*excess/(total_excess-expected_total))
peak_period<-annual %>% group_by(code) %>% slice_max(excess,n=1,with_ties=FALSE) %>% ungroup()
raw_hashes<-tools::md5sum(names(d$quality$hashes))
stopifnot(identical(unname(raw_hashes),unname(d$quality$hashes)))
result<-list(audit_summary=audit_summary,name_pairs=name_pairs,conflicts=conflicts,
 audit_counts=bind_rows(audit_counts),audit_results=audit_results,audit_stability=audit_stability,
 national_windows=national_windows,relative_comparison=relative_comparison,case_comparison=case_comparison,
 case_shares=case_shares,robust_cases=robust_cases,decomposition=decomposition,annual=annual,
 drivers=drivers,peak_period=peak_period,overlap=bind_rows(overlap))
saveRDS(result,'outputs/evidence_audit.rds')
for(item in c('audit_summary','name_pairs','conflicts','audit_stability','case_comparison','robust_cases','decomposition','annual'))
 write.csv(result[[item]],paste0('outputs/',item,'.csv'),row.names=FALSE,na='')
print(audit_summary);print(conflicts);print(audit_stability);print(national_windows)
print(case_comparison %>% select(township,variable,local_2024,local_2025,share_2024,share_2025,relative_growth))
print(robust_cases %>% select(case_label,township,region,ratio,precise_ratio,category,largest_incident_pct))
print(drivers);print(peak_period);print(result$overlap)
