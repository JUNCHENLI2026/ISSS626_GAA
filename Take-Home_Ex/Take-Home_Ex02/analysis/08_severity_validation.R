# Validate the previously selected cases without selecting replacement examples.
# Every randomisation and deletion scenario retains the full 330-test BH family.
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep)})
d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds')
h<-readRDS('outputs/evidence_audit.rds')
g<-d$townships;n<-nrow(g)
for(expr in parse('analysis/04_extensions.R')) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
    as.character(expr[[2]])=='severity_test')eval(expr)
}
seeds<-c(20265001L,20265002L,20265003L)
B<-1999L
samples<-list('All assigned events'=d$assigned,
 'Precision code 1'=filter(d$assigned,geo_precision==1))
repeated<-list();deleted<-list();deletion_audit<-list()
for(sample_name in names(samples)){
 records<-samples[[sample_name]]
 for(seed in seeds){
  message('Severity seed check: ',sample_name,', seed ',seed)
  out<-severity_test(records,sample_name,B=B,seed=seed)
  repeated[[paste(sample_name,seed)]]<-out %>% mutate(seed=seed,
   supported=q<.05 & observed>expected)
 }
 for(code in h$robust_cases$code){
  ids<-sort(unique(c(match(code,g$TS_PCODE),s$nb[[match(code,g$TS_PCODE)]])))
  candidates<-records %>% filter(township_id%in%g$TS_PCODE[ids]) %>%
   arrange(desc(fatalities),event_id_cnty)
  stopifnot(nrow(candidates)>0,candidates$fatalities[1]>0)
  # Remove one record, using its ID only in memory to break tied maxima.
  reduced<-records %>% filter(event_id_cnty!=candidates$event_id_cnty[1])
  stopifnot(nrow(reduced)==nrow(records)-1L)
  message('Largest-incident deletion: ',sample_name,', focal code ',code)
  out<-severity_test(reduced,sample_name,B=B,seed=seeds[1])
  deleted[[paste(sample_name,code)]]<-out %>% mutate(deletion_case=.env$code,
   supported=q<.05 & observed>expected)
  deletion_audit[[paste(sample_name,code)]]<-tibble(scenario=sample_name,
   deletion_case=code,removed_records=1L,removed_fatalities=candidates$fatalities[1],
   before_records=nrow(records),after_records=nrow(reduced),
   before_fatalities=sum(records$fatalities),after_fatalities=sum(reduced$fatalities))
 }
}
repeated<-bind_rows(repeated);deleted<-bind_rows(deleted)
deletion_audit<-bind_rows(deletion_audit)
cases<-h$robust_cases %>% select(code,township,region)
seed_cases<-repeated %>% inner_join(cases,by='code') %>%
 group_by(code,township,scenario) %>% summarise(runs=n(),supported_runs=sum(supported),
 ratio=first(ratio),p_min=min(p),p_max=max(p),q_min=min(q),q_max=max(q),.groups='drop')
seed_families<-repeated %>% group_by(scenario,seed) %>% summarise(
 supported_neighbourhoods=sum(supported),.groups='drop')
deletion_cases<-deleted %>% filter(code==deletion_case) %>% inner_join(cases,by='code') %>%
 left_join(deletion_audit,by=c('scenario','deletion_case')) %>%
 select(code,township,scenario,observed,expected,ratio,p,q,supported,removed_fatalities)
stopifnot(nrow(seed_cases)==6,nrow(deletion_cases)==6,
 identical(unname(tools::md5sum(names(d$quality$hashes))),unname(d$quality$hashes)))
saveRDS(list(seeds=seeds,B=B,repeated=repeated,seed_cases=seed_cases,
 seed_families=seed_families,deleted=deleted,deletion_cases=deletion_cases,
 deletion_audit=deletion_audit),'outputs/severity_validation.rds')
write.csv(seed_cases,'outputs/severity_seed_checks.csv',row.names=FALSE)
write.csv(deletion_cases,'outputs/severity_deletion_checks.csv',row.names=FALSE)
print(seed_cases);print(seed_families);print(deletion_cases)
