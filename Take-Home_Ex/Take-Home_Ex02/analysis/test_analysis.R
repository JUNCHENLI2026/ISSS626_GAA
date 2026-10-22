# Run from the exercise directory: Rscript --vanilla analysis/test_analysis.R
source('analysis/report_helpers.R')
test_months <- as.Date(c('2023-11-01', '2025-04-01', '2025-06-01', '2022-10-01', NA))
expected_months <- c('November 2023', 'April 2025', 'June 2025', 'October 2022', NA_character_)
stopifnot(identical(format_month_en(test_months), expected_months))
local({
  original_locale <- Sys.getlocale('LC_TIME')
  on.exit(Sys.setlocale('LC_TIME', original_locale))
  for (candidate in c('C', 'Chinese (Simplified)_China.utf8', 'zh_CN.UTF-8')) {
    changed <- suppressWarnings(Sys.setlocale('LC_TIME', candidate))
    if (nzchar(changed)) stopifnot(identical(format_month_en(test_months), expected_months))
  }
})
suppressPackageStartupMessages({library(spdep);library(sf);library(dplyr)})
for(expr in parse('analysis/02_statistics.R')) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
    as.character(expr[[2]])=='classify_pattern') eval(expr)
}
test_pattern<-function(z,expected,tau=.3,p=.2){
 pp<-ifelse(z==0,1,.01)
 stopifnot(identical(classify_pattern(z,pp,tau,p),expected))
}
for(sign in c(1,-1)) {
 label<-if(sign==1)'hot spot' else 'cold spot'
 test_pattern(sign*c(rep(0,9),1),paste('New',label))
 test_pattern(sign*c(rep(0,8),1,1),paste('Consecutive',label))
 test_pattern(sign*c(0,1,rep(0,7),1),paste('Sporadic',label))
 test_pattern(sign*c(-1,rep(0,8),1),paste('New',label))
 test_pattern(sign*c(1,-1,rep(0,7),1),paste('Oscillating',label))
 test_pattern(sign*rep(1,10),paste('Persistent',label))
 test_pattern(sign*c(0,rep(1,9)),paste('Intensifying',label),tau=sign*.3,p=.01)
 test_pattern(sign*c(0,rep(1,9)),paste('Diminishing',label),tau=-sign*.3,p=.01)
 test_pattern(sign*c(rep(1,9),0),paste('Historical',label))
 test_pattern(sign*c(rep(1,8),0,0),'No pattern detected')
 test_pattern(sign*c(0,1,rep(0,8)),'No pattern detected')
}
test_pattern(rep(0,10),'No pattern detected')
stopifnot(classify_pattern(rep(1,10),rep(.05,10),0,1)=='No pattern detected')

d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds');e<-readRDS('outputs/extensions.rds')
stopifnot(nrow(d$panel)==330*57,!anyDuplicated(d$panel[c('township_id','month')]),
 sum(d$panel$events)==nrow(d$assigned),sum(d$panel$fatalities)==sum(d$assigned$fatalities),
 identical(unname(tools::md5sum(names(d$quality$hashes))),unname(d$quality$hashes)))
for(v in c('events','fatalities')) {
 gi<-s$gis[s$gis$variable==v,]
 stopifnot(identical(gi$township_id,d$panel$township_id),identical(gi$month,d$panel$month))
 a<-s$ehsa[s$ehsa$variable==v,]
 hot<-a$classification=='Sporadic hot spot';cold<-a$classification=='Sporadic cold spot'
 stopifnot(all(a$final_p[hot|cold]<.05),all(a$final_gi[hot]>0),all(a$final_gi[cold]<0),
 all(a$mk_p[a$map_class!='Not displayed']<.05))
}
stopifnot(all(e$checks$trend_supported==(e$checks$mk_p<.05 & e$checks$mk_q<.05 &
 e$checks$block_3<.05 & e$checks$block_6<.05 & e$checks$block_12<.05)),
 nrow(e$cases)==4,!anyDuplicated(e$cases$code),nrow(e$case_trends)==4*57*2,
 all(e$severity$p>=.001 & e$severity$p<=1),all(e$severity$q>=e$severity$p-1e-12))
stopifnot(all.equal(e$severity$q,p.adjust(e$severity$p,'BH')))

# Small-graph random-labelling checks: exchangeable identical marks give p = 1.
test_env<-new.env(parent=globalenv())
test_env$n<-3L;test_env$g<-data.frame(TS_PCODE=c('A','B','C'))
test_env$s<-list(nb=cell2nb(1,3,type='rook'))
for(expr in parse('analysis/04_extensions.R')) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
    as.character(expr[[2]])=='severity_test') eval(expr,envir=test_env)
}
tiny<-data.frame(township_id=c('A','B','C'),month=as.Date('2021-01-01'),event_type='Battles',fatalities=2)
identical_marks<-test_env$severity_test(tiny,'Test',B=99L)
stopifnot(all(identical_marks$p==1),all(identical_marks$observed==identical_marks$expected))
tiny$fatalities<-c(0,0,9)
reference_marks<-test_env$severity_test(tiny,'Test',B=9999L)
# Each end neighbourhood contains two of three events: P(sum >= 9) = 2/3.
stopifnot(reference_marks$p[1]==1,reference_marks$p[2]==1,
 abs(reference_marks$p[3]-2/3)<.03,
 isTRUE(all.equal(reference_marks$expected,c(6,9,6))))

# Independent small-graph check of the conditional Gi* sampler.
source('analysis/fast_gi.R')
x<-c(1,2,4,8,16,32,64,128)
w<-nb2listw(include.self(cell2nb(2,4,type='rook')),style='W')
a<-fast_space_time_gi(x,w,nsim=19999L)
exact<-sapply(seq_along(x),function(i){
 ids<-w$neighbours[[i]];k<-length(ids);observed<-sum(x[ids])/k
 vals<-colSums(combn(x[-i],k-1))/k+x[i]/k
 min(1,2*min(mean(vals<=observed+1e-9),mean(vals>=observed-1e-9)))
})
stopifnot(max(abs(a$p-exact))<.03,isTRUE(all.equal(a$z,as.numeric(localG(x,w)))))
cat('PASS: classification boundaries, final-bin rule, raw hashes, cube totals, ordering, map masks, robustness and conditional sampler.\n')

# Evidence-audit contracts: reconcile records, set arithmetic and benchmarks.
h<-readRDS('outputs/evidence_audit.rds')
stopifnot(sum(h$audit_summary$records)==nrow(d$assigned),
 sum(h$name_pairs$records)==sum(h$audit_summary$records[h$audit_summary$audit_class!='Exact name match']),
 sum(h$audit_summary$records[h$audit_summary$audit_class!='Exact name match' &
 h$audit_summary$audit_class!='Missing township name'])==d$quality$name_check$exact_name_disagreement,
 all(h$audit_stability$baseline_hot==h$audit_stability$retained+h$audit_stability$lost),
 all(h$audit_stability$alternative_hot==h$audit_stability$retained+h$audit_stability$gained),
 all(h$audit_stability$Jaccard>=0 & h$audit_stability$Jaccard<=1))
for(v in c('events','fatalities')){
 a<-h$relative_comparison %>% filter(variable==v)
 stopifnot(abs(sum(a$share_2024)-100)<1e-8,abs(sum(a$share_2025)-100)<1e-8,
 abs(sum(a$expected_2025)-unique(a$national_2025))<1e-8,
 sum(a$local_2024)==unique(a$national_2024),sum(a$local_2025)==unique(a$national_2025))
}
stopifnot(nrow(h$robust_cases)==3,all(h$robust_cases$excess & h$robust_cases$precise_excess),
 all(h$overlap$Jaccard<=.5),nrow(h$case_shares)==4*19*2)
for(code in h$robust_cases$code){
 a<-h$robust_cases %>% filter(.data$code==.env$code)
 b<-h$decomposition %>% filter(.data$code==.env$code)
 c<-h$annual %>% filter(.data$code==.env$code)
 stopifnot(abs(sum(b$observed)-a$observed)<1e-8,abs(sum(b$expected)-a$expected)<1e-8,
 abs(sum(c$observed)-a$observed)<1e-8,abs(sum(c$expected)-a$expected)<1e-8)
}
cat('PASS: name-audit reconciliation, hotspot set arithmetic, national benchmarks and severity decompositions.\n')

# Validation keeps the original selected cases and the full BH family.
j<-readRDS('outputs/severity_validation.rds')
stopifnot(nrow(j$repeated)==330*2*3,nrow(j$deleted)==330*2*3,
 nrow(j$seed_cases)==6,nrow(j$deletion_cases)==6,
 setequal(j$seed_cases$code,h$robust_cases$code),
 all(j$deletion_audit$before_records-j$deletion_audit$after_records==1),
 all(j$deletion_audit$before_fatalities-j$deletion_audit$after_fatalities==
     j$deletion_audit$removed_fatalities))
for(sample_name in unique(j$repeated$scenario)){
 for(seed in j$seeds){
  z<-j$repeated %>% filter(scenario==sample_name,.data$seed==.env$seed)
  stopifnot(nrow(z)==330,!anyDuplicated(z$code),
   isTRUE(all.equal(z$q,p.adjust(z$p,'BH'))),
   all(z$p>=1/(j$B+1)),all(z$supported==(z$q<.05 & z$observed>z$expected)))
 }
 for(code in h$robust_cases$code){
  z<-j$deleted %>% filter(scenario==sample_name,deletion_case==.env$code)
  stopifnot(nrow(z)==330,isTRUE(all.equal(z$q,p.adjust(z$p,'BH'))))
 }
}
cat('PASS: repeated-seed and largest-record-deletion families, p-values and deletion accounting.\n')

# Supplementary global tests reuse exactly the same local-analysis geography.
gm <- readRDS('outputs/global_moran.rds')
stopifnot(nrow(gm$results)==2,nrow(gm$scatter)==660,
 all(gm$results$permutations==9999L),all(gm$results$n==330),
 all(gm$results$permutation_p>0 & gm$results$permutation_p<=1),
 isTRUE(all.equal(gm$results$holm_p,p.adjust(gm$results$permutation_p,'holm'))))
gw <- spdep::nb2listw(s$nb,style='W',zero.policy=FALSE)
for(v in c('events','fatalities')) {
 a <- gm$scatter[gm$scatter$variable==v,]
 r <- gm$results[gm$results$variable==v,]
 x <- d$townships[[v]]
 centered <- x-mean(x)
 manual_I <- length(x)/spdep::Szero(gw)*
   sum(centered*spdep::lag.listw(gw,centered))/sum(centered^2)
 stopifnot(identical(as.character(a$code),as.character(d$townships$TS_PCODE)),
   max(abs(a$z-as.numeric(scale(x))))<1e-12,
   max(abs(a$lag_z-spdep::lag.listw(gw,a$z)))<1e-12,
   abs(r$I-manual_I)<1e-12,abs(r$slope-r$I)<1e-12,
   abs(r$null_expectation+1/329)<1e-12)
 old <- s$static[s$static$variable==v,]
 old <- old[match(a$code,old$code),]
 stopifnot(isTRUE(all.equal(a$local_q,unname(old$moran_q))),
   all(a$quadrant[old$moran_p<.05]==old$lisa[old$moran_p<.05]))
 set.seed(r$seed)
 repeat_test <- spdep::moran.mc(x,gw,nsim=9999,alternative='two.sided')
 null <- repeat_test$res[seq_len(9999)]
 repeated_p <- min(1,2*min((1+sum(null>=manual_I-1e-14))/10000,
                         (1+sum(null<=manual_I+1e-14))/10000))
 stopifnot(abs(repeated_p-r$permutation_p)<1e-12,
   abs(as.numeric(repeat_test$statistic)-r$I)<1e-12)
}
cat('PASS: global Moran formula, township order, scatterplot slope, local support and repeatable permutations.\n')

gs <- gm$sensitivity
stopifnot(nrow(gs)==8,!anyDuplicated(gs[c('variable','graph','transform')]),
 all(gs$permutations==9999),all(gs$permutation_p>=.0002 & gs$permutation_p<=1),
 isTRUE(all.equal(gs$holm_p_8,p.adjust(gs$permutation_p,'holm'))),
 all(gs$positive_supported==(gs$I>0 & gs$holm_p_8<.05)),
 identical(gm$sensitivity_graphs$Queen,s$nb))
xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(d$townships)))
expected_kn <- spdep::make.sym.nb(spdep::knn2nb(spdep::knearneigh(xy,k=6),row.names=d$townships$TS_PCODE))
actual_kn <- gm$sensitivity_graphs[['Symmetric 6-NN']]
# Construction-call attributes differ between independent code paths; compare
# every edge and P-code rather than the stored language expression.
stopifnot(identical(attr(actual_kn,'region.id'),attr(expected_kn,'region.id')),
 all(vapply(seq_along(expected_kn),function(i)identical(actual_kn[[i]],expected_kn[[i]]),logical(1))))
for(k in seq_len(nrow(gs))) {
 row <- gs[k,]
 nb <- gm$sensitivity_graphs[[row$graph]]
 stopifnot(all(spdep::card(nb)>0),
  all(vapply(seq_along(nb),function(i) !(i %in% nb[[i]]),logical(1))))
 if(row$graph=='Symmetric 6-NN')stopifnot(all(vapply(seq_along(nb),
   function(i)all(vapply(nb[[i]],function(j)i %in% nb[[j]],logical(1))),logical(1))))
 lw <- spdep::nb2listw(nb,style='W')
 x <- d$townships[[row$variable]]
 if(row$transform=='log1p')x <- log1p(x)
 centered <- x-mean(x)
 direct <- length(x)/spdep::Szero(lw)*sum(centered*spdep::lag.listw(lw,centered))/sum(centered^2)
 stopifnot(abs(row$I-direct)<1e-12,
  row$min_neighbours==min(spdep::card(nb)),row$max_neighbours==max(spdep::card(nb)))
 set.seed(row$seed)
 rerun <- spdep::moran.mc(x,lw,nsim=9999,alternative='two.sided')
 null <- rerun$res[seq_len(9999)]
 observed <- as.numeric(rerun$statistic)
 pv <- min(1,2*min((1+sum(null>=observed))/10000,(1+sum(null<=observed))/10000))
 stopifnot(abs(pv-row$permutation_p)<1e-12)
 if(row$graph=='Queen' && row$transform=='Raw') {
  original <- gm$results[gm$results$variable==row$variable,]
  stopifnot(abs(original$I-row$I)<1e-12,original$permutation_p==row$permutation_p)
 }
}
cat('PASS: eight-setting global sensitivity, graph construction, Holm family and repeated permutations.\n')
