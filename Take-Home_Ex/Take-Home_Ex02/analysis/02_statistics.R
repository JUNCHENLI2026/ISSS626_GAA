# Local association and a one-month-lag space-time Gi* analysis.
suppressPackageStartupMessages({library(spdep);library(Kendall);library(dplyr);library(tidyr);library(sf)})
d <- readRDS('outputs/prepared.rds')
g <- d$townships
n <- nrow(g); nt <- length(d$months)
set.seed(20261018)
spdep::set.coresOption(1L)
nb0 <- poly2nb(g,queen=TRUE,row.names=g$TS_PCODE,snap=10)
xy <- st_coordinates(st_point_on_surface(g))
isolates <- which(card(nb0)==0L)
nb <- nb0
if(length(isolates)) {
  kn <- knn2nb(knearneigh(xy,k=1))
  for(i in isolates) {j <- kn[[i]][1];nb[[i]]<-as.integer(j);nb[[j]]<-sort(unique(as.integer(c(nb[[j]][nb[[j]]>0],i))))}
}
lw <- nb2listw(nb,style='W',zero.policy=TRUE)
gstar_w <- nb2listw(include.self(nb),style='W',zero.policy=TRUE)
kn6 <- make.sym.nb(knn2nb(knearneigh(xy,k=6),row.names=g$TS_PCODE))
weights_audit <- data.frame(township=g$TS,code=g$TS_PCODE,queen_neighbours=card(nb0),analysis_neighbours=card(nb),isolate_link_added=seq_len(n)%in%isolates)
local_results <- function(x,w,star,seed=20261018,nsim=4999L) {
  mi <- localmoran_perm(x,w,nsim=nsim,alternative='two.sided',iseed=seed,no_repeat_in_row=TRUE,zero.policy=TRUE)
  gi <- localG_perm(x,star,nsim=nsim,alternative='two.sided',iseed=seed+1L,no_repeat_in_row=TRUE,zero.policy=TRUE)
  z <- as.numeric(scale(x)); lagz <- lag.listw(w,z,zero.policy=TRUE)
  pmi <- mi[,'Pr(z != E(Ii)) Sim']; pgi <- attr(gi,'internals')[,'Pr(z != E(Gi)) Sim']
  tibble(Ii=mi[,'Ii'],moran_p=pmi,moran_q=p.adjust(pmi,'BH'),gi_z=as.numeric(gi),gi_p=pgi,gi_q=p.adjust(pgi,'BH'),
    lisa=ifelse(pmi<.05,ifelse(z>=0,ifelse(lagz>=0,'High-High','High-Low'),ifelse(lagz>=0,'Low-High','Low-Low')),'Not significant'),
    hotspot=ifelse(pgi<.05,ifelse(as.numeric(gi)>0,'Hot spot','Cold spot'),'Not significant'))
}
static <- bind_rows(lapply(c('events','fatalities'),function(v) bind_cols(tibble(code=g$TS_PCODE,variable=v),local_results(g[[v]],lw,gstar_w))))

# Time-major ordering: each monthly block contains the same 330 townships.
stnb <- vector('list',n*nt)
for(t in seq_len(nt)) for(i in seq_len(n)) {
  spatial <- sort(unique(c(i,nb[[i]][nb[[i]]>0])))
  here <- (t-1L)*n+i
  stnb[[here]] <- sort(as.integer(c((t-1L)*n+spatial,if(t>1L)(t-2L)*n+spatial else integer())))
}
class(stnb)<-'nb';attr(stnb,'region.id')<-as.character(seq_along(stnb));attr(stnb,'type')<-'space-time';attr(stnb,'sym')<-FALSE;attr(stnb,'self.included')<-TRUE
stw <- nb2listw(stnb,style='W',zero.policy=TRUE)
stopifnot(all(vapply(seq_along(stnb),function(i)i%in%stnb[[i]],logical(1))))

# Explicit, sign-specific pattern rules. No classification depends on place names.
classify_pattern <- function(z,p,tau,mkp) {
  hot<-p<.05 & z>0;cold<-p<.05 & z<0;m<-length(z)
  one_side <- function(a,opp,sign) {
    count<-sum(a);last<-a[m];share<-mean(a)
    if(count==0)return(NULL)
    if(last && count==1)return(paste('New',sign))
    if(share>=.9) {
      if(!last)return(paste('Historical',sign))
      if(mkp<.05) {
        growing <- if(sign=='hot spot') tau>0 else tau<0
        return(paste(if(growing)'Intensifying' else 'Diminishing',sign))
      }
      return(paste('Persistent',sign))
    }
    if(last) {
      run<-rle(rev(a))$lengths[1]
      if(run>=2 && run==count)return(paste('Consecutive',sign))
      if(any(opp))return(paste('Oscillating',sign))
    }
    # ArcGIS Pro 3.4 definitions require a same-sign significant final bin.
    if(last && !any(opp))return(paste('Sporadic',sign))
    NULL
  }
  a<-one_side(hot,cold,'hot spot');if(!is.null(a))return(a)
  a<-one_side(cold,hot,'cold spot');if(!is.null(a))return(a)
  'No pattern detected'
}
ehsa <- list(); gis <- list()
source('analysis/fast_gi.R')
for(v in c('events','fatalities')) {
  message('Space-time Gi*: ',v)
  x<-d$panel[[v]]
  gi<-fast_space_time_gi(x,stw,nsim=999L,seed=20261018L)
  zz<-gi$z; pp<-gi$p
  stopifnot(all(is.finite(zz)),all(is.finite(pp)),length(zz)==n*nt)
  gis[[v]]<-d$panel %>% select(township_id,month) %>% mutate(variable=v,gi_z=zz,p=pp,q=p.adjust(pp,'BH'))
  zm<-matrix(zz,nrow=n); pm<-matrix(pp,nrow=n)
  r<-lapply(seq_len(n),function(i) {
    mk<-MannKendall(zm[i,]); tau<-as.numeric(mk$tau); p<-as.numeric(mk$sl)
    if(!is.finite(p)){p<-1;tau<-0}
    # Three-month-block permutation preserves within-block temporal ordering.
    set.seed(20261018+i)
    blocks<-split(seq_len(nt),rep(seq_len(nt/3),each=3))
    obs<-tau
    null<-replicate(999,as.numeric(MannKendall(zm[i,unlist(blocks[sample(length(blocks))])])$tau))
    finite_null<-null[is.finite(null)]
    block_p<-if(length(finite_null)) (1+sum(abs(finite_null)>=abs(obs)))/(length(finite_null)+1) else 1
    tibble(code=g$TS_PCODE[i],variable=v,tau=tau,mk_p=p,block_p=block_p,hot_months=sum(pm[i,]<.05 & zm[i,]>0),cold_months=sum(pm[i,]<.05 & zm[i,]<0),
      final_gi=zm[i,nt],final_p=pm[i,nt],classification=classify_pattern(zm[i,],pm[i,],tau,p))
  })
  ehsa[[v]]<-bind_rows(r) %>% mutate(mk_q=p.adjust(mk_p,'BH'),map_class=ifelse(mk_p<.05 & classification!='No pattern detected',classification,'Not displayed'))
}
ehsa<-bind_rows(ehsa);gis<-bind_rows(gis)

# Sensitivity: alternative spatial neighbours, geolocation precision, and near-boundary records.
precise<-d$assigned %>% filter(geo_precision==1) %>% count(township_id)
x_precise<-replace_na(precise$n[match(g$TS_PCODE,precise$township_id)],0)
near<-d$unassigned %>% filter(match_count==0,nearest_m<=1000) %>% count(nearest_id)
x_near<-g$events+replace_na(near$n[match(g$TS_PCODE,near$nearest_id)],0)
sensitivity<-bind_rows(
 bind_cols(tibble(scenario='Six-nearest-neighbour graph',code=g$TS_PCODE),local_results(g$events,nb2listw(kn6,style='W'),nb2listw(include.self(kn6),style='W'))),
 bind_cols(tibble(scenario='Precisely located events only',code=g$TS_PCODE),local_results(x_precise,lw,gstar_w)),
 bind_cols(tibble(scenario='Add records within 1 km of boundary',code=g$TS_PCODE),local_results(x_near,lw,gstar_w)))
baseline<-static %>% filter(variable=='events')
sensitivity_summary<-sensitivity %>% group_by(scenario) %>% summarise(lisa_agreement=mean(lisa==baseline$lisa),hotspot_agreement=mean(hotspot==baseline$hotspot),hotspots=sum(hotspot=='Hot spot'),.groups='drop')
inference<-static %>% group_by(variable) %>% summarise(lisa_raw=sum(moran_p<.05),lisa_BH=sum(moran_q<.05),gi_raw=sum(gi_p<.05),gi_BH=sum(gi_q<.05),.groups='drop')
ehsa_summary<-ehsa %>% group_by(variable) %>% summarise(MK_raw=sum(mk_p<.05),MK_BH=sum(mk_q<.05),block_supported=sum(mk_p<.05 & block_p<.05),displayed=sum(map_class!='Not displayed'),.groups='drop')
saveRDS(list(static=static,ehsa=ehsa,gis=gis,weights=weights_audit,nb=nb,queen_components=n.comp.nb(nb0)$nc,analysis_components=n.comp.nb(nb)$nc,
 sensitivity=sensitivity,sensitivity_summary=sensitivity_summary,inference=inference,ehsa_summary=ehsa_summary), 'outputs/statistics.rds')
write.csv(static,'outputs/local_statistics.csv',row.names=FALSE)
write.csv(ehsa,'outputs/ehsa_statistics.csv',row.names=FALSE)
print(inference);print(ehsa_summary);print(sensitivity_summary)
print(table(ehsa$variable,ehsa$map_class))
