# Adapt the Hands-on 5(1) global test and 5(2) centred scatterplot workflow.
# Reuse the existing township order and repaired queen graph; do not change EHSA.
suppressPackageStartupMessages({library(spdep);library(sf);library(dplyr);library(ggplot2)})
d <- readRDS('outputs/prepared.rds')
s <- readRDS('outputs/statistics.rds')
g <- d$townships
stopifnot(identical(as.character(attr(s$nb,'region.id')),as.character(g$TS_PCODE)),
          all(card(s$nb)>0))
w <- nb2listw(s$nb,style='W',zero.policy=FALSE)
n <- nrow(g)
B <- 9999L
rows <- points <- list()
tests <- list()
for(v in c('events','fatalities')) {
  x <- g[[v]]
  stopifnot(length(x)==n,all(is.finite(x)),sd(x)>0)
  z <- as.numeric(scale(x))
  lag_z <- lag.listw(w,z,zero.policy=FALSE)
  seed <- if(v=='events') 2026102101L else 2026102102L
  set.seed(seed)
  mc <- moran.mc(x,w,nsim=B,alternative='two.sided',zero.policy=FALSE)
  # Explicit inclusive tails with a +1 correction avoid a zero Monte Carlo p.
  simulated <- mc$res[seq_len(B)]
  observed <- as.numeric(mc$statistic)
  permutation_p <- min(1,2*min((1+sum(simulated>=observed))/(B+1),
                              (1+sum(simulated<=observed))/(B+1)))
  # Independent algebraic check: S0 = n for these non-isolated W weights.
  slope <- sum(z*lag_z)/sum(z*z)
  direct_I <- (n/Szero(w))*slope
  stopifnot(abs(Szero(w)-n)<1e-10,abs(direct_I-as.numeric(mc$statistic))<1e-12)
  local <- s$static %>% filter(variable==v)
  q <- local$moran_q[match(g$TS_PCODE,local$code)]
  quadrant <- case_when(z>0 & lag_z>0~'High-High',z<0 & lag_z<0~'Low-Low',
                        z>0 & lag_z<0~'High-Low',z<0 & lag_z>0~'Low-High',TRUE~'On mean')
  rows[[v]] <- data.frame(variable=v,n=n,I=direct_I,null_expectation=-1/(n-1),
                          permutation_p=permutation_p,permutations=B,seed=seed,slope=slope)
  points[[v]] <- data.frame(code=g$TS_PCODE,variable=v,z=z,lag_z=lag_z,
                            quadrant=quadrant,local_q=q,
                            local_support=ifelse(q<.05,'Local BH q < 0.05','Other townships'))
  tests[[v]] <- mc
}
results <- bind_rows(rows)
results$holm_p <- p.adjust(results$permutation_p,method='holm')
# Cross two representations with two graphs for each outcome. Keep all eight
# comparisons in one Holm family, including the original raw/queen results.
xy <- st_coordinates(st_point_on_surface(g))
kn6 <- make.sym.nb(knn2nb(knearneigh(xy,k=6),row.names=g$TS_PCODE))
graphs <- list(Queen=s$nb,`Symmetric 6-NN`=kn6)
sensitivity_rows <- list()
for(v in c('events','fatalities')) for(graph_name in names(graphs)) for(transform in c('Raw','log1p')) {
  x <- if(transform=='Raw') g[[v]] else log1p(g[[v]])
  sw <- nb2listw(graphs[[graph_name]],style='W',zero.policy=FALSE)
  base <- results[results$variable==v,]
  # The same seed within an outcome pairs the random reallocations across settings.
  set.seed(base$seed)
  mc <- moran.mc(x,sw,nsim=B,alternative='two.sided',zero.policy=FALSE)
  null <- mc$res[seq_len(B)]
  obs <- as.numeric(mc$statistic)
  pv <- min(1,2*min((1+sum(null>=obs))/(B+1),(1+sum(null<=obs))/(B+1)))
  centered <- x-mean(x)
  direct <- n/Szero(sw)*sum(centered*lag.listw(sw,centered))/sum(centered^2)
  stopifnot(abs(obs-direct)<1e-12)
  sensitivity_rows[[length(sensitivity_rows)+1L]] <- data.frame(variable=v,
    graph=graph_name,transform=transform,I=obs,permutation_p=pv,
    permutations=B,seed=base$seed,min_neighbours=min(card(graphs[[graph_name]])),
    max_neighbours=max(card(graphs[[graph_name]])))
}
sensitivity <- bind_rows(sensitivity_rows)
sensitivity$holm_p_8 <- p.adjust(sensitivity$permutation_p,'holm')
sensitivity$positive_supported <- sensitivity$I>0 & sensitivity$holm_p_8<.05
scatter <- bind_rows(points)
labels <- setNames(sprintf('%s\nGlobal I = %.3f',c('Event counts','Reported fatalities'),results$I),results$variable)
palette <- c('High-High'='#b2182b','Low-Low'='#2166ac','High-Low'='#ef8a62','Low-High'='#67a9cf','On mean'='#777777')
p <- ggplot(scatter,aes(z,lag_z))+
  geom_hline(yintercept=0,colour='grey65',linetype='dashed')+
  geom_vline(xintercept=0,colour='grey65',linetype='dashed')+
  geom_abline(data=results,aes(slope=slope,intercept=0),colour='#222222',linewidth=.65)+
  geom_point(aes(colour=quadrant,shape=local_support),size=2.3,alpha=.8)+
  facet_wrap(~variable,labeller=as_labeller(labels))+
  scale_colour_manual(values=palette)+
  scale_shape_manual(values=c('Local BH q < 0.05'=17,'Other townships'=1))+
  labs(title='Overall spatial association and local quadrants',
       subtitle='330 townships; January 2021 to September 2025; unchanged queen-based W weights',
       x='Standardised township total (z)',y='Spatial lag of the same z',colour='Descriptive quadrant',shape='Existing local Moran test',
       caption='Lines pass through the origin with slope I. Colour alone does not imply significance.\nTriangles use existing outcome-specific BH-adjusted local tests, not the global p-values.')+
  theme_minimal(base_size=12)+theme(legend.position='bottom',panel.grid.minor=element_blank(),
    plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0))+
  guides(colour=guide_legend(nrow=1,order=1),shape=guide_legend(nrow=1,order=2))
dir.create('outputs/figures',recursive=TRUE,showWarnings=FALSE)
ggsave('outputs/figures/global_moran_scatter.png',p,width=11,height=7,dpi=160,bg='white')
saveRDS(list(results=results,scatter=scatter,tests=tests,sensitivity=sensitivity,
            sensitivity_graphs=graphs),'outputs/global_moran.rds')
write.csv(results,'outputs/global_moran_summary.csv',row.names=FALSE)
write.csv(sensitivity,'outputs/global_moran_sensitivity.csv',row.names=FALSE)
print(results)
print(sensitivity)
