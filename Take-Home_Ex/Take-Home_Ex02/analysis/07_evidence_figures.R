suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(sf)})
h<-readRDS('outputs/evidence_audit.rds');d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
 plot.title=element_text(face='bold'),plot.title.position='plot',plot.caption.position='plot',
 plot.background=element_rect(fill='white',colour=NA),legend.position='bottom'))
save_plot<-function(p,name,w=10,h=7)ggsave(paste0('outputs/figures/',name,'.png'),p,width=w,height=h,dpi=160,bg='white')
for(v in c('events','fatalities')){
 z<-h$case_shares %>% filter(variable==v)
 p<-ggplot(z,aes(quarter,share_pct))+geom_line(colour='#235677',linewidth=.7)+geom_point(size=1,colour='#235677')+
  facet_wrap(~township,ncol=2,scales='fixed')+scale_x_date(date_breaks='1 year',date_labels='%Y')+
  labs(title=paste('Township share of national recorded',if(v=='events')'events' else 'fatalities'),
   subtitle='Quarterly totals divided by the corresponding national total',x=NULL,y='National share (%)',
   caption='All 19 quarters are complete. The series ends in Q3 2025.\nShares describe relative concentration, not population risk or an adjusted EHSA test.')
 save_plot(p,paste0('national_share_',v),10,6)
}
labels<-setNames(paste(h$robust_cases$case_label,h$robust_cases$township,sep=': '),h$robust_cases$code)
z<-h$decomposition %>% mutate(case=labels[code],event_type=ifelse(event_type=='Explosions/Remote violence','Explosions / remote violence',event_type))
p<-ggplot(z,aes(excess,event_type,fill=excess>0))+geom_col(width=.65)+geom_vline(xintercept=0,colour='#555555')+
 facet_wrap(~case,ncol=1)+scale_fill_manual(values=c('TRUE'='#287480','FALSE'='#b96e45'),guide='none')+
 labs(title='Which event types contribute to conditional excess fatalities?',
  subtitle='Observed minus expected reported fatalities in each selected neighbourhood',
  x='Observed minus expected fatalities',y=NULL,
  caption='Expected marks are calculated within month and event type nationally.\nBars are an additive descriptive decomposition, not separately tested effects.')
save_plot(p,'robust_severity_types',10,8)
z<-h$annual %>% mutate(case=labels[code],period=ifelse(year==2025,'2025 Jan-Sep',as.character(year)))
p<-ggplot(z,aes(period,excess))+geom_col(fill='#287480',width=.65)+geom_hline(yintercept=0,colour='#555555')+
 facet_wrap(~case,ncol=1)+labs(title='When did the conditional excess accumulate?',
  subtitle='Descriptive contributions to each neighbourhood total',x=NULL,y='Observed minus expected fatalities',
  caption='2021-2024 cover full calendar years; 2025 covers January-September only.\nThese contributions should not be compared as equal-duration annual rates.')
save_plot(p,'robust_severity_periods',10,8)
g<-d$townships;maps<-list()
for(i in seq_len(nrow(h$robust_cases))){
 code<-h$robust_cases$code[i];focal<-match(code,g$TS_PCODE);ids<-unique(c(focal,s$nb[[focal]]))
 maps[[i]]<-g %>% mutate(case=labels[[code]],area=case_when(TS_PCODE==code~'Focal township',row_number()%in%ids~'Included neighbour',TRUE~'Context'))
}
m<-bind_rows(maps);bbox<-st_bbox(m %>% filter(area!='Context'))
points<-g %>% filter(TS_PCODE%in%h$robust_cases$code) %>% st_point_on_surface() %>% mutate(case=labels[TS_PCODE])
p<-ggplot(m)+geom_sf(aes(fill=area),colour='#999999',linewidth=.15)+
 geom_sf_text(data=points,aes(label=TS),size=3,fontface='bold')+
 scale_fill_manual(values=c('Focal township'='#27727c','Included neighbour'='#b9d4d5','Context'='#f4f4f4'))+
 facet_wrap(~case,nrow=1)+coord_sf(xlim=c(bbox['xmin']-30000,bbox['xmax']+30000),ylim=c(bbox['ymin']-30000,bbox['ymax']+30000),expand=FALSE)+
 theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank(),legend.title=element_blank())+
 labs(title='Supported severity cases refer to neighbourhoods',subtitle='Local extent in Sagaing; the focal polygon is darker',
 caption='Neighbourhoods overlap. Their totals cannot be added and the three cases are not independent replications.')
save_plot(p,'robust_severity_neighbours',11,6)
