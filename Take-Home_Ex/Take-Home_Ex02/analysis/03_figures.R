suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(sf);library(scales)})
d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds');g<-d$townships
dir.create('outputs/figures',recursive=TRUE,showWarnings=FALSE)
theme_set(theme_minimal(base_size=12)+theme(plot.title=element_text(face='bold',size=16),plot.subtitle=element_text(size=11),panel.grid.minor=element_blank(),legend.position='bottom',plot.background=element_rect(fill='white',colour=NA)))
save_plot<-function(p,name,w=10,h=7){ggsave(file.path('outputs/figures',paste0(name,'.png')),p,width=w,height=h,dpi=160,bg='white')}
save_slide_map<-function(p,name){
  slide_plot<-p+theme(plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),legend.text=element_text(size=17),legend.key.size=grid::unit(.5,'cm'),plot.margin=margin(2,2,2,2))+guides(fill=guide_legend(ncol=1))
  save_plot(slide_plot,paste0(name,'_slide'),5.5,7)
}
map_theme<-theme(axis.title=element_blank(),axis.text=element_blank(),panel.grid.major=element_blank(),legend.title=element_blank())
# Selected administrative anchors orient readers without labelling every township.
regions<-g %>% group_by(ST) %>% summarise(.groups='drop')
anchors<-regions %>% filter(ST%in%c('Sagaing','Mandalay','Shan (North)','Kayin')) %>%
 st_point_on_surface()
geographic_anchors<-list(
 geom_sf(data=regions,fill=NA,colour='#666666',linewidth=.22,inherit.aes=FALSE),
 geom_sf_label(data=anchors,aes(label=ST),size=2.7,fill='white',
  colour='#252525',label.padding=grid::unit(.10,'lines'),inherit.aes=FALSE))
var_names<-c(events='Recorded armed-conflict events',fatalities='Reported fatalities')
desc<-g %>% select(TS_PCODE,events,fatalities) %>% pivot_longer(c(events,fatalities),names_to='variable',values_to='value')
p<-ggplot(desc)+geom_sf(aes(fill=value),colour='white',linewidth=.07)+facet_wrap(~variable,labeller=as_labeller(var_names))+scale_fill_viridis_c(trans='log1p',breaks=c(0,10,100,1000),labels=label_number(big.mark=','),name='Total (log-scaled colours)',guide=guide_colourbar(barwidth=grid::unit(10,'cm')))+map_theme+theme(legend.title=element_text(size=10))+labs(title='Conflict burden varies strongly across townships',subtitle='January 2021 to September 2025. Counts are not population-adjusted risks.')
save_plot(p,'burden',10,8)
monthly<-d$panel %>% group_by(month) %>% summarise(events=sum(events),fatalities=sum(fatalities),.groups='drop')
p<-monthly %>% pivot_longer(c(events,fatalities)) %>% ggplot(aes(month,value))+geom_line(colour='#214e75',linewidth=.8)+facet_wrap(~name,ncol=1,scales='free_y',labeller=as_labeller(var_names))+scale_y_continuous(labels=label_number(big.mark=','))+labs(x=NULL,y=NULL,title='Monthly recorded conflict and reported fatalities',subtitle='Complete monthly bins. The 2025 series ends in September.')
save_plot(p,'monthly',10,6)
pal<-c('High-High'='#b2182b','Low-Low'='#2166ac','High-Low'='#ef8a62','Low-High'='#67a9cf','Hot spot'='#b2182b','Cold spot'='#2166ac')
for(v in c('events','fatalities'))for(kind in c('lisa','hotspot')){
  m<-g %>% left_join(s$static %>% filter(variable==v),by=c('TS_PCODE'='code'))
  active<-m[m[[kind]]!='Not significant',]
  p<-ggplot()+geom_sf(data=g,fill='#ededed',colour='white',linewidth=.08)+geom_sf(data=active,aes(fill=.data[[kind]]),colour='white',linewidth=.08)+scale_fill_manual(values=pal,drop=TRUE)+map_theme+labs(title=paste(if(kind=='lisa')'Local Moran clusters and outliers' else 'Getis-Ord Gi* hot and cold spots'),subtitle=paste(var_names[v], '| two-sided permutation p < 0.05'),caption='Grey: not significant at 0.05. Main maps are unadjusted exploratory tests.')
  save_plot(p+geographic_anchors,paste(v,kind,sep='_'),7,8)
  if((v=='events' && kind=='lisa') || (v=='fatalities' && kind=='hotspot'))save_slide_map(p,paste(v,kind,sep='_'))
}
ehpal<-c('New hot spot'='#67000d','Consecutive hot spot'='#a50f15','Intensifying hot spot'='#cb181d','Persistent hot spot'='#ef3b2c','Diminishing hot spot'='#fb6a4a','Sporadic hot spot'='#fc9272','Oscillating hot spot'='#990099','Historical hot spot'='#fcbba1',
 'New cold spot'='#08306b','Consecutive cold spot'='#08519c','Intensifying cold spot'='#2171b5','Persistent cold spot'='#4292c6','Diminishing cold spot'='#6baed6','Sporadic cold spot'='#9ecae1','Oscillating cold spot'='#756bb1','Historical cold spot'='#c6dbef')
for(v in c('events','fatalities')){
  m<-g %>% left_join(s$ehsa %>% filter(variable==v),by=c('TS_PCODE'='code'))
  active<-m %>% filter(map_class!='Not displayed')
  p<-ggplot()+geom_sf(data=g,fill='#ededed',colour='white',linewidth=.08)+geom_sf(data=active,aes(fill=map_class),colour='white',linewidth=.08)+scale_fill_manual(values=ehpal)+guides(fill=guide_legend(ncol=2))+map_theme+labs(title='Emerging hot and cold spot patterns',subtitle=paste(var_names[v], '| Mann-Kendall p < 0.05'),caption='Pattern bins use Gi* permutation p < 0.05.\nGrey: no displayed significant trend-pattern combination.')
  save_plot(p+geographic_anchors,paste0(v,'_ehsa'),8,9)
  save_slide_map(p,paste0(v,'_ehsa'))
}
top<-st_drop_geometry(g) %>% slice_max(events,n=20,with_ties=FALSE) %>% arrange(events)
heat<-s$gis %>% filter(variable=='events',township_id %in% top$TS_PCODE) %>% mutate(township=factor(township_id,levels=top$TS_PCODE,labels=paste(top$TS,top$ST,sep=', ')))
p<-ggplot(heat,aes(month,township,fill=gi_z))+geom_tile()+scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',midpoint=0,name='Gi* z')+labs(x=NULL,y=NULL,title='Space-time clustering changes within high-burden townships',subtitle='Twenty largest event totals. Colours show Gi* strength, not significance.')+theme(panel.grid=element_blank(),axis.text.y=element_text(size=9))
save_plot(p,'gi_heatmap',11,7)
robust<-s$static %>% filter(variable=='events') %>% mutate(label=case_when(gi_p>=.05~'Not significant',gi_q<.05~'Survives BH adjustment',TRUE~'Nominal only'))
m<-g %>% left_join(robust,by=c('TS_PCODE'='code'))
p<-ggplot()+geom_sf(data=g,fill='#ededed',colour='white',linewidth=.08)+geom_sf(data=m %>% filter(label!='Not significant'),aes(fill=label),colour='white',linewidth=.08)+scale_fill_manual(values=c('Survives BH adjustment'='#154360','Nominal only'='#f4b183'))+map_theme+labs(title='How much evidence survives multiple testing?',subtitle='Event-count Gi*: nominal p < 0.05 versus BH q < 0.05')
save_plot(p,'multiplicity',7,8)
save_slide_map(p,'multiplicity')
top10<-g %>% st_drop_geometry() %>% slice_max(events,n=10,with_ties=FALSE) %>% arrange(events)
p<-ggplot(top10,aes(events,reorder(paste(TS,ST,sep=', '),events)))+geom_col(fill='#214e75')+scale_x_continuous(labels=label_number(big.mark=','))+labs(x='Recorded events',y=NULL,title='Ten largest township event totals',subtitle='Absolute burden over 57 months; not a ranking of individual risk.')
save_plot(p,'top_townships',10,6)
saveRDS(list(monthly=monthly,top10=top10), 'outputs/figure_data.rds')
