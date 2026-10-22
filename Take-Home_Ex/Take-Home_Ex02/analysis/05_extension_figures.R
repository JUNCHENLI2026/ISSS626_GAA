suppressPackageStartupMessages({library(ggplot2);library(sf);library(dplyr);library(tidyr)})
d<-readRDS('outputs/prepared.rds');s<-readRDS('outputs/statistics.rds');e<-readRDS('outputs/extensions.rds')
g<-d$townships
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='bottom',
 plot.title=element_text(face='bold'),plot.title.position='plot',plot.caption.position='plot',
 plot.background=element_rect(fill='white',colour=NA)))
save_plot<-function(p,name,w=10,h=7)ggsave(paste0('outputs/figures/',name,'.png'),p,width=w,height=h,dpi=160,bg='white')
map_theme<-theme(axis.title=element_blank(),axis.text=element_blank(),panel.grid=element_blank(),legend.title=element_blank())
regions<-g %>% group_by(ST) %>% summarise(.groups='drop')
case_points<-st_point_on_surface(g %>% inner_join(e$cases %>% select(code,case_id),by=c('TS_PCODE'='code')))

m<-g %>% left_join(e$checks,by=c('TS_PCODE'='code'))
p<-ggplot(m)+geom_sf(aes(fill=evidence),colour='white',linewidth=.08)+
 geom_sf(data=regions,fill=NA,colour='#666666',linewidth=.15)+
 facet_wrap(~variable,labeller=as_labeller(c(events='Events',fatalities='Reported fatalities')))+
 scale_fill_manual(values=c('No displayed trend-pattern'='#eeeeee','Sensitivity-dependent pattern'='#e5a665','BH and all block checks'='#166675'))+
 guides(fill=guide_legend(ncol=1))+map_theme+labs(title='Strength of evidence for displayed EHSA patterns',
 subtitle='Same corrected class definitions. Trend checks use 999 block permutations.',
 caption='Additional support means MK BH q < .05 and p < .05 for 3-, 6- and 12-month blocks.\nThis is a sensitivity screen, not a new joint significance test.')
save_plot(p,'trend_robustness',10,8)

m<-g %>% left_join(e$concordance,by=c('TS_PCODE'='code'),suffix=c('','.derived'))
p<-ggplot(m)+geom_sf(aes(fill=category),colour='white',linewidth=.08)+
 geom_sf(data=regions,fill=NA,colour='#666666',linewidth=.2)+
 geom_sf_text(data=case_points,aes(label=case_id),size=3,colour='#111111',fontface='bold')+
 scale_fill_manual(values=c('Both outcomes'='#773f8f','Events only'='#277f8e','Fatalities only'='#d57c36','Neither outcome'='#eeeeee'))+
 map_theme+labs(title='Frequency and fatality hot spots overlap only partly',
 subtitle='Positive Gi* with BH q < 0.05, adjusted separately for each outcome',
 caption='Case labels identify the four township profiles. Grey does not imply safety.')
save_plot(p,'outcome_concordance',8,9)

m<-g %>% left_join(e$severity,by=c('TS_PCODE'='code')) %>% mutate(label=case_when(
 excess & precise_excess~'Excess in both location subsets',excess~'Excess in full data only',TRUE~'No full-data excess at q < .05'))
p<-ggplot(m)+geom_sf(aes(fill=label),colour='white',linewidth=.08)+
 geom_sf(data=regions,fill=NA,colour='#666666',linewidth=.2)+
 scale_fill_manual(values=c('Excess in both location subsets'='#853251','Excess in full data only'='#e6b17b','No full-data excess at q < .05'='#eeeeee'))+
 guides(fill=guide_legend(ncol=1))+map_theme+labs(title='Neighbourhood fatality severity\nconditional on event composition',
 subtitle='999 random-label permutations within month and event type\nBH q < .05',
 caption='Expected totals condition on incident locations, months and types.\nA highlighted township refers to its self-including queen neighbourhood.')
save_plot(p,'conditional_severity',8,9)
save_plot(p+theme(plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),
 legend.text=element_text(size=14)), 'conditional_severity_slide',5.5,7)

# Two aligned count panels and two Gi* panels per case, without dual axes.
for(id in e$cases$case_id){
 a<-e$case_trends %>% filter(case_id==id)
 count<-a %>% transmute(month,value,significant=FALSE,
   measure=ifelse(variable=='events','1. Monthly events','2. Monthly reported fatalities'))
 z<-a %>% transmute(month,value=gi_z,significant=p<.05,
   measure=ifelse(variable=='events','3. Event Gi* z','4. Fatality Gi* z'))
 series<-bind_rows(count,z)
 p<-ggplot(series,aes(month,value))+geom_line(colour='#235677',linewidth=.6)+
   geom_point(data=filter(series,significant),colour='#b34542',size=1.3)+
   facet_wrap(~measure,ncol=2,scales='free_y')+
   scale_x_date(date_breaks='1 year',date_labels='%Y')+
   labs(x=NULL,y=NULL,title=paste0(id,': ',a$township[1],', ',a$region[1]),
     subtitle='Red points: space-time Gi* permutation p < .05. Each panel uses its own vertical scale.',
     caption='All 57 months appear, including zero-event months. Lines are descriptive, not forecasts.')
 save_plot(p,paste0('case_',gsub('Case ','',id)),11,6.5)
}

# A locator supplements the unlabelled national statistical maps.
loc<-st_transform(regions,4326)
centres<-st_point_on_surface(regions)
case_xy<-st_coordinates(case_points)
leaders<-data.frame(x=case_xy[,1],y=case_xy[,2],xend=case_xy[,1]+80000,yend=case_xy[,2])
p<-ggplot()+geom_sf(data=regions,fill='#f3f3f3',colour='#858585',linewidth=.3)+
 geom_sf(data=g %>% filter(TS_PCODE %in% e$cases$code),fill='#c46648',colour='white',linewidth=.15)+
 geom_sf_text(data=centres %>% filter(ST %in% setdiff(unique(e$cases$region),'Mandalay')),aes(label=ST),size=3,colour='#354d60')+
 geom_sf_text(data=centres %>% filter(ST=='Mandalay'),aes(label=ST),size=3,colour='#354d60',nudge_x=230000,nudge_y=-70000)+
 geom_segment(data=leaders,aes(x=x,y=y,xend=xend,yend=yend),colour='#555555',linewidth=.25)+
 geom_sf_label(data=case_points,aes(label=case_id),size=3,linewidth=.1,nudge_x=80000)+
 map_theme+labs(title='Township case-study locations',subtitle='State and region outlines give geographical context.')
save_plot(p,'case_locator',7,9)
