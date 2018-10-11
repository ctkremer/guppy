

# (needs other data from guppy_glmm_paper_figures.R)


# Load Longevity data:

dat<-read.csv("/Users/colin/Research/Active/guppy/FocalFish_fitness.csv")
dat<-dat[-which(dat$FishID_nodash=='FCA1O6G8B1005' & dat$GenGrp=='F1xN'),]  #drop incorrect row entry
dat3<-dat[!is.na(dat$Longevity),]
dat3$cohort.fctr<-as.factor(dat3$cohort)
dat3$hindex2<-dat3$hindex^2

# shift longevity variable to a [0,inf) range, to make it easier to check for inflation (zero inflation)
dat3$Longevity<-dat3$Longevity-1
dat3C<-dat3[dat3$Stream=='Caigual',]
dat3C$zero<-ifelse(dat3C$Longevity==0,0,1)
dat3T<-dat3[dat3$Stream=='Taylor',]
dat3T$zero<-ifelse(dat3T$Longevity==0,0,1)




# generate final plots:
tmp<-c3CT2[c3CT2$stream=='Caigual',]
tmp.tab<-tab32[tab32$stream=='Caigual',]
tmp.tab2<-tmp.tab
tmp.tab2$fit<-14.5
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$fit<-14.5
pA<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_point(data=dat3C,aes(x=hindex,y=LRS,colour=Sex),size=0.8,alpha=0.4)+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0,alpha=0.2,size=3)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_line(data=tmp.tab2,aes(x=hindex,y=fit,colour=Sex),linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0.01,0))+
  scale_y_continuous('Longevity (months)',expand=c(0.01,0),limits=c(0,15.5))+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  scale_shape_discrete(guide=FALSE)+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype=1,colour='black',size=0.3),
        plot.margin = margin(2, 15, 0, 2, "pt"),
        legend.position = 'none')+
  #legend.position = c(0.1,0.8),
  #legend.key.size = unit(0.5,'cm'))+
  ggtitle('A.')
pA



tmp<-c3CT2[c3CT2$stream=='Taylor' & c3CT2$parameter=='mu',]
tmp.tab<-tab32[tab32$stream=='Taylor'& tab32$parameter=='mu',]
tmp.tab2<-tmp.tab[tmp.tab$Sex=='M',]
tmp.tab2$fit<-14.5
tmp.tab2<-rbind(tmp.tab[tmp.tab$Sex=='M',],tmp.tab2)
tmp.tab$fit<-14.5
pB2<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_point(data=dat3T,aes(x=hindex,y=LRS,colour=Sex),size=0.8,alpha=0.4)+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=fit),colour=gray(0.7),height=0,alpha=0.2,size=3)+
  geom_line(aes(colour=Sex),size=0.8)+
  #  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
  #  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.35,size=0.4)+
  geom_line(data=tmp.tab2,aes(x=hindex,y=fit),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=fit),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0.01,0))+
  scale_y_continuous('Longevity (months)',expand=c(0.01,0),limits=c(0,15.5))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(2, 15, 0, 2, "pt"),
        legend.position = 'none')+
  ggtitle('B.')
pB2




# Load LRS data:

dat<-read.csv("./FocalFish_fitness.csv")
dat<-dat[-which(dat$FishID_nodash=='FCA1O6G8B1005' & dat$GenGrp=='F1xN'),] #drop incorrect row entry
dat2<-dat[!is.na(dat$LRS),]
dat2$cohort.fctr<-as.factor(dat2$cohort)
dat2$hindex2<-dat2$hindex^2
dat2$reprodQ<-ifelse(dat2$LRS>0,1,0)

dat2C2<-dat2[dat2$Stream=='Caigual' & dat2$cohort<=10,]
dat2T2<-dat2[dat2$Stream=='Taylor' & dat2$cohort<=10,]




### Panel C: LRS by h-index in Caigual

tmp<-c3CT[c3CT$par=='mu' & c3CT$stream=='Caigual',]
tmp.tab<-tab3[tab3$par=='mu.max' & tab3$stream=='Caigual',]
tmp.tab2<-tmp.tab
tmp.tab2$mn<-46
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$mn<-46
pC<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_point(data=dat2C2,aes(x=hindex,y=LRS),colour='purple',size=0.8,alpha=0.4)+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=mn),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(colour='purple',size=0.8)+
  geom_line(data=tmp.tab2,aes(x=hs,y=mn),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=mn),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0.01,0))+
  scale_y_continuous('Lifetime reproductive success',expand=c(0.01,0),limits=c(0,50))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(2, 15, 0, 2, "pt"))+
  ggtitle('C.')
pC


### Panel D: LRS by h-index in Taylor

tmp<-c3CT[c3CT$par=='mu' & c3CT$stream=='Taylor',]
tmp.tab<-tab3[tab3$par=='mu.max' & tab3$stream=='Taylor',]
tmp.tab2<-tmp.tab
tmp.tab2$mn<-46
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$mn<-46
pD<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_point(data=dat2T2,aes(x=hindex,y=LRS),colour='purple',size=0.8,alpha=0.4)+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=mn),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(colour='purple',size=0.8)+
  #  geom_point(data=tmp.tab,shape=18,size=2,colour='black')+ 
  geom_line(data=tmp.tab2,aes(x=hs,y=mn),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=mn),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0.01,0))+
  scale_y_continuous('Lifetime reproductive success',expand=c(0.01,0),limits=c(0,50))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(2, 15, 0, 2, "pt"))+
  ggtitle('D.')
pD




pCombO<-arrangeGrob(pA,pB2,pC,pD,nrow=2)

scl<-0.8
ggsave("./figures/Fig_2_with_obs_v1.pdf",pCombO,width=scl*8.5,height=scl*7.5)


