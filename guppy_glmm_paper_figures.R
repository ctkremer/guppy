
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Code for generating figures associated with analysis of guppy fitness data

# for Fitzpatrick et al. manuscript

# developed by CT Kremer, 10/2018

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Load tools & data: ####

library(dplyr)
library(mgcv)
library(ggplot2)
library(bbmle)
library(glmmTMB)
library(gridExtra)

### Load data:

# LRS
c3C<-read.csv("./data/LRS_Caigual_bootstrapped_mzinbC.1_CI_bands.csv")
c3T<-read.csv("./data/LRS_Taylor_bootstrapped_mzinbT.6_CI_bands.csv")

m1<-read.csv("./data/LRS_Caigual_bootstrapped_mzinbC.1_CIs.csv")
m2<-read.csv("./data/LRS_Taylor_bootstrapped_mzinbT.6_CIs.csv")

pd1<-read.csv("./data/LRS_pd1.csv")
pd2.mu<-read.csv("./data/LRS_pd2_mu.csv")
pd2.p<-read.csv("./data/LRS_pd2_p.csv")

# Longevity
c3C2<-read.csv("./data/longevity_Caigual_bootstrapped_mzinbC.4_CI_bands.csv")
c3C2$parameter<-'mu'
c3T2<-read.csv("./data/longevity_Taylor_bootstrapped_mzinbT.9_CI_bands.csv")
tab32<-read.csv("./data/longevity_quad_max_table.csv")

# confidence band information:
c3CT<-rbind(c3C,c3T)
c3CT2<-rbind(c3C2,c3T2)

# confidence intervals around quadratic maxima
m1$par<-"mu.max"
tab1<-m1 %>% group_by(par) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab2<-m2 %>% group_by(par) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab1$stream<-'Caigual'
tab2$stream<-'Taylor'

# Enrich this with x-axis heights...
bobC.mu<-pd1 %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))
bobT.mu<-pd2.mu %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))
bobT.p<-pd2.p %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))

tab3<-data.frame(rbind(tab1,tab2))
tab3$mn<-c(bobC.mu$fit,bobT.mu$fit,bobT.p$fit)
tab3$hs<-c(bobC.mu$hindex,bobT.mu$hindex,bobT.p$hindex)


# raw data: LRS
dat<-read.csv("./FocalFish_fitness.csv")
dat<-dat[-which(dat$FishID_nodash=='FCA1O6G8B1005' & dat$GenGrp=='F1xN'),] #drop incorrect row entry
dat2<-dat[!is.na(dat$LRS),]
dat2$cohort.fctr<-as.factor(dat2$cohort)
dat2$hindex2<-dat2$hindex^2
dat2$reprodQ<-ifelse(dat2$LRS>0,1,0)
head(dat2)

dat2C<-dat2[dat2$Stream=='Caigual' & dat2$cohort<=10,]
dat2T<-dat2[dat2$Stream=='Taylor' & dat2$cohort<=10,]
head(dat2T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Combined graphics: ####

### Panel A: Longevity in Caigual

# generate final plots:
tmp<-c3CT2[c3CT2$stream=='Caigual',]
tmp.tab<-tab32[tab32$stream=='Caigual',]
tmp.tab2<-tmp.tab
tmp.tab2$fit<-14.5
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$fit<-14.5
pA<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0,alpha=0.2,size=3)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_line(data=tmp.tab2,aes(x=hindex,y=fit,colour=Sex),linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Longevity (months)',expand=c(0,0),limits=c(0,15.5))+
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


### Panel B: Longevity in Taylor

tmp<-c3CT2[c3CT2$stream=='Taylor' & c3CT2$parameter=='mu',]
tmp.tab<-tab32[tab32$stream=='Taylor'& tab32$parameter=='mu',]

pB<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.15,size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Longevity (months)',expand=c(0,0),limits=c(0,5))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('B.')
pB

tmp<-c3CT2[c3CT2$stream=='Taylor' & c3CT2$parameter=='mu',]
tmp.tab<-tab32[tab32$stream=='Taylor'& tab32$parameter=='mu',]
tmp.tab2<-tmp.tab[tmp.tab$Sex=='M',]
tmp.tab2$fit<-14.5
tmp.tab2<-rbind(tmp.tab[tmp.tab$Sex=='M',],tmp.tab2)
tmp.tab$fit<-14.5
pB2<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=fit),colour=gray(0.7),height=0,alpha=0.2,size=3)+
  geom_line(aes(colour=Sex),size=0.8)+
#  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
#  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.35,size=0.4)+
  geom_line(data=tmp.tab2,aes(x=hindex,y=fit),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=fit),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Longevity (months)',expand=c(0,0),limits=c(0,15.5))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(2, 15, 0, 2, "pt"),
        legend.position = 'none')+
  ggtitle('B.')
pB2


### Panel C: LRS by h-index in Caigual

tmp<-c3CT[c3CT$par=='mu' & c3CT$stream=='Caigual',]
tmp.tab<-tab3[tab3$par=='mu.max' & tab3$stream=='Caigual',]
tmp.tab2<-tmp.tab
tmp.tab2$mn<-46
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$mn<-46
pC<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=mn),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(colour='purple',size=0.8)+
#  geom_point(data=tmp.tab,shape=3,size=2,colour='purple')+
#  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci),height=1,colour='purple',size=0.4)+
  geom_line(data=tmp.tab2,aes(x=hs,y=mn),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=mn),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Lifetime reproductive success',expand=c(0,0),limits=c(0,50))+
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
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=mn),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(colour='purple',size=0.8)+
#  geom_point(data=tmp.tab,shape=18,size=2,colour='black')+ 
  geom_line(data=tmp.tab2,aes(x=hs,y=mn),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=mn),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Lifetime reproductive success',expand=c(0,0),limits=c(0,50))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(2, 15, 0, 2, "pt"))+
  ggtitle('D.')
pD

grid.arrange(pA,pB2,pC,pD,nrow=2)

pComb<-arrangeGrob(pA,pB2,pC,pD,nrow=2)

scl<-0.8
ggsave("./figures/Fig_2_v1.pdf",pComb,width=scl*8.5,height=scl*7.5)





#### Supplemental figure: zero inflation in Taylor, LRS and Longevity  ####


# Longevity

tmp<-c3CT2[c3CT2$stream=='Taylor' & c3CT2$parameter=='z',]
tmp.tab<-tab32[tab32$stream=='Taylor' & tab32$parameter=='z' & tab32$Sex=='M',]
tmp.tab2<-tmp.tab
tmp.tab2$fit<-0.9
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$fit<-0.9

pSA<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=fit),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_line(data=tmp.tab2,aes(x=hindex,y=fit),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=fit),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Zero-inflation probability',expand=c(0,0),limits=c(0,1))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(10, 15, 0, 2, "pt"))+
  ggtitle('A. Longevity (months)')
pSA

# LRS:

tmp<-c3CT[c3CT$par=='z' & c3CT$stream=='Taylor',]
tmp.tab<-tab3[tab3$par=='z.max' & tab3$stream=='Taylor',]
tmp.tab2<-tmp.tab
tmp.tab2$mn<-0.9
tmp.tab2<-rbind(tmp.tab,tmp.tab2)
tmp.tab$mn<-0.9

pSB<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,y=mn),colour=gray(0.7),height=0,alpha=0.4,size=3)+
  geom_line(colour='purple',size=0.8)+
  geom_line(data=tmp.tab2,aes(x=hs,y=mn),colour='black',linetype=2,size=0.2)+
  geom_point(data=tmp.tab,aes(y=mn),colour='black',shape=18,size=2)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Zero-inflation probability',limits=c(0,1),expand=c(0,0))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(10, 15, 0, 2, "pt"))+
  ggtitle('B. Lifetime reproductive success')
pSB

grid.arrange(pSA,pSB,nrow=2)

pCombS<-arrangeGrob(pSA,pSB,nrow=2)

scl<-0.8
ggsave("./figures/Fig_Supp_v1.pdf",pCombS,width=scl*8.5/2,height=scl*7.5)


