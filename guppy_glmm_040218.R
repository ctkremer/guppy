#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Code analyzing lifetime reproductive success of guppies

# for Fitzpatrick et al. manuscript

# developed by CT Kremer, 3/2018

# data from Caigual and Taylor will be analyzed separately

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Load tools & functions: ####

library(dplyr)
library(mgcv)
library(ggplot2)
library(bbmle)
library(glmmTMB)

# function for computing the maxima of a quadratic function on a finite interval (x1,x2):
getmax<-function(a,b,c,x1,x2){
  mx <- -b/(2*a)
  y1<-a*x1*x1+b*x1+c
  y2<-a*x2*x2+b*x2+c
  
  if(a<0 & mx >= x1 & mx <=x2){  # if fnx is concave down w/internal max
    res<-mx
  }else{
    res<-c(x1,x2)[which(c(y1,y2)==max(c(y1,y2)))]
  }
  return(unlist(res))
}

# helper function for calculating confidence intervals:
get.q<-function(x,p){
  as.vector(quantile(x,probs=c(p)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Load & format data: ####

dat<-read.csv("/Users/colin/Research/Active/guppy/FocalFish_fitness.csv")
dat<-dat[-which(dat$FishID_nodash=='FCA1O6G8B1005' & dat$GenGrp=='F1xN'),] #drop incorrect row entry
dat2<-dat[!is.na(dat$LRS),]
dat2$cohort.fctr<-as.factor(dat2$cohort)
dat2$hindex2<-dat2$hindex^2
dat2$reprodQ<-ifelse(dat2$LRS>0,1,0)
head(dat2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Caigual model fitting ####
dat2C<-dat2[dat2$Stream=='Caigual',]
dat2C<-dat2C[dat2C$cohort<=10,]  # omit last 3 cohorts, as their reproductive success is underestimated
head(dat2C)


# interecept model
mzinbC.0 <- glmmTMB(LRS~hindex+(1|cohort.fctr),
                    family=nbinom2,data=dat2C)
summary(mzinbC.0)

# current best model:
mzinbC.1 <- glmmTMB(LRS~hindex+I(hindex^2)+(1|cohort.fctr),
                    family=nbinom2,data=dat2C)
summary(mzinbC.1)

mzinbC.2 <- glmmTMB(LRS~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat2C)
summary(mzinbC.2)

mzinbC.3 <- glmmTMB(LRS~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat2C)
summary(mzinbC.3)

mzinbC.4 <- glmmTMB(LRS~hindex+I(hindex^2)+(1|cohort.fctr),
                    ziformula = ~1,
                    family=nbinom2,data=dat2C)
summary(mzinbC.4)

mzinbC.5 <- glmmTMB(LRS~hindex+I(hindex^2)+(1|cohort.fctr),
                    ziformula = ~hindex,
                    family=nbinom2,data=dat2C)
summary(mzinbC.5)

# Model comparison
# mxinbC.1 has lowest AICc
AICctab(mzinbC.0,mzinbC.1,mzinbC.2,mzinbC.3,mzinbC.4,mzinbC.5,base=T)

# Visualize fit of best model
newdat1<-data.frame(hindex=seq(0,1,0.01),cohort.fctr='new')
newdat1$cohort.fctr<-as.factor(newdat1$cohort.fctr)

pd1<-predict(mzinbC.1,newdata=newdat1,allow.new.levels=T,se.fit=T)
pd1<-data.frame(pd1)
pd1<-data.frame(newdat1,pd1)

# Save these predictions/error bands for later:
write.csv(pd1,"/Users/colin/Research/Active/guppy/data/LRS_pd1.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Taylor model fitting ####
dat2T<-dat2[dat2$Stream=='Taylor',]
dat2T<-dat2T[dat2T$cohort<=10,]
head(dat2T)

# fit models
mzinbT.0 <- glmmTMB(LRS~hindex+(1|cohort.fctr),
                    family=nbinom2,data=dat2T)
summary(mzinbT.0)

mzinbT.1 <- glmmTMB(LRS~(hindex+I(hindex^2))+(1|cohort.fctr),
                    family=nbinom2,data=dat2T)
summary(mzinbT.1)

mzinbT.2 <- glmmTMB(LRS~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat2T)
summary(mzinbT.2)

mzinbT.3 <- glmmTMB(LRS~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat2T)
summary(mzinbT.3)

mzinbT.4 <- glmmTMB(LRS~(hindex+I(hindex^2))+(1|cohort.fctr),
                    ziformula=~1,
                    family=nbinom2,data=dat2T)
summary(mzinbT.4)

mzinbT.5 <- glmmTMB(LRS~(hindex+I(hindex^2))+(1|cohort.fctr),
                    ziformula=~hindex,
                    family=nbinom2,data=dat2T)
summary(mzinbT.5)

# current best-fit model
mzinbT.6 <- glmmTMB(LRS~(hindex+I(hindex^2))+(1|cohort.fctr),
                    ziformula=~(hindex+I(hindex^2)),
                    family=nbinom2,data=dat2T)
summary(mzinbT.6)

mzinbT.7 <- glmmTMB(LRS~(hindex+I(hindex^2))+(1|cohort.fctr),
                    ziformula=~(hindex+I(hindex^2))+Sex,
                    family=nbinom2,data=dat2T)
summary(mzinbT.7)

# Model comparison
AICctab(mzinbT.0,mzinbT.1,mzinbT.2,mzinbT.3,mzinbT.4,mzinbT.5,mzinbT.6,mzinbT.7,base=T)

# Generate predictions and plots:
newdat2<-data.frame(hindex=seq(0,1,0.01),cohort.fctr='new')
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)

# For the conditional mean (conditioned on zero inflation)
pd2.mu<-predict(mzinbT.6,newdata=newdat2,allow.new.levels=T,zitype='conditional',se.fit=T)
pd2.mu<-data.frame(pd2.mu)
pd2.mu<-data.frame(newdat2,pd2.mu)
head(pd2.mu)

# and zero inflation model:
pd2.p<-predict(mzinbT.6,newdata=newdat2,allow.new.levels=T,zitype='zprob',se.fit=T)
pd2.p<-data.frame(pd2.p)
pd2.p<-data.frame(newdat2,pd2.p)
head(pd2.p)

# save predictions
write.csv(pd2.mu,"/Users/colin/Research/Active/guppy/data/LRS_pd2_mu.csv",row.names=F)
write.csv(pd2.p,"/Users/colin/Research/Active/guppy/data/LRS_pd2_p.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###### Set up Monte Carlo simulaion from the best model fits, to assess uncertainty #####

# Best Caigual model: 
summary(mzinbC.1)

# Best Taylor model:
summary(mzinbT.6)

# neither model suggested significant differences by sex...

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Focus on Caigual ####

summary(mzinbC.1)

# Find position of maxima
cfs<-unlist(fixef(mzinbC.1))
getmax(cfs[3],cfs[2],cfs[1],0,1)

# Set up data structures and plotting:
nsims<-10000
covars2<-na.omit(dat2C[,c('hindex','cohort.fctr')])
sim.opts3<-matrix(NA,nrow=nsims,ncol=2)

simC<-simulate(mzinbC.1,nsim=nsims)
sim.dat<-data.frame(covars2,simC)

hs<-seq(0,1,0.01)
newdat2<-data.frame(hindex=hs,cohort.fctr='new')
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)
preds.mn<-matrix(NA,nrow=length(hs),ncol=nsims)

for(i in 1:nsims){
  print(i/nsims)
  
  # fit model to simulated data, and extract location of quadratic peak...
  mzinbC.1.sim <- glmmTMB(sim.dat[,2+i]~(hindex+I(hindex^2))+
                            (1|cohort.fctr),family=nbinom2,data=sim.dat)
  
  # Find position of maxima
  cfs<-unlist(fixef(mzinbC.1.sim))
  
  # Male or female:
  a1 <- cfs[3]
  b1 <- cfs[2]
  c1 <- cfs[1]
  mx1 <- getmax(a1,b1,c1,0,1)
  
  # stash optima:
  sim.opts3[i,]<-c(i,mx1)
  
  # save predictions for assembling simulated confidence bands:
  pd.mn<-predict(mzinbC.1.sim,newdata=newdat2,allow.new.levels=T)
  preds.mn[,i]<-pd.mn
}
write.csv(preds.mn,"/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CI_bands_mu_raw.csv",row.names=F)
write.csv(sim.opts3,"/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CIs_raw.csv",row.names=F)
head(sim.opts3)

# format output:
sim.opts3<-data.frame(sim.opts3)
names(sim.opts3)<-c('rep','MF')
#m1<-melt(sim.opts3,id.var='rep')
m1<-sim.opts3
names(m1)<-c('rep','hindex')
head(m1)

# actual bootstrapped confidence intervals:
m1 %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))

#        lw.ci up.ci
#  1 0.6629158     1

# what fraction of the bootstrapped replicates is 1 or higher?
test<-function(x) length(x[x>=1])/length(x)
m1 %>% summarise(p.1=test(hindex))
#     p.1
# 1 0.094

# save output
write.csv(m1,"/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CIs.csv",row.names=F)

### Compute confidence bands:
lwr.mn<-apply(preds.mn,1,FUN=get.q,0.025)
upr.mn<-apply(preds.mn,1,FUN=get.q,0.975)
cbandC<-data.frame(hs,lwr.mn,upr.mn)
cbandC$mn.mn<-pd1$fit
head(cbandC)

c2C<-melt(cbandC,id.vars='hs')
c2C$par<-'mu'
c2C$type<-ifelse(grepl(pattern = 'lwr.mn',c2C$variable),'lwr',NA)
c2C$type<-ifelse(grepl(pattern = 'upr.mn',c2C$variable),'upr',c2C$type)
c2C$type<-ifelse(grepl(pattern = 'mn.mn',c2C$variable),'mn',c2C$type)
c2C$stream<-'Caigual'
c3C<-dcast(c2C[,-2],stream+hs+par~type)
head(c3C)

ggplot(c3C,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_line(colour='purple')+
  theme_bw()

write.csv(c3C,"/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CI_bands.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#### Repeat for Taylor ####

# Best model:
summary(mzinbT.6)

# Find position of maxima
cfs<-unlist(fixef(mzinbT.6))
getmax(cfs[3],cfs[2],cfs[1],0,1)
getmax(cfs[6],cfs[5],cfs[4],0,1)


# Set up data structures and plotting:
nsims<-10000
covars<-na.omit(dat2T[,c('hindex','cohort.fctr')])
sim.opts4<-matrix(NA,nrow=nsims,ncol=3)

simT<-simulate(mzinbT.6,nsim=nsims)
sim.dat<-data.frame(covars,simT)

hs<-seq(0,1,0.01)
newdat2<-data.frame(hindex=hs,cohort.fctr='new')
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)
preds.mu<-matrix(NA,nrow=length(hs),ncol=nsims)
preds.z<-matrix(NA,nrow=length(hs),ncol=nsims)

for(i in 1:nsims){
  print(i/nsims)
  
  # fit model to simulated data, and extract location of quadratic peak...
  mzinbT.6.sim <- glmmTMB(sim.dat[,2+i]~(hindex+I(hindex^2))+(1|cohort.fctr),
                          ziformula = ~hindex+I(hindex^2),family=nbinom2,data=sim.dat)
  
  # Find position of maxima
  cfs<-unlist(fixef(mzinbT.6.sim))
  
  # LRS (M/F)
  a1 <- cfs[3]
  b1 <- cfs[2]
  c1 <- cfs[1]
  mx1 <- getmax(a1,b1,c1,0,1)
  
  # LRS zi (M/F)
  a2 <- cfs[6]
  b2 <- cfs[5]
  c2 <- cfs[4]
  mx2 <- getmax(a2,b2,c2,0,1)
  
  # stash optima:
  sim.opts4[i,]<-c(i,mx1,mx2)
  
  # save predictions for assembling simulated confidence bands:
  pd2.mu<-predict(mzinbT.6.sim,newdata=newdat2,allow.new.levels=T,zitype='conditional')
  preds.mu[,i]<-pd2.mu
  
  pd2.z<-predict(mzinbT.6.sim,newdata=newdat2,allow.new.levels=T,zitype='zprob')
  preds.z[,i]<-pd2.z
}
write.csv(preds.mu,"/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CI_bands_mu_raw.csv",row.names=F)
write.csv(preds.z,"/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CI_bands_z_raw.csv",row.names=F)
write.csv(sim.opts4,"/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CIs_raw.csv",row.names=F)
head(sim.opts4)

# format output:
sim.opts4<-data.frame(sim.opts4)
names(sim.opts4)<-c('rep','mu.max','z.max')
m2<-melt(sim.opts4,id.var='rep')
names(m2)<-c('rep','par','hindex')

ggplot(m2,aes(x=hindex))+
  geom_histogram(aes(fill=par))+
  scale_fill_manual(values=c('red','blue'))+
  scale_x_continuous(limits=c(0,1.01))+
  facet_wrap(~par)

# actual bootstrapped confidence intervals:
m2 %>% group_by(par) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))

# A tibble: 2 x 3
#       par     lw.ci     up.ci
#  1 mu.max 0.5854822 0.9456893
#  2  z.max 0.0000000 0.3093356

# what fraction of the bootstrapped replicates is 1 or higher?
test<-function(x) length(x[x>=1])/length(x)
m2 %>% group_by(par) %>% summarise(p.1=test(hindex))
# 1 mu.max 0.0173

test2<-function(x) length(x[x==0])/length(x)
m2 %>% group_by(par) %>% summarise(p.1=test2(hindex))
# 2  z.max 0.3772

# save output
write.csv(m2,"/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CIs.csv",row.names=F)

### Compute confidence bands:
lwr.mu<-apply(preds.mu,1,FUN=get.q,0.025)
upr.mu<-apply(preds.mu,1,FUN=get.q,0.975)
lwr.z<-apply(preds.z,1,FUN=get.q,0.025)
upr.z<-apply(preds.z,1,FUN=get.q,0.975)
cbandT<-data.frame(hs,lwr.mu,upr.mu,lwr.z,upr.z)
cbandT$mn.mu<-pd2.mu$fit
cbandT$mn.z<-pd2.p$fit
head(cbandT)

c2T<-melt(cbandT,id.vars='hs')
c2T$par<-ifelse(grepl(pattern = 'mu',c2T$variable),'mu','z')
c2T$type<-ifelse(grepl(pattern = 'lwr',c2T$variable),'lwr',NA)
c2T$type<-ifelse(grepl(pattern = 'upr',c2T$variable),'upr',c2T$type)
c2T$type<-ifelse(grepl(pattern = 'mn',c2T$variable),'mn',c2T$type)
c2T$stream<-'Taylor'
c3T<-dcast(c2T[,-2],stream+hs+par~type)
head(c3T)

ggplot(c3T,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_line(colour='purple')+
  facet_wrap(~par,scales='free_y')+
  theme_bw()

write.csv(c3T,"/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CI_bands.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Combined graphics:

# quick load data:
c3C<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CI_bands.csv")

# file missing:
#c3T<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CI_bands.csv")
m1<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_Caigual_bootstrapped_mzinbC.1_CIs.csv")

# file missing:
#m2<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_Taylor_bootstrapped_mzinbT.6_CIs.csv")
pd1<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_pd1.csv")
pd2.mu<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_pd2_mu.csv")
pd2.p<-read.csv("/Users/colin/Research/Active/guppy/data/LRS_pd2_p.csv")

# confidence band information:
c3CT<-rbind(c3C,c3T)

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

tmp<-c3CT[c3CT$par=='mu' & c3CT$stream=='Caigual',]
tmp.tab<-tab3[tab3$par=='mu.max' & tab3$stream=='Caigual',]
p1C<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_line(colour='purple',size=0.8)+
  geom_point(data=tmp.tab,shape=3,size=2,colour='purple')+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci),height=1,colour='purple',size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Reproductive success',expand=c(0,0),limits=c(0,50))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('B.')

tmp<-c3CT[c3CT$par=='mu' & c3CT$stream=='Taylor',]
tmp.tab<-tab3[tab3$par=='mu.max' & tab3$stream=='Taylor',]
p1T<-ggplot(tmp,aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_line(colour='purple',size=0.8)+
  geom_point(data=tmp.tab,shape=3,size=2,colour='purple')+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci),height=1,colour='purple',size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Reproductive success',expand=c(0,0),limits=c(0,50))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('A.')

p2<-ggplot(c3CT[c3CT$par=='z',],aes(x=hs,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr),fill='purple',alpha=0.2)+
  geom_line(colour='purple',size=0.8)+
  geom_point(data=tab3[tab3$par=='z.max',],shape=3,size=2,colour='purple')+
  geom_errorbarh(data=tab3[tab3$par=='z.max',],aes(xmin=lw.ci,xmax=up.ci),height=0.05,colour='purple',size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Zero-inflation probability',limits=c(0,1),expand=c(0,0))+
  facet_wrap(~stream)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('C.')

grid.arrange(p1T,p1C,p2,nrow=2)

pcomb<-arrangeGrob(p1T,p1C,p2,nrow=2)

scl<-0.8
ggsave("/Users/colin/Research/Active/guppy/figures/LRS_by_hindex_and_stream_newCIs.pdf",pcomb,width=scl*7.75,height=scl*7.5)

