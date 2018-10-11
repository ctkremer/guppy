#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Code analyzing longevity of guppies

# for Fitzpatrick et al. manuscript

# developed by CT Kremer, 3/2018

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#### Load tools & functions: ####

library(dplyr)
library(mgcv)
library(ggplot2)
library(bbmle)
library(glmmTMB)
library(reshape2)
library(gridExtra)

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

#### Load data: ####
dat<-read.csv("./FocalFish_fitness.csv")
dat<-dat[-which(dat$FishID_nodash=='FCA1O6G8B1005' & dat$GenGrp=='F1xN'),]  #drop incorrect row entry
dat3<-dat[!is.na(dat$Longevity),]
dat3$cohort.fctr<-as.factor(dat3$cohort)
dat3$hindex2<-dat3$hindex^2
head(dat3)

# shift longevity variable to a [0,inf) range, to make it easier to check for inflation (zero inflation)
dat3$Longevity<-dat3$Longevity-1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Fit models: Caigual #####

dat3C<-dat3[dat3$Stream=='Caigual',]
dat3C$zero<-ifelse(dat3C$Longevity==0,0,1)
head(dat3C)

# Use negative binomial family because observations are discrete approximations (months) of an underlying, continuous variable (time). This also makes it possible to consider 0 observations. Alternatively, we could have left Longevity un-shifted and used a gamma distribution, but that makes checking for zero inflation a less obvious process.

mzinbC.0 <- glmmTMB(Longevity~hindex+(1|cohort.fctr),
                    family=nbinom2,data=dat3C)
summary(mzinbC.0)

mzinbC.1 <- glmmTMB(Longevity~hindex+I(hindex^2)+(1|cohort.fctr),
                    family=nbinom2,data=dat3C)
summary(mzinbC.1)

mzinbC.2 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat3C)
summary(mzinbC.2)

mzinbC.3 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat3C)
summary(mzinbC.3)

mzinbC.4 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~1,
                    family=nbinom2,data=dat3C)
summary(mzinbC.4)

mzinbC.5 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~Sex,
                    family=nbinom2,data=dat3C)
summary(mzinbC.5)

mzinbC.6 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2),
                    family=nbinom2,data=dat3C)
summary(mzinbC.6)

mzinbC.7 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~1,
                    family=nbinom2,data=dat3C)
summary(mzinbC.7)

mzinbC.8 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2),
                    family=nbinom2,data=dat3C)
summary(mzinbC.8)

mzinbC.9 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2)+Sex,
                    family=nbinom2,data=dat3C)
summary(mzinbC.9)

AICctab(mzinbC.0,mzinbC.1,mzinbC.2,mzinbC.3,mzinbC.4,mzinbC.5,mzinbC.6,mzinbC.7,mzinbC.8,mzinbC.9,base=T)

# best model is:
summary(mzinbC.4)

newdat1<-rbind(data.frame(hindex=seq(0,1,0.01),Sex='M',cohort.fctr='new'),data.frame(hindex=seq(0,1,0.01),Sex='F',cohort.fctr='new'))
newdat1$Sex<-as.factor(newdat1$Sex)
newdat1$cohort.fctr<-as.factor(newdat1$cohort.fctr)

pd1.mu<-predict(mzinbC.4,newdata=newdat1,allow.new.levels=T,zitype='conditional',se.fit=T)
pd1.mu<-data.frame(pd1.mu)
pd1.mu<-data.frame(newdat1,pd1.mu)
head(pd1.mu)

# and zero inflation model:
pd1.p<-predict(mzinbC.4,newdata=newdat1,allow.new.levels=T,zitype='zprob',se.fit=T)
pd1.p<-data.frame(pd1.p)
pd1.p<-data.frame(newdat2,pd1.p)
head(pd1.p)

bob<-pd1.mu %>% group_by(Sex) %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Fit models: Taylor #####

dat3T<-dat3[dat3$Stream=='Taylor',]
dat3T$zero<-ifelse(dat3T$Longevity==0,0,1)
head(dat3T)

mzinbT.0 <- glmmTMB(Longevity~hindex+(1|cohort.fctr),
                    family=nbinom2,data=dat3T)
summary(mzinbT.0)

mzinbT.1 <- glmmTMB(Longevity~hindex+I(hindex^2)+(1|cohort.fctr),
                    family=nbinom2,data=dat3T)
summary(mzinbT.1)

mzinbT.2 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat3T)
summary(mzinbT.2)

mzinbT.3 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    family=nbinom2,data=dat3T)
summary(mzinbT.3)

mzinbT.4 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~1,
                    family=nbinom2,data=dat3T)
summary(mzinbT.4)

mzinbT.5 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~Sex,
                    family=nbinom2,data=dat3T)
summary(mzinbT.5)

mzinbT.6 <- glmmTMB(Longevity~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2),
                    family=nbinom2,data=dat3T)
summary(mzinbT.6)

mzinbT.7 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~1,
                    family=nbinom2,data=dat3T)
summary(mzinbT.7)

mzinbT.8 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2),
                    family=nbinom2,data=dat3T)
summary(mzinbT.8)

mzinbT.9 <- glmmTMB(Longevity~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                    ziformula = ~hindex+I(hindex^2)+Sex,
                    family=nbinom2,data=dat3T)
summary(mzinbT.9)

# best model: mzinbT.9
AICctab(mzinbT.0,mzinbT.1,mzinbT.2,mzinbT.3,mzinbT.4,mzinbT.5,mzinbT.6,mzinbT.7,mzinbT.8,mzinbT.9,base=T)

# prep graphics:
newdat2<-rbind(data.frame(hindex=seq(0,1,0.01),Sex='M',cohort.fctr='new'),data.frame(hindex=seq(0,1,0.01),Sex='F',cohort.fctr='new'))
newdat2$Sex<-as.factor(newdat2$Sex)
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)

# mean value, conditioned on zero inflation term
pd2<-predict(mzinbT.9,newdata=newdat2,allow.new.levels=T,se.fit=T,zitype='conditional')
pd2<-data.frame(pd2)
pd2<-data.frame(newdat2,pd2)
head(pd2)

bob2.mu<-pd2 %>% group_by(Sex) %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))

ggplot(pd2,aes(x=hindex,y=fit))+
  geom_line(aes(colour=Sex))+
  scale_y_continuous(limits=c(0,5))

# and for the zero inflation sub-model:
pd2.z<-predict(mzinbT.9,newdata=newdat2,allow.new.levels=T,zitype='zprob',se.fit=T)
pd2.z<-data.frame(pd2.z)
pd2.z<-data.frame(newdat2,pd2.z)
head(pd2.z)

bob2.z<-pd2.z %>% group_by(Sex) %>% summarize(hindex=hindex[which(fit==max(fit))],fit=max(fit))

ggplot(pd2.z,aes(x=hindex,y=fit))+
  geom_line(aes(colour=Sex))+
  geom_point(data=dat3T,aes(y=zero),alpha=0.1)+
  scale_y_continuous(limits=c(0,1))

#regplot(dat3T$zero,dat3T$hindex,n=10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##### Combined predictions and maxima:  #####

# needed for adding error bars to final graphs:
bob$stream<-'Caigual'
bob$parameter<-'mu'
bob2.mu$stream<-'Taylor'
bob2.mu$parameter<-'mu'
bob2.z$stream<-'Taylor'
bob2.z$parameter<-'z'
bob3<-rbind(bob,bob2.mu,bob2.z)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###### Simulate from these model fits to determine confidence bands/intervals #####

#### First for Caigual ####

# best model
summary(mzinbC.4)

# Set up data structures and plotting:
nsims<-10000
covars2<-na.omit(dat3C[,c('Sex','hindex','cohort.fctr')])
sim.opts3<-matrix(NA,nrow=nsims,ncol=3)

# simulate from model
simC<-simulate(mzinbC.4,nsim=nsims)
sim.dat<-data.frame(covars2,simC)

# mock data frame for predicted curves
newdat2<-rbind(data.frame(hindex=seq(0,1,0.01),Sex='M',cohort.fctr='new'),data.frame(hindex=seq(0,1,0.01),Sex='F',cohort.fctr='new'))
newdat2$Sex<-as.factor(newdat2$Sex)
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)

# storage matrix
predsC.mu<-matrix(NA,nrow=nrow(newdat2),ncol=nsims)
predsC.z<-matrix(NA,nrow=nrow(newdat2),ncol=nsims)

i<-1
for(i in 1:nsims){
  print(i/nsims)
  
  # fit model to simulated data, and extract location of quadratic peak...
  mzinbC.4.sim <- glmmTMB(sim.dat[,3+i]~(hindex+I(hindex^2))*Sex+(1|cohort.fctr),
                          ziformula = ~1,family=nbinom2,data=sim.dat)
  
  # Find position of maxima
  cfs<-unlist(fixef(mzinbC.4.sim))
  
  # Females, mu
  a1 <- cfs[3]
  b1 <- cfs[2]
  c1 <- cfs[1]
  mx1 <- getmax(a1,b1,c1,0,1)

  # Males, mu
  a2 <- cfs[3]+cfs[6]
  b2 <- cfs[2]+cfs[5]
  c2 <- cfs[1]+cfs[4]
  mx2 <- getmax(a2,b2,c2,0,1)

  # stash optima:
  sim.opts3[i,]<-c(i,mx1,mx2)
  
  # save predictions for assembling simulated confidence bands:
  pd2.mu<-predict(mzinbC.4.sim,newdata=newdat2,allow.new.levels=T,zitype='conditional')
  predsC.mu[,i]<-pd2.mu
  
  # not really necessary; z parameter is independent of hindex, sex.
  pd2.z<-predict(mzinbC.4.sim,newdata=newdat2,allow.new.levels=T,zitype='zprob')
  predsC.z[,i]<-pd2.z
}
write.csv(predsC.mu,"./data/longevity_Caigual_bootstrapped_mzinbC.4_CI_bands_mu_raw.csv",row.names=F)
write.csv(predsC.z,"./data/longevity_Caigual_bootstrapped_mzinbC.4_CI_bands_z_raw.csv",row.names=F) # this output was incomplete on initial execution; re-run?
write.csv(sim.opts3,"./data/longevity_Caigual_bootstrapped_mzinbC.4_CIs_raw.csv",row.names=F)
head(sim.opts3)

# format output:
sim.opts3<-data.frame(sim.opts3)
names(sim.opts3)<-c('rep','F','M')
m1<-melt(sim.opts3,id.var='rep')
names(m1)<-c('rep','Sex','hindex')

# distribution of quadratic maxima across simulations:
gC<-ggplot(m1,aes(x=hindex))+
  geom_histogram(aes(fill=Sex))+
  scale_fill_manual(values=c('red','blue'))+
  scale_x_continuous(limits=c(0,1.01))+
  facet_wrap(~Sex)
gC

# actual bootstrapped confidence intervals:
m1 %>% group_by(Sex) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))

# A tibble: 2 x 3
#      Sex     lw.ci     up.ci
# 1      F 0.4278453 1.0000000
# 2      M 0.3018971 0.4834537

# Point estimates:
cfs<-unlist(fixef(mzinbC.4))

# Females
a1 <- cfs[3]
b1 <- cfs[2]
c1 <- cfs[1]
getmax(a1,b1,c1,0,1)

# Males
a2 <- cfs[3]+cfs[6]
b2 <- cfs[2]+cfs[5]
c2 <- cfs[1]+cfs[4]
getmax(a2,b2,c2,0,1)

# what fraction of the bootstrapped replicates is 1 or higher?
test<-function(x) length(x[x>=1])/length(x)
m1 %>% group_by(Sex) %>% summarise(p.1=test(hindex))
# 1      F 0.0924
# 2      M 0.0001

# save output
write.csv(m1,"./data/longevity_Caigual_bootstrapped_mzinbC.4_CIs.csv",row.names=F)


### Confidence bands:

# compute:
lwr.mn<-apply(predsC.mu,1,FUN=get.q,0.025)
upr.mn<-apply(predsC.mu,1,FUN=get.q,0.975)
cbandC<-data.frame(newdat2[,1:2],mn.mn=pd1.mu$fit,lwr.mn,upr.mn)
head(cbandC)

# re-arrange data:
c2C<-melt(cbandC,id.vars=c('Sex','hindex'))
c2C$type<-ifelse(grepl(pattern = 'lwr.mn',c2C$variable),'lwr',NA)
c2C$type<-ifelse(grepl(pattern = 'upr.mn',c2C$variable),'upr',c2C$type)
c2C$type<-ifelse(grepl(pattern = 'mn.mn',c2C$variable),'mn',c2C$type)
c2C$stream<-'Caigual'
c3C<-dcast(c2C[,-3],stream+Sex+hindex~type)
head(c3C)

# visualize:
ggplot(c3C,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex))+
  scale_fill_manual(values=c('red','blue'))+
  scale_colour_manual(values=c('red','blue'))+
  theme_bw()

write.csv(c3C,"./data/longevity_Caigual_bootstrapped_mzinbC.4_CI_bands.csv",row.names=F)

# In future: add processing for zi terms....

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Repeat for Taylor ####

# best model:
summary(mzinbT.9)

# Set up data structures and plotting:
nsims<-10000
covars<-na.omit(dat3T[,c('Sex','hindex','cohort.fctr')])
sim.opts4<-matrix(NA,nrow=nsims,ncol=5)

simT<-simulate(mzinbT.9,nsim=nsims)
sim.dat<-data.frame(covars,simT)

# mock data frame for predicted curves
newdat2<-rbind(data.frame(hindex=seq(0,1,0.01),Sex='M',cohort.fctr='new'),data.frame(hindex=seq(0,1,0.01),Sex='F',cohort.fctr='new'))
newdat2$Sex<-as.factor(newdat2$Sex)
newdat2$cohort.fctr<-as.factor(newdat2$cohort.fctr)

# storage matrix
predsT.mu<-matrix(NA,nrow=nrow(newdat2),ncol=nsims)
predsT.z<-matrix(NA,nrow=nrow(newdat2),ncol=nsims)

i<-1
for(i in 1:nsims){
  print(i/nsims)
  
  # fit model to simulated data, and extract location of quadratic peak...
  mzinbT.9.sim <- glmmTMB(sim.dat[,3+i]~(hindex+I(hindex^2))+Sex+(1|cohort.fctr),
                          ziformula = ~hindex+I(hindex^2)+Sex,
                          family=nbinom2,data=sim.dat)
  
  # Find position of maxima
  cfs<-unlist(fixef(mzinbT.9.sim))
  
  # Females, mu
  a1 <- cfs[3]
  b1 <- cfs[2]
  c1 <- cfs[1]
  mx1 <- getmax(a1,b1,c1,0,1)
  
  # Males, mu
  a2 <- cfs[3]
  b2 <- cfs[2]
  c2 <- cfs[1]+cfs[4]
  mx2 <- getmax(a2,b2,c2,0,1)
  
  # Females, z
  a3 <- cfs[7]
  b3 <- cfs[6]
  c3 <- cfs[5]
  mx3 <- getmax(a3,b3,c3,0,1)
  
  # Males, z
  a4 <- cfs[7]
  b4 <- cfs[6]
  c4 <- cfs[5]+cfs[8]
  mx4 <- getmax(a4,b4,c4,0,1)
  
  # stash optima:
  sim.opts4[i,]<-c(i,mx1,mx2,mx3,mx4)
  
  # save predictions for assembling simulated confidence bands:
  pd2.mu<-predict(mzinbT.9.sim,newdata=newdat2,allow.new.levels=T,zitype='conditional')
  predsT.mu[,i]<-pd2.mu
  
  pd2.z<-predict(mzinbT.9.sim,newdata=newdat2,allow.new.levels=T,zitype='zprob')
  predsT.z[,i]<-pd2.z
}
write.csv(predsT.mu,"./data/longevity_Taylor_bootstrapped_mzinbT.9_CI_bands_mu_raw.csv",row.names=F)
write.csv(predsT.z,"./data/longevity_Taylor_bootstrapped_mzinbT.9_CI_bands_zi_raw.csv",row.names=F)
write.csv(sim.opts4,"./data/longevity_Taylor_bootstrapped_mzinbT.9_CIs_raw.csv",row.names=F)
head(sim.opts4)

# format output:
sim.opts4<-data.frame(sim.opts4)
names(sim.opts4)<-c('rep','mu.F','mu.M','zi.F','zi.M')
m2<-melt(sim.opts4,id.var='rep')
m2$parameter<-ifelse(grepl(pattern = 'mu',m2$variable),'mu','z')
m2$Sex<-ifelse(grepl(pattern = 'M',m2$variable),'M','F')
names(m2)[3]<-'hindex'
names(m2)<-c('rep','Sex','hindex')
m2<-m2[,c('rep','parameter','Sex','hindex')]
head(m2)

# distribution of quadratic maxima
gT<-ggplot(m2,aes(x=hindex))+
  geom_histogram(aes(fill=Sex))+
  scale_fill_manual(values=c('red','blue'))+
  scale_x_continuous(limits=c(0,1.01))+
  #scale_y_continuous(limits=c(0,150))+
  facet_grid(parameter~Sex)
gT

# actual bootstrapped confidence intervals:
m2 %>% group_by(parameter) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))

# point estimates:
cfs<-unlist(fixef(mzinbT.9))

# Females, mu
a1 <- cfs[3]
b1 <- cfs[2]
c1 <- cfs[1]
getmax(a1,b1,c1,0,1)

# Females, z
a3 <- cfs[7]
b3 <- cfs[6]
c3 <- cfs[5]
getmax(a3,b3,c3,0,1)

# what fraction of the bootstrapped replicates is 1 or higher?
test<-function(x) length(x[x>=1])/length(x)
m2 %>% group_by(parameter) %>% summarise(p.1=test(hindex))
# A tibble: 2 x 2
#    parameter    p.1
#  1        mu 0.0105
#  2         z 0.0003

# save output
write.csv(m2,"./data/longevity_Taylor_bootstrapped_mzinbT.9_CIs.csv",row.names=F)


### Confidence bands:

# compute for conditional mean:
lwr.mn<-apply(predsT.mu,1,FUN=get.q,0.025)
upr.mn<-apply(predsT.mu,1,FUN=get.q,0.975)
cbandT<-data.frame(newdat2[,1:2],mn.mn=pd2$fit,lwr.mn,upr.mn)
head(cbandT)

# re-arrange data:
c2T<-melt(cbandT,id.vars=c('Sex','hindex'))
c2T$type<-ifelse(grepl(pattern = 'lwr.mn',c2T$variable),'lwr',NA)
c2T$type<-ifelse(grepl(pattern = 'upr.mn',c2T$variable),'upr',c2T$type)
c2T$type<-ifelse(grepl(pattern = 'mn.mn',c2T$variable),'mn',c2T$type)
c2T$stream<-'Taylor'
c3T<-dcast(c2T[,-3],stream+Sex+hindex~type)
head(c3T)

# compute for zero inflation
lwr.mn<-apply(predsT.z,1,FUN=get.q,0.025)
upr.mn<-apply(predsT.z,1,FUN=get.q,0.975)
cbandT<-data.frame(newdat2[,1:2],mn.mn=pd2.z$fit,lwr.mn,upr.mn)
head(cbandT)

c2T<-melt(cbandT,id.vars=c('Sex','hindex'))
c2T$type<-ifelse(grepl(pattern = 'lwr.mn',c2T$variable),'lwr',NA)
c2T$type<-ifelse(grepl(pattern = 'upr.mn',c2T$variable),'upr',c2T$type)
c2T$type<-ifelse(grepl(pattern = 'mn.mn',c2T$variable),'mn',c2T$type)
c2T$stream<-'Taylor'
c3T.z<-dcast(c2T[,-3],stream+Sex+hindex~type)
c3T.z$parameter<-'z'
c3T$parameter<-'mu'

# combine:
c3T<-rbind(c3T,c3T.z)

# visualize
ggplot(c3T,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex))+
  scale_fill_manual(values=c('red','blue'))+
  scale_colour_manual(values=c('red','blue'))+
  facet_wrap(~parameter,nrow=2,scales='free_y')+
  theme_bw()

write.csv(c3T,"./data/longevity_Taylor_bootstrapped_mzinbT.9_CI_bands.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Combine data on confidence intervals around quadratic maxima:
tab1<-m1 %>% group_by(Sex) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab2<-m2 %>% group_by(Sex) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab1$stream<-'Caigual'
tab2$stream<-'Taylor'

tab3<-data.frame(rbind(tab1,tab2))
tab3<-merge(tab3,bob3,by=c('Sex','stream'))

write.csv(tab3,"./data/longevity_quad_max_table.csv",row.names=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Computations for horizontal error bars ####

m1<-read.csv("./data/longevity_Taylor_bootstrapped_mzinbC.4_CIs.csv")
m1$parameter<-'mu'

m2<-read.csv("./data/longevity_Taylor_bootstrapped_mzinbT.9_CIs.csv")

# combine data on quad maximum's position:
tab1<-m1 %>% group_by(Sex,parameter) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab2<-m2 %>% group_by(Sex,parameter) %>% summarise(lw.ci=quantile(hindex,probs=0.025),up.ci=quantile(hindex,probs=0.975))
tab1$stream<-'Caigual'
tab2$stream<-'Taylor'

# combine estimated CI's from different streams:
tab3<-data.frame(rbind(tab1,tab2))

# merge with corresponding point estimates:
tab3<-merge(tab3,bob3,by=c('Sex','stream','parameter'))

write.csv(tab3,"./data/longevity_quad_max_table.csv",row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Final plotting, both streams: ####

# load data for plotting
c3C<-read.csv("./data/longevity_Caigual_bootstrapped_mzinbC.4_CI_bands.csv")
c3C$parameter<-'mu'
c3T<-read.csv("./data/longevity_Taylor_bootstrapped_mzinbT.9_CI_bands.csv")
tab3<-read.csv("./data/longevity_quad_max_table.csv")

# confidence band information:
c3CT<-rbind(c3C,c3T)

# generate final plots:
tmp<-c3CT[c3CT$stream=='Caigual',]
tmp.tab<-tab3[tab3$stream=='Caigual',]
p1L<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.3,size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Longevity',expand=c(0,0),limits=c(0,15))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('B.')
p1L

tmp<-c3CT[c3CT$stream=='Taylor' & c3CT$parameter=='mu',]
tmp.tab<-tab3[tab3$stream=='Taylor'& tab3$parameter=='mu',]
p2L<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.15,size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Longevity',expand=c(0,0),limits=c(0,5))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('A.')
p2L

tmp<-c3CT[c3CT$stream=='Taylor' & c3CT$parameter=='z',]
tmp.tab<-tab3[tab3$stream=='Taylor' & tab3$parameter=='z',]
p3L<-ggplot(tmp,aes(x=hindex,y=mn))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Sex),alpha=0.2)+
  geom_line(aes(colour=Sex),size=0.8)+
  geom_point(data=tmp.tab,aes(colour=Sex,y=fit),shape=3,size=1.5)+
  geom_errorbarh(data=tmp.tab,aes(xmin=lw.ci,xmax=up.ci,colour=Sex,y=fit),height=0.025,size=0.4)+
  scale_x_continuous('Hybrid index',expand=c(0,0))+
  scale_y_continuous('Zero-inflation probability',expand=c(0,0),limits=c(0,1))+
  facet_wrap(~stream)+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('C.')
p3L

pcombL<-arrangeGrob(p2L,p1L,p3L,nrow=2)

scl<-0.8
ggsave("./figures/Longevity_by_hindex_and_stream_newCIs_v2.pdf",pcombL,width=scl*8.5,height=scl*7.5)

