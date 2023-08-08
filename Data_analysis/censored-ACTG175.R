rm(list=ls())
options(scipen=999)
library(nnet)
library(broom)
library(survival)
library(speff2trial)
library(tidyr)



data(ACTG175)

mydata <- ACTG175 %>% drop_na(cd496)
summary(mydata)###删失率75.6%

attach(mydata)

sd(strat)
sd(cd40)
sd(cd80)
sd(offtrt)
sd(cd496)
sd(cd820)
sd(cens)
sd(arms)
sd(days)

mul <-multinom(arms~strat+cd40+cd80+offtrt+cd496+cd820,data=mydata)
g <-tidy(mul)

print(g,n=21)

in1<-g[[3]][[1]]
a11<-g[[3]][[2]]
a12<-g[[3]][[3]]
a13<-g[[3]][[4]]
a14<-g[[3]][[5]]
a15<-g[[3]][[6]]
a16<-g[[3]][[7]]
in2<-g[[3]][[8]]
a21<-g[[3]][[9]]
a22<-g[[3]][[10]]
a23<-g[[3]][[11]]
a24<-g[[3]][[12]]
a25<-g[[3]][[13]]
a26<-g[[3]][[14]]
in3<-g[[3]][[15]]
a31<-g[[3]][[16]]
a32<-g[[3]][[17]]
a33<-g[[3]][[18]]
a34<-g[[3]][[19]]
a35<-g[[3]][[20]]
a36<-g[[3]][[21]]


##################删失时间生存函数估计
tr0 <- data.frame(dplyr::filter(mydata,arms==0))
cen0 <- 1- tr0$cens
kg0 <- survfit(Surv(days,cen0)~1,tr0)
append0 <- c(min(kg0$time),1)          
kextract0 <- data.frame(time=kg0$time,surv = kg0$surv)
dd0 <- rbind(append0,kextract0)

################################outcome回归系数估计
coefest0 <- data.frame(dplyr::filter(tr0,cens==1))
lmexpres0 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=coefest0)
opera0 <- tidy(lmexpres0)
alpha0 <- opera0[[2]][[1]]
beta10 <- opera0[[2]][[2]]
beta20 <- opera0[[2]][[3]]
beta30 <- opera0[[2]][[4]]
beta40 <- opera0[[2]][[5]]
beta50 <- opera0[[2]][[6]]
beta60 <- opera0[[2]][[7]]


tr1 <- data.frame(dplyr::filter(mydata,arms==1))
cen1 <- 1- tr1$cens
kg1 <- survfit(Surv(days,cen1)~1,tr1)
append1 <- c(min(kg1$time),1)          
kextract1 <- data.frame(time=kg1$time,surv = kg1$surv)
dd1 <- rbind(append1,kextract1)


coefest1 <- data.frame(dplyr::filter(tr1,cens==1))
lmexpres1 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=coefest1)
opera1 <- tidy(lmexpres1)
alpha1 <- opera1[[2]][[1]]
beta11 <- opera1[[2]][[2]]
beta21 <- opera1[[2]][[3]]
beta31 <- opera1[[2]][[4]]
beta41 <- opera1[[2]][[5]]
beta51 <- opera1[[2]][[6]]
beta61 <- opera1[[2]][[7]]


tr2 <- data.frame(dplyr::filter(mydata,arms==2))
cen2 <- 1- tr2$cens
kg2 <- survfit(Surv(days,cen2)~1,tr2)
append2 <- c(min(kg2$time),1)          
kextract2 <- data.frame(time=kg2$time,surv = kg2$surv)
dd2 <- rbind(append2,kextract2)


coefest2 <- data.frame(dplyr::filter(tr2,cens==1))
lmexpres2 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=coefest2)
opera2 <- tidy(lmexpres2)
alpha2 <- opera2[[2]][[1]]
beta12 <- opera2[[2]][[2]]
beta22 <- opera2[[2]][[3]]
beta32 <- opera2[[2]][[4]]
beta42 <- opera2[[2]][[5]]
beta52 <- opera2[[2]][[6]]
beta62 <- opera2[[2]][[7]]


tr3 <- data.frame(dplyr::filter(mydata,arms==3))
cen3 <- 1- tr3$cens
kg3 <- survfit(Surv(days,cen3)~1,tr3)
append3 <- c(min(kg3$time),1)          
kextract3 <- data.frame(time=kg3$time,surv = kg3$surv)
dd3 <- rbind(append3,kextract3)


coefest3 <- data.frame(dplyr::filter(tr3,cens==1))
lmexpres3 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=coefest3)
opera3 <- tidy(lmexpres3)
alpha3 <- opera3[[2]][[1]]
beta13 <- opera3[[2]][[2]]
beta23 <- opera3[[2]][[3]]
beta33 <- opera3[[2]][[4]]
beta43 <- opera3[[2]][[5]]
beta53 <- opera3[[2]][[6]]
beta63 <- opera3[[2]][[7]]

par(mfrow=c(2,2))
plot(surv ~ time,data = dd0,type = "s",main = "K-M for arms=0")
plot(surv ~ time,data = dd1,type = "s",main = "K-M for arms=1")
plot(surv ~ time,data = dd2,type = "s",main = "K-M for arms=2")
plot(surv ~ time,data = dd3,type = "s",main = "K-M for arms=3")
par(mfrow=c(1,1))

#cenal <- 1-cens
#kgal <- survfit(Surv(days,cens)~1,mydata)
#appendal <- c(min(kgal$time),1)          
#kextractal <- data.frame(time=kgal$time,surv = kgal$surv)
#ddal <- rbind(appendal,kextractal)
#plot(surv ~ time,data = ddal,type = "s",main = "K-M for all")

#####IPW
n<-1342
XC1<-rep(0,times=n)
XC2<-rep(0,times=n)
XC3<-rep(0,times=n)

phat0<-rep(0,times=n)
phat1<-rep(0,times=n)
phat2<-rep(0,times=n)
phat3<-rep(0,times=n)

ave0 <- rep(0,times=n)
ave1 <- rep(0,times=n)
ave2 <- rep(0,times=n)
ave3 <- rep(0,times=n)

r0<-rep(0,times=n)
r1<-rep(0,times=n)
r2<-rep(0,times=n)
r3<-rep(0,times=n)

chy <- rep(0,times=n)

for (i in 1:nrow(mydata)){

X <- model.matrix(~strat+cd40+cd80+offtrt+cd496+cd820,data=mydata)
C1 <- matrix(c(in1,a11,a12,a13,a14,a15,a16),ncol=1)
C2 <- matrix(c(in2,a21,a22,a23,a24,a25,a26),ncol=1)
C3 <- matrix(c(in3,a31,a32,a33,a34,a35,a36),ncol=1)

XC1[i] <- exp(X[i,]%*%C1)
XC2[i] <- exp(X[i,]%*%C2)
XC3[i] <- exp(X[i,]%*%C3)

phat0[i] <- 1/(1+XC1[i]+XC2[i]+XC3[i])
phat1[i] <- (XC1[i])*phat0[i]
phat2[i] <- (XC2[i])*phat0[i]
phat3[i] <- (XC3[i])*phat0[i]

coef0 <- c(alpha0,beta10,beta20,beta30,beta40,beta50,beta60)
coef1 <- c(alpha1,beta11,beta21,beta31,beta41,beta51,beta61)
coef2 <- c(alpha2,beta12,beta22,beta32,beta42,beta52,beta62)
coef3 <- c(alpha3,beta13,beta23,beta33,beta43,beta53,beta63)


if (arms[i]==0){
 if(cens[i]==1){
   chy[i] <- X[i,]%*%coef0
   for(k in 1:nrow(dd0)){
      if (dd0$time[k]==days[i]){
   	r0[i] <- (phat0[i]*dd0$surv[k])^(-1)
   	ave0[i] <- (chy[i]/(phat0[i]*dd0$surv[k]))
	}
   }
 }
}
else if (arms[i]==1){
 if(cens[i]==1){
   chy[i] <- X[i,]%*%coef1
   for(k in 1:nrow(dd1)){
      if (dd1$time[k]==days[i]){
   	r1[i] <- (phat1[i]*dd1$surv[k])^(-1)
   	ave1[i] <- (chy[i]/(phat1[i]*dd1$surv[k]))
	}
   }
 }
}
else if (arms[i]==2){
 if(cens[i]==1){
   chy[i] <- X[i,]%*%coef2
   for(k in 1:nrow(dd2)){
      if (dd2$time[k]==days[i]){
   	r2[i] <- (phat2[i]*dd2$surv[k])^(-1)
   	ave2[i] <- (chy[i]/(phat2[i]*dd2$surv[k]))
	}
   }
 }
}
else if (arms[i]==3){
 if(cens[i]==1){
   chy[i] <- X[i,]%*%coef3
   for(k in 1:nrow(dd3)){
      if (dd3$time[k]==days[i]){
   	r3[i] <- (phat3[i]*dd3$surv[k])^(-1)
   	ave3[i] <- (chy[i]/(phat3[i]*dd3$surv[k]))
	}
   }
 }
}


}

est0 <- sum(ave0)/sum(r0)
est1 <- sum(ave1)/sum(r1)
est2 <- sum(ave2)/sum(r2)
est3 <- sum(ave3)/sum(r3)

ATE10 <- est1-est0
ATE10

ATE20 <- est2-est0
ATE20

ATE30 <- est3-est0
ATE30

ATE21 <- est2-est1
ATE21

ATE31 <- est3-est1
ATE31

ATE32 <- est3-est2
ATE32



