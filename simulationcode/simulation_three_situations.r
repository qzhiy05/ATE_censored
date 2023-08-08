
rm(list=ls()); options(scipen=999)##科学计数法
library(nnet)
library(broom)
library(tidyr)
library(survival)

# Simulate
m <- 1000

simipw <- function(x1,x2,x3){

n = 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n,1,0.3)

p1 <- exp(0.3+x1+x2+x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))
p2 <- exp(0.8+2*x1+2*x2+2*x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))
p3 <- exp(0.1+3*x1+3*x2+3*x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))


A <- 1:n;
  for (i in 1:n){ 
    p <- c(p1[i],p2[i],p3[i])
    A[i]<- sample(1:3, 1, replace=T, prob=p)
  }


wheat <- data.frame(A,x1,x2,x3)




mod.fitr<-multinom(A ~ x1 + x2 + x3, data = wheat)
#summary(mod.fitr)


g<-tidy(mod.fitr)

a1<-g[[3]][[1]]
a2<-g[[3]][[2]]
a3<-g[[3]][[3]]
a4<-g[[3]][[4]]
a5<-g[[3]][[5]]
a6<-g[[3]][[6]]
a7<-g[[3]][[7]]
a8<-g[[3]][[8]]


XC1<-rep(0,times=n)
XC2<-rep(0,times=n)
phat1<-rep(0,times=n)
phat2<-rep(0,times=n)
phat3<-rep(0,times=n)

err <- rnorm(n)

y <- 1:n;

ave1 <- rep(0,times=n)
ave2 <- rep(0,times=n)
ave3 <- rep(0,times=n)
r1<-rep(0,times=n)
r2<-rep(0,times=n)
r3<-rep(0,times=n)


c1ave1 <- rep(0,times=n)
c1ave2 <- rep(0,times=n)
c1ave3 <- rep(0,times=n)
c1r1<-rep(0,times=n)
c1r2<-rep(0,times=n)
c1r3<-rep(0,times=n)

c2ave1 <- rep(0,times=n)
c2ave2 <- rep(0,times=n)
c2ave3 <- rep(0,times=n)
c2r1<-rep(0,times=n)
c2r2<-rep(0,times=n)
c2r3<-rep(0,times=n)



###################################################

X <- model.matrix(~x1+x2+x3)
C1 <- matrix(c(a1,a2,a3,a4),ncol=1)
C2 <- matrix(c(a5,a6,a7,a8),ncol=1)


I<-matrix(NA,nrow=3,ncol=n)

for (i in 1:n){


XC1[i] <- exp(X[i,]%*%C1)
XC2[i] <- exp(X[i,]%*%C2)


phat1[i] <- 1/(1+XC1[i]+XC2[i])
phat2[i] <- (XC1[i])*phat1[i]
phat3[i] <- (XC2[i])*phat1[i]



I[,i]<-rmultinom(1,1,c(phat1[i],phat2[i],phat3[i]))

alpha <- c(0,1,2)
beta1 <- c(1,0.5,1)
beta2 <- c(2,1,2)
beta3 <- c(3,2,1)

y[i] <- t(I[,i])%*%alpha+t(I[,i])%*%beta1*x1[i]+t(I[,i])%*%beta2*x2[i]+t(I[,i])%*%beta3*x3[i]+err[i]


}


##################################NO CENSORED
comalda <- data.frame(y,t(I))
for (i in 1:nrow(comalda)){
if (comalda$X1[i]==1){
   r1[i] <- phat1[i]^(-1)
   ave1[i] <- (y[i]/phat1[i])

}
else if (comalda$X2[i]==1){
   r2[i] <- phat2[i]^(-1)
   ave2[i] <- (y[i]/phat2[i])

}
else if (comalda$X3[i]==1){
   r3[i] <- phat3[i]^(-1)
   ave3[i] <- (y[i]/phat3[i])

}

}

est1 <- sum(ave1)/sum(r1)
est2 <- sum(ave2)/sum(r2)
est3 <- sum(ave3)/sum(r3)

ATE21 <- est2-est1
ATE21

ATE32 <- est3-est2
ATE32

ATE31 <- est3-est1
ATE31


censortime1 <- runif(n,0,11)
delta1 <- numeric(n)
delta1[censortime1 > y] = 1
T = pmin(censortime1,y) 
#expodata = data.frame(T,delta1) 
#summary(expodata)###删失率


censortime2 <- runif(n,0,5)
delta2 <- numeric(n)
delta2[censortime2 > y] = 1
TT = pmin(censortime2,y)
#expodataa = data.frame(TT,delta2) 
#summary(expodataa)###删失率



############# 估计删失时间的生存函数


####约20%删失的情况
comalda1 <- data.frame(y,delta1,censortime1,t(I),T)   ####I转置后，三个列向量名分别为X1,X2,X3

u1 <- data.frame(dplyr::filter(comalda1,X1==1))
cen1 <- 1-u1$delta1
k1 <- survfit(Surv(y,cen1)~1,u1)
append1 <- c(min(k1$time),1)          
kextract1 <- data.frame(time=k1$time,surv = k1$surv)
dd1 <- rbind(append1,kextract1)

u2 <- data.frame(dplyr::filter(comalda1,X2==1))
cen2 <- 1-u2$delta1
k2 <- survfit(Surv(y,cen2)~1,u2)
append2 <- c(min(k2$time),1)          
kextract2 <- data.frame(time=k2$time,surv = k2$surv)
dd2 <- rbind(append2,kextract2)

u3 <- data.frame(dplyr::filter(comalda1,X3==1))
cen3 <- 1-u3$delta1
k3 <- survfit(Surv(y,cen3)~1,u3)
append3 <- c(min(k3$time),1)          
kextract3 <- data.frame(time=k3$time,surv = k3$surv)
dd3 <- rbind(append3,kextract3)

####约40%删失的情况
comalda2 <- data.frame(y,delta2,censortime2,t(I),TT)   ####I转置后，三个列向量名分别为X1,X2,X3

uu1 <- data.frame(dplyr::filter(comalda2,X1==1))
ccen1 <- 1-uu1$delta2
kk1 <- survfit(Surv(y,ccen1)~1,uu1)
aappend1 <- c(min(kk1$time),1)          
kkextract1 <- data.frame(time=kk1$time,surv = kk1$surv)
sdd1 <- rbind(aappend1,kkextract1)

uu2 <- data.frame(dplyr::filter(comalda2,X2==1))
ccen2 <- 1-uu2$delta2
kk2 <- survfit(Surv(y,ccen2)~1,uu2)
aappend2 <- c(min(kk2$time),1)          
kkextract2 <- data.frame(time=kk2$time,surv = kk2$surv)
sdd2 <- rbind(aappend2,kkextract2)

uu3 <- data.frame(dplyr::filter(comalda2,X3==1))
ccen3 <- 1-uu3$delta2
kk3 <- survfit(Surv(y,ccen3)~1,uu3)
aappend3 <- c(min(kk3$time),1)          
kkextract3 <- data.frame(time=kk3$time,surv = kk3$surv)
sdd3 <- rbind(aappend3,kkextract3)

################################20%CENCORED
for( i in 1:nrow(comalda1)){

if (comalda1$X1[i]==1){
if(delta1[i]==1){
for(k in 1:nrow(dd1)){
      if (dd1$time[k]==T[i]){

   	c1r1[i] <- (phat1[i]*dd1$surv[k])^(-1)
   	c1ave1[i] <- (T[i]/(phat1[i]*dd1$surv[k]))
	
	}
   }
 }
}
else if (comalda1$X2[i]==1){
if(delta1[i]==1){
for(k in 1:nrow(dd2)){
      if (dd2$time[k]==T[i]){

   	c1r2[i] <- (phat2[i]*dd2$surv[k])^(-1)
   	c1ave2[i] <- (T[i]/(phat2[i]*dd2$surv[k]))
	
	}
   }
 }
}
else if (comalda1$X3[i]==1){
if(delta1[i]==1){
for(k in 1:nrow(dd3)){
      if (dd3$time[k]==T[i]){

   	c1r3[i] <- (phat3[i]*dd3$surv[k])^(-1)
   	c1ave3[i] <- (T[i]/(phat3[i]*dd3$surv[k]))
	
	}
   }
 }
}


}

fest1 <- sum(c1ave1)/sum(c1r1)
fest2 <- sum(c1ave2)/sum(c1r2)
fest3 <- sum(c1ave3)/sum(c1r3)


fATE21 <- fest2-fest1
fATE21

fATE31 <- fest3-fest1
fATE31

fATE32 <- fest3-fest2
fATE32

##############################40%CENCORED
for( i in 1:nrow(comalda2)){

if (comalda2$X1[i]==1){
if(delta2[i]==1){
for(k in 1:nrow(sdd1)){
      if (sdd1$time[k]==TT[i]){

   	c2r1[i] <- (phat1[i]*sdd1$surv[k])^(-1)
   	c2ave1[i] <- (TT[i]/(phat1[i]*sdd1$surv[k]))
	
	}
   }
 }
}
else if (comalda2$X2[i]==1){
if(delta2[i]==1){
for(k in 1:nrow(sdd2)){
      if (sdd2$time[k]==TT[i]){

   	c2r2[i] <- (phat2[i]*sdd2$surv[k])^(-1)
   	c2ave2[i] <- (TT[i]/(phat2[i]*sdd2$surv[k]))
	
	}
   }
 }
}
else if (comalda2$X3[i]==1){
if(delta2[i]==1){
for(k in 1:nrow(sdd3)){
      if (sdd3$time[k]==TT[i]){

   	c2r3[i] <- (phat3[i]*sdd3$surv[k])^(-1)
   	c2ave3[i] <- (TT[i]/(phat3[i]*sdd3$surv[k]))
	
	}
   }
 }
}


}

sest1 <- sum(c2ave1)/sum(c2r1)
sest2 <- sum(c2ave2)/sum(c2r2)
sest3 <- sum(c2ave3)/sum(c2r3)


sATE21 <- sest2-sest1
sATE21

sATE31 <- sest3-sest1
sATE31

sATE32 <- sest3-sest2
sATE32


ARES<-matrix(c(ATE21,ATE31,ATE32,fATE21,fATE31,fATE32,sATE21,sATE31,sATE32),nrow=1)

return(ARES)
}

resultp <- matrix(NA, nrow=m, ncol=9) #creat a matrix to hold outcome
for(i in 1:m) resultp[i,] <- simipw(x1,x2,x3)


mean(resultp[,1]) ##ATE21--IPW
var(resultp[,1]) ##ATE21--IPW

mean(resultp[,2]) ##ATE31--IPW
var(resultp[,2]) ##ATE31--IPW

mean(resultp[,3]) ##ATE32--IPW
var(resultp[,3]) ##ATE32--IPW

mean(resultp[,4]) ##ATE21--C1IPW
var(resultp[,4]) ##ATE21--C1IPW

mean(resultp[,5]) ##ATE31--C1IPW
var(resultp[,5]) ##ATE31--C1IPW

mean(resultp[,6]) ##ATE32--C1IPW
var(resultp[,6]) ##ATE32--C1IPW

mean(resultp[,7]) ##ATE21--C2IPW
var(resultp[,7]) ##ATE21--C2IPW

mean(resultp[,8]) ##ATE31--C2IPW
var(resultp[,8]) ##ATE31--C2IPW

mean(resultp[,9]) ##ATE32--C2IPW
var(resultp[,9]) ##ATE32--C2IPW


bia21 <- mean(resultp[,1]-0.7) ## ATE21--bias
bia21

bia31 <- mean(resultp[,2]-1.4) ## ATE31--bias
bia31

bia32 <- mean(resultp[,3]-0.7) ## ATE32--bias
bia32

bia1c21 <- mean(resultp[,4]-0.7) ## ATE21--biasC1
bia1c21

bia1c31 <- mean(resultp[,5]-1.4) ## ATE31--biasC1
bia1c31

bia1c32 <- mean(resultp[,6]-0.7) ## ATE32--biasC1
bia1c32

bia2c21 <- mean(resultp[,7]-0.7) ## ATE21--biasC2
bia2c21

bia2c31 <- mean(resultp[,8]-1.4) ## ATE31--biasC2
bia2c31

bia2c32 <- mean(resultp[,9]-0.7) ## ATE32--biasC2
bia2c32


