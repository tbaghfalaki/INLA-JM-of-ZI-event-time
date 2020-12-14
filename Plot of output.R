Data=read.table("simulated_data_ZINB_PL.txt",header=TRUE)
attach(Data)
head(Data)

y=Data[,1:10]
y=as.matrix(y)
n2=length(y[,1])
m=length(y[1,])
###################################################
rateobs=rep(m,n2)

ID=X1=X2=Time=matrix(NA,n2,m)
for(i in 1:n2){ 
  X1[i,1:rateobs[i]]=rep(x1[i],rateobs[i])
  X2[i,1:rateobs[i]]= rep(x2[i],rateobs[i])
  Time[i,1:rateobs[i]]=c(seq(1:m)/m)[1:rateobs[i]]
  ID[i,1:rateobs[i]]=rep(i,rateobs[i])
} 

y=as.numeric(t(y))
ID=as.numeric(t(ID))
timevar=as.numeric(t(Time))
X1=as.numeric(t(X1))
X2=as.numeric(t(X2))




index=1:length(y)
indobs=index[is.na(y)==FALSE]
n1=length(indobs)
timevare =  timevar[indobs] 
X1e =  X1[indobs] 
X2e = X2[indobs] 
IDe =  ID[indobs] 
ye=y[indobs]


timevar =  timevare 
X1 =  X1e
X2 = X2e 
ID =  IDe 
y=ye



longdat1=survdat1=list()
longdat1$y=y
longdat1$TIME =timevar
longdat1$SEX =X1
longdat1$PREVOI =X2
longdat1$ID=ID




knots=quantile(timevar,prob=c(.25,.50,.75)) # B-spline
SS=bsp(timevar, k=15)$Z
#####################################
survdat1$CENSOR=status
survdat1$SURVTIME=surt
survdat1$SEX=x1

longdat=longdat1
survdat=survdat1
n1 <- length(longdat1$y)
n2 <- length(survdat1$CENSOR)

################# defining Z and its covariates ########################
Z=rep(0,length(y))
Z[y==0]=1
y[y==0]=NA
longdat$Z=Z



################## ################## ################## ###############
z.long <- c(rep(NA, n1), longdat$Z,  rep(NA, n2))
y.long <- c(longdat$y, rep(NA, n1),  rep(NA, n2))
y.surv <- inla.surv(time = c(rep(NA, 2*n1), survdat$SURVTIME), event = c(rep(NA, 2*n1), survdat$CENSOR))
Yjoint <- list(y.long, z.long, y.surv)



linear.covariate <- data.frame(mu = as.factor(c(rep(1, n1), rep(2, n1), rep(3, n2))), 
                               l.TIME = c(longdat$TIME, rep(0, n1), rep(0, n2)), 
                               l.SEX = c(longdat$SEX, rep(0, n1), rep(0, n2)), 
                               l.PREVOI = c(longdat$PREVOI, rep(0, n1), rep(0, n2)),
                               s.SEX = c(rep(0, n1), rep(0, n1), survdat$SEX))

linear.covariate$z.TIME=c(rep(0, n1), longdat$TIME, rep(0, n2))
linear.covariate$z.SEX=c(rep(0, n1), longdat$SEX, rep(0, n2))

SS0=SS*0
SS2=matrix(0,n2,length(SS[1,]))

linear.covariate$l.SS=rbind(SS,SS*0,SS2)
linear.covariate$z.SS=rbind(SS*0,SS,SS2)


ntime <- length(unique(longdat$TIME))
random.covariate <- list(U11 = c(longdat1$ID, rep(NA, n1), rep(NA, n2)),
                         U21 = c( rep(NA, n1),longdat1$ID+n2, rep(NA, n2)), 
                         U12 = c(rep(NA, n1), rep(NA, n1), 1:n2),
                         U22 = c(rep(NA, n1), rep(NA, n1), n2+(1:n2)),
                         U3 = c(rep(NA, n1), rep(NA, n1), 1:n2))


random.covariate$U11z = c(rep(NA, n1), longdat1$ID, rep(NA, n2))
random.covariate$U21z = c(rep(NA, n1),longdat1$ID+n2, rep(NA, n2))
random.covariate$U12z = c(rep(NA, n1), rep(NA, n1), 1:n2)
random.covariate$U22z = c(rep(NA, n1), rep(NA, n1), n2+(1:n2))


joint.data <- c(linear.covariate,random.covariate)
joint.data$Y <- Yjoint



formula = Y ~ mu-1 +  l.SEX +l.PREVOI+  s.SEX  + z.SEX+l.SS+z.SS+
  f(U11 , model="iid2d",n=2*n2) + 
  f(U21,  copy="U11") + f(U12, copy="U11", fixed= FALSE, param=c(1,1)) + 
  f(U22, copy="U11",fixed = FALSE, param = c(1,1))


time00 = Sys.time()
joint.inla <- inla(formula, family = c("zeroinflatednbinomial0","binomial","weibullsurv"),
                   data = joint.data, control.compute=list(dic=TRUE),control.family = list(
                     list(hyper = list(prob = list(initial = -10,
                                                   fixed = TRUE))),list(),list()))


summary(joint.inla)
#################################### Plot ##################################

round(joint.inla$summary.fixed, 4)[1,1]
AAA=round(joint.inla$summary.fixed, 4)[-c(1:7),]
tabi=length(AAA[,1])/2
tabi
betal=round(joint.inla$summary.fixed, 4)[-c(1:7),][1:tabi,1]
betaz=round(joint.inla$summary.fixed, 4)[-c(1:7),][-(1:tabi),1]
DD=cbind(c(ID,ID+100,c(1:n2)+200),joint.data$l.SS)
DD[DD[,1]==1,]
numid=as.numeric(matrix(table(ID)))
TT=seq(1:m)/m

index=1:n2
i=min(index[numid==m])
TTnew=c()
for(j in 1:numid[i]){ 
  TTnew[j]=joint.data$l.SS[1:n1,][ID==i,][j,]%*%betal
}

TTnewlambda=round(joint.inla$summary.fixed, 4)[1,1]+TTnew


TTnew=c()
for(j in 1:numid[i]){ 
  TTnew[j]=joint.data$l.SS[1:n1,][ID==i,][j,]%*%betaz
}


TTnewpi=round(joint.inla$summary.fixed, 4)[2,1]+TTnew




par(mfrow=c(1,2))
opar <- par()
par(bg="aliceblue", mar=c(4,4.5,1.7,.25))
muyt=c()
for(j in 1:m){  
  muyt[j]=sin(2*pi*TT[j])
}

plot(TT, TTnewlambda, ylim = c(-2,2), type = "l",xlab="Time",main="(a)",
     ylab=expression(paste(hat(g[1]))),lwd=2,lty=3)
points(TT,muyt,pch=20,col="darkorchid4")
######################################################
muyt=c()
for(j in 1:m){  
  muyt[j]=-cos(2*pi*TT[j])
}


plot(TT, TTnewpi, ylim = c(-2,2), type = "l",xlab="Time",,main="(b)",
     ylab=expression(paste(hat(g[2]))),lwd=2,lty=3)
points(TT,muyt,pch=20,col="darkorchid4")



