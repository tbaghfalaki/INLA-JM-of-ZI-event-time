rm(list=ls())
library(mvtnorm);library(countreg)

  n2=300;m=10
  Nobs=n2*m
  x1=rnorm(n2)
  x2=rbinom(n2,1,.5)
  D=matrix(c(1,.5,.5,1),2,2)
  solve(D)
  
  b=rmvnorm(n2,c(0,0),D)

  y=z=muy=mu1=matrix(0,n2,m)
  B=c(2,-1,1,-1)
  B1=c(1,1,-1)
  TT=seq(1:m)/m
  for(i in 1:n2){  
    for(j in 1:m){  
      muy[i,j]=B[1]+B[2]*TT[j]+B[3]*x1[i]+B[4]*x2[i]+b[i,1] 
      mu1[i,j]=B1[1]+B1[2]*j/m+B1[3]*x1[i]+ b[i,2]
      
      z[i,j]=rbinom(1,1,exp(mu1[i,j])/(1+exp(mu1[i,j])))
    }}
  theta=3
  
  for(i in 1:n2){  
    for(j in 1:m){  
      if(z[i,j]==1)(y[i,j]=0)
      if(z[i,j]==0)(y[i,j]=rztnbinom(1, mu=exp(muy[i,j]), theta=theta))
    }}
  table(z)
  r1=1
  r2=-1
  beta2=c(.5,-.5)
  mut=surt=rep(0,n2)
  r=1
  for(i in 1:n2){
    mut[i]<-(beta2[1]+beta2[2]*x1[i]+r1*b[i, 1]+r2*b[i, 2])}
  
  for (i in 1:n2) {
    surt[i] = rweibull(1,r,exp(-mut[i]/r))
  }
  status=rep(1,n2)
  quan=quantile(surt,.8)
  for (i in 1:n2) {
    
    u=runif(1)
    if(surt[i] >quan)((surt[i]=quan) & (status[i]=0))
  }
  
  for (i in 1:n2) {
    if((surt[i]>TT[2]) & (surt[i]<TT[3]))(y[i,3:m]=NA)
    if((surt[i]>TT[3]) & (surt[i]<TT[4]))(y[i,4:m]=NA)
    if((surt[i]>TT[4]) & (surt[i]<TT[5]))(y[i,5:m]=NA)
    if((surt[i]>TT[5]) & (surt[i]<TT[6]))(y[i,6:m]=NA)
    if((surt[i]>TT[6]) & (surt[i]<TT[7]))(y[i,7:m]=NA)
    if((surt[i]>TT[7]) & (surt[i]<TT[8]))(y[i,8:m]=NA)
    if((surt[i]>TT[8]) & (surt[i]<TT[9]))(y[i,9:m]=NA)
    if((surt[i]>TT[9]) & (surt[i]<TT[10]))(y[i,m]=NA)
  }
  rateobs=rep(m,n2)
  for (i in 1:n2) {
    rateobs[i]=length(y[i,][is.na(y[i,])==FALSE])
  }
  
   colnames(y)<-c("y1","y2","y3","y4","y5","y6","y7","y8","y9","y10")
  Data=cbind(y,surt,status,x1,x2,rateobs)
  head(Data)
  write.table(Data,"simulated_data_ZINB.txt")
  