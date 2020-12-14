  rm(list=ls())
  library(R2OpenBUGS)
  Data=read.table("simulated_data_ZINB.txt",header=TRUE)
  attach(Data)
  head(Data)
  
  y=Data[,1:10]
  y=as.matrix(y)
  n=length(y[,1])
  m=length(y[1,])
  
  betamu1=c(0,0,0,0)
  betamu2=c(0,0,0)
  betamu3=c(0,0)
  Sigma3=structure(.Data=c(0.01,0,
                           0,0.01),.Dim=c(2,2))
  
  Sigma2=structure(.Data=c(0.01,0,0,
                           0,0.01,0,
                           0,0,0.01),.Dim=c(3,3))
  
  Sigma1=structure(.Data=c(0.01,0,0,0,
                           0,0.01,0,0,
                           0,0,0.01,0,
                           0,0,0,0.01),.Dim=c(4,4))
  
  
  Sigma3=diag(2)
  Sigma2=diag(3)
  Sigma1=diag(4)
  U0=c(0,0)
  R=structure(.Data=c(1, 0, 0,1), .Dim=c(2,2))
  
  #****************************************
  sink("model.file.txt")        
  
  
  cat("model{
  K<-1000
  for (i in 1:n) {
    for (j in 1:rateobs[i]) {
      zeros[i,j]<-0
      zeros[i,j]~dpois(phi[i,j])
      phi[i,j]<-  - ll[i,j]+K				
      ll[i,j]<-z[i,j]*log(pi[i,j]) +
 (1-z[i,j])*(
log(1-pi[i,j])+loggam(r+y[i,j])-loggam(r)-loggam(y[i,j]+1)+r*log(r/(r+lambda[i, j]))+y[i,j]*log(lambda[i, j]/(lambda[i, j]+r))-log(1-pow(r/(r+lambda[i, j]),r))   )
               
      log(lambda[i, j])<-beta1[1]+beta1[2]*t[j]+beta1[3]*x1[i]+beta1[4]*x2[i]+U[i,1]
      logit(pi[i, j])<-beta2[1]+beta2[2]*t[j]+beta2[3]*x1[i]+U[i,2]
    }
    surt[i] ~ dweib(p,mut[i]) I(surt.cen[i],)      

    log(mut[i])<-beta3[1]+beta3[2]*x1[i]+r1*U[i, 1]+r2*U[i, 2]
    U[i,1:2] ~ dmnorm(U0[],tau[,])   
  }  
  r~ dgamma(1,1) 
  p ~ dgamma(1,1) 
  sigma[1:2,1:2]<-inverse(tau[,])
  
  #priors
  tau[1:2,1:2] ~ dwish(R[,], 2)
  beta1[1:4]~dmnorm(betamu1[],Sigma1[,])
  beta2[1:3]~dmnorm(betamu2[],Sigma2[,])
  beta3[1:2]~dmnorm(betamu3[],Sigma3[,])
  r1~dnorm(0, 0.0001)
  r2~dnorm(0, 0.0001)
}

", fill = TRUE)
  
  sink()
  #****************************************
  t=seq(1:m)/m

  delta=1-status
  surt.cen=rep(0,n)
  for(i in 1:n){
    if(status[i]==0)((surt.cen[i]=surt[i]) & (surt[i]=NA))
    if(status[i]==1)((surt.cen[i]=0) )
  }
  
  z=matrix(0,n,m)
  z[y==0]=1
  data <- list ("n","betamu1","betamu2","betamu3","t","x1","x2","R","rateobs",
                "y","z","surt","surt.cen","Sigma1","Sigma2","Sigma3","U0")
  U=matrix(0,n,2)
  inits <- function(){
    list(beta1=rep(0,4),beta2=rep(0,3),beta3=rep(0,2),r1=1,r2=-1,p=1,r=1,
         U=U)
  }
  parameters <- c("beta1","beta2","beta3","r1","r2","tau","sigma","p","r")
  time00 = Sys.time()
  my.sim <- bugs(data, inits, parameters, model.file="model.file.txt",
                 n.chains=2, n.iter=10000)#,debug = T)
  
  
  print(my.sim,digits=3)

  
#################################### Results ##################################
  Inference for Bugs model at "model.file.txt", 
  Current: 2 chains, each with 10000 iterations (first 5000 discarded)
  Cumulative: n.sims = 10000 iterations saved
  mean    sd       2.5%        25%        50%        75%      97.5%  Rhat n.eff
  beta1[1]    1.882e+00 0.136  1.618e+00  1.794e+00  1.883e+00  1.977e+00  2.134e+00 1.001  2700
  beta1[2]   -1.021e+00 0.150 -1.309e+00 -1.121e+00 -1.024e+00 -9.200e-01 -7.190e-01 1.003   700
  beta1[3]    1.026e+00 0.091  8.630e-01  9.580e-01  1.025e+00  1.091e+00  1.201e+00 1.004  5400
  beta1[4]   -9.540e-01 0.147 -1.265e+00 -1.050e+00 -9.490e-01 -8.490e-01 -6.900e-01 1.023    76
  beta2[1]    9.070e-01 0.125  6.670e-01  8.200e-01  9.060e-01  9.850e-01  1.161e+00 1.006   330
  beta2[2]    1.123e+00 0.198  7.480e-01  9.850e-01  1.119e+00  1.257e+00  1.524e+00 1.005   920
  beta2[3]   -1.063e+00 0.096 -1.258e+00 -1.126e+00 -1.062e+00 -9.980e-01 -8.790e-01 1.011   150
  beta3[1]    4.430e-01 0.109  2.390e-01  3.690e-01  4.390e-01  5.120e-01  6.700e-01 1.002  3400
  beta3[2]   -4.640e-01 0.110 -6.920e-01 -5.350e-01 -4.610e-01 -3.870e-01 -2.570e-01 1.003   910
  r1          1.051e+00 0.242  6.670e-01  8.780e-01  1.018e+00  1.191e+00  1.627e+00 1.004  1200
  r2         -1.067e+00 0.259 -1.691e+00 -1.209e+00 -1.032e+00 -8.860e-01 -6.660e-01 1.008   210
  tau[1,1]    1.429e+00 0.276  9.910e-01  1.234e+00  1.395e+00  1.589e+00  2.080e+00 1.001 10000
  tau[1,2]   -6.340e-01 0.239 -1.201e+00 -7.640e-01 -6.040e-01 -4.670e-01 -2.570e-01 1.005  2500
  tau[2,1]   -6.340e-01 0.239 -1.201e+00 -7.640e-01 -6.040e-01 -4.670e-01 -2.570e-01 1.005  2500
  tau[2,2]    1.298e+00 0.288  8.570e-01  1.090e+00  1.257e+00  1.462e+00  1.973e+00 1.006   410
  sigma[1,1]  9.320e-01 0.145  6.810e-01  8.290e-01  9.200e-01  1.023e+00  1.237e+00 1.001  3100
  sigma[1,2]  4.470e-01 0.129  2.180e-01  3.580e-01  4.390e-01  5.300e-01  7.180e-01 1.001  5000
  sigma[2,1]  4.470e-01 0.129  2.180e-01  3.580e-01  4.390e-01  5.300e-01  7.180e-01 1.001  5000
  sigma[2,2]  1.036e+00 0.189  7.040e-01  9.010e-01  1.022e+00  1.154e+00  1.437e+00 1.004   480
  p           1.005e+00 0.082  8.580e-01  9.470e-01  9.990e-01  1.056e+00  1.178e+00 1.006   290
  r           3.145e+00 0.421  2.387e+00  2.852e+00  3.122e+00  3.407e+00  4.043e+00 1.001 10000
  deviance    4.705e+06 0.000  4.705e+06  4.705e+06  4.705e+06  4.705e+06  4.705e+06 1.000     1
  
  For each parameter, n.eff is a crude measure of effective sample size,
  and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
  
  DIC info (using the rule, pD = Dbar-Dhat)
  pD = 339.400 and DIC = 4705000.000
  DIC is an estimate of expected predictive error (lower deviance is better).