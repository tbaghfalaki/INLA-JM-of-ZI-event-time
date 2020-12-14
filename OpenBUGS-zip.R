rm(list=ls())
library(R2OpenBUGS)
Data=read.table("simulated_data_ZIP.txt",header=TRUE)
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
      ll[i,j]<-z[i,j]*log(pi[i,j]) + (1-z[i,j])*(log(1-pi[i,j])+y[i,j]*log(lambda[i,j])- 
                lambda[i,j] - loggam(y[i,j]+1)-log(1-exp(-lambda[i,j])))
               
      log(lambda[i, j])<-beta1[1]+beta1[2]*t[j]+beta1[3]*x1[i]+beta1[4]*x2[i]+U[i,1]
      logit(pi[i, j])<-beta2[1]+beta2[2]*t[j]+beta2[3]*x1[i]+U[i,2]
    }
    surt[i] ~ dweib(p,mut[i]) I(surt.cen[i],)      

    log(mut[i])<-beta3[1]+beta3[2]*x1[i]+r1*U[i, 1]+r2*U[i, 2]
    U[i,1:2] ~ dmnorm(U0[],tau[,])   
  }  
  
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
  list(beta1=rep(0,4),beta2=rep(0,3),beta3=rep(0,2),r1=1,r2=-1,p=1,
       U=U)
}
parameters <- c("beta1","beta2","beta3","r1","r2","tau","sigma","p")
time00 = Sys.time()
my.sim <- bugs(data, inits, parameters, model.file="model.file.txt",
               n.chains=2, n.iter=10000)#,debug = T)


print(my.sim,digits=3)


#################################### Results ##################################
Inference for Bugs model at "model.file.txt", 
Current: 2 chains, each with 10000 iterations (first 5000 discarded)
Cumulative: n.sims = 10000 iterations saved
mean    sd       2.5%        25%        50%        75%      97.5%  Rhat n.eff
beta1[1]    1.992e+00 0.160  1.606e+00  1.912e+00  2.013e+00  2.102e+00  2.242e+00 1.196    19
beta1[2]   -1.006e+00 0.057 -1.121e+00 -1.043e+00 -1.006e+00 -9.700e-01 -8.900e-01 1.019   170
beta1[3]    1.023e+00 0.115  7.990e-01  9.460e-01  1.028e+00  1.089e+00  1.246e+00 1.746     5
beta1[4]   -7.440e-01 0.221 -1.023e+00 -9.220e-01 -7.970e-01 -6.410e-01 -2.400e-01 1.833     4
beta2[1]    1.083e+00 0.127  8.310e-01  9.960e-01  1.083e+00  1.168e+00  1.341e+00 1.008   240
beta2[2]    9.060e-01 0.199  5.170e-01  7.710e-01  9.030e-01  1.038e+00  1.307e+00 1.007   670
beta2[3]   -1.024e+00 0.097 -1.214e+00 -1.091e+00 -1.023e+00 -9.570e-01 -8.410e-01 1.015   110
beta3[1]    6.470e-01 0.145  3.810e-01  5.440e-01  6.390e-01  7.410e-01  9.520e-01 1.016   100
beta3[2]   -5.340e-01 0.129 -7.910e-01 -6.170e-01 -5.320e-01 -4.450e-01 -2.910e-01 1.128    17
r1          1.473e+00 0.330  8.890e-01  1.237e+00  1.450e+00  1.677e+00  2.183e+00 1.018   580
r2         -1.228e+00 0.337 -1.937e+00 -1.444e+00 -1.206e+00 -9.880e-01 -6.420e-01 1.006   510
tau[1,1]    1.590e+00 0.287  1.095e+00  1.387e+00  1.563e+00  1.765e+00  2.218e+00 1.004  4400
tau[1,2]   -9.390e-01 0.283 -1.574e+00 -1.111e+00 -9.090e-01 -7.400e-01 -4.640e-01 1.006   310
tau[2,1]   -9.390e-01 0.283 -1.574e+00 -1.111e+00 -9.090e-01 -7.400e-01 -4.640e-01 1.006   310
tau[2,2]    1.651e+00 0.374  1.030e+00  1.390e+00  1.604e+00  1.872e+00  2.498e+00 1.020   160
sigma[1,1]  9.830e-01 0.137  7.510e-01  8.860e-01  9.710e-01  1.067e+00  1.289e+00 1.002  1500
sigma[1,2]  5.580e-01 0.127  3.240e-01  4.700e-01  5.510e-01  6.370e-01  8.230e-01 1.001  2800
sigma[2,1]  5.580e-01 0.127  3.240e-01  4.700e-01  5.510e-01  6.370e-01  8.230e-01 1.001  2800
sigma[2,2]  9.650e-01 0.189  6.420e-01  8.310e-01  9.460e-01  1.079e+00  1.387e+00 1.016   280
p           1.056e+00 0.101  8.760e-01  9.870e-01  1.050e+00  1.119e+00  1.272e+00 1.007   950
deviance    4.762e+06 0.000  4.762e+06  4.762e+06  4.762e+06  4.762e+06  4.762e+06 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = Dbar-Dhat)
pD = 347.600 and DIC = 4762000.000
DIC is an estimate of expected predictive error (lower deviance is better).