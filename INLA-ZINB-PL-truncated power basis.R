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




SS=tp(timevar,degree=3, k=15)$Z # truncated power function
#SS
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
#################################### Results ##################################
Call:
  c("inla(formula = formula, family = c(\"zeroinflatednbinomial0\", ", " \"binomial\", 
   \"weibullsurv\"), data = joint.data, control.compute = list(dic = TRUE), ", " control.family = 
   list(list(hyper = list(prob = list(initial = -10, ", " fixed = TRUE))), list(), list()))") 
Time used:
  Pre = 30.5, Running = 84.7, Post = 2.47, Total = 118 
Fixed effects:
  mean     sd 0.025quant 0.5quant 0.975quant    mode kld
mu1        0.873  0.182      0.509    0.875      1.224   0.880   0
mu2       -0.558  0.142     -0.838   -0.557     -0.281  -0.556   0
mu3        0.593  0.118      0.367    0.591      0.828   0.588   0
l.SEX      0.955  0.133      0.698    0.953      1.221   0.950   0
l.PREVOI  -1.047  0.205     -1.454   -1.046     -0.647  -1.043   0
s.SEX     -0.451  0.120     -0.692   -0.449     -0.222  -0.444   0
z.SEX     -0.943  0.107     -1.158   -0.942     -0.739  -0.938   0
l.SS1      2.007  1.218     -0.363    1.999      4.421   1.983   0
l.SS2     -6.942  5.157    -17.160   -6.913      3.107  -6.857   0
l.SS3      2.480  9.786    -16.724    2.471     21.714   2.455   0
l.SS4     12.741 12.654    -12.047   12.718     37.634  12.673   0
l.SS5    -14.240 14.124    -41.957  -14.246     13.489 -14.258   0
l.SS6      4.376 14.587    -24.309    4.391     32.955   4.421   0
l.SS7     -2.437 14.607    -31.114   -2.440     26.222  -2.442   0
l.SS8      2.675 20.414    -37.409    2.675     42.717   2.679   0
l.SS9      2.749 26.776    -49.823    2.749     55.273   2.751   0
l.SS10     1.064 31.003    -59.806    1.063     61.884   1.064   0
z.SS1      5.461  1.261      3.001    5.455      7.951   5.443   0
z.SS2    -14.757  5.255    -25.125  -14.740     -4.496 -14.705   0
z.SS3     11.852  9.761     -7.307   11.850     31.005  11.847   0
z.SS4     -4.866 11.987    -28.358   -4.881     18.691  -4.911   0
z.SS5      7.426 12.541    -17.146    7.409     32.073   7.375   0
z.SS6     -5.342 12.675    -30.216   -5.346     19.533  -5.354   0
z.SS7     -4.478 12.808    -29.624   -4.478     20.646  -4.477   0
z.SS8      7.296 20.142    -32.249    7.295     46.811   7.294   0
z.SS9      5.680 26.539    -46.426    5.679     57.743   5.680   0
z.SS10     2.126 30.972    -58.684    2.125     62.884   2.126   0

Random effects:
  Name	  Model
U11 IID2D model
U21 Copy
U12 Copy
U22 Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant   mode
size for nbinomial zero-inflated observations  0.887 0.179      0.582    0.872      1.284  0.843
alpha parameter for weibullsurv[3]             1.058 0.081      0.923    1.049      1.239  1.023
Precision for U11 (component 1)                0.739 0.128      0.523    0.726      1.027  0.699
Precision for U11 (component 2)                1.245 0.274      0.783    1.219      1.855  1.171
Rho1:2 for U11                                 0.449 0.115      0.214    0.451      0.663  0.451
Beta for U12                                   0.738 0.144      0.469    0.732      1.035  0.711
Beta for U22                                  -0.723 0.199     -1.128   -0.716     -0.346 -0.694

Expected number of effective parameters(stdev): 247.93(11.81)
Number of equivalent replicates : 12.79 

Marginal log-Likelihood:  -9897.30 
Posterior marginals for the linear predictor and
the fitted values are computed

