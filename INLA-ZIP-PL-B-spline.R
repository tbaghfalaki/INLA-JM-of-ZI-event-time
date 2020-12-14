Data=read.table("simulated_data_ZIP_PL.txt",header=TRUE)
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
joint.inla <- inla(formula, family = c("zeroinflatedpoisson0","binomial","weibullsurv"),
                   data = joint.data, control.compute=list(dic=TRUE),control.family = list(
                     list(hyper = list(prob = list(initial = -10,
                                                   fixed = TRUE))),list(),list()))


summary(joint.inla)
#################################### Results ##################################



Call:
  c("inla(formula = formula, family = c(\"zeroinflatedpoisson0\", \"binomial\", ", " \"weibullsurv\"), 
   data = joint.data, control.compute = list(dic = TRUE), ", " control.family = list(list(hyper = 
   list(prob = list(initial = -10, ", " fixed = TRUE))), list(), list()))") 
Time used:
  Pre = 8.14, Running = 99.8, Post = 10.9, Total = 119 
Fixed effects:
  mean     sd 0.025quant 0.5quant 0.975quant   mode kld
mu1      -0.048  1.924     -3.825   -0.048      3.726 -0.048   0
mu2      -0.215  1.922     -3.989   -0.215      3.555 -0.215   0
mu3       0.527  0.123      0.296    0.523      0.784  0.515   0
l.SEX     0.927  0.110      0.715    0.926      1.147  0.923   0
l.PREVOI -1.053  0.172     -1.393   -1.052     -0.719 -1.050   0
s.SEX    -0.436  0.128     -0.693   -0.434     -0.191 -0.430   0
z.SEX    -0.964  0.103     -1.170   -0.963     -0.766 -0.960   0
l.SS1    -2.972 23.016    -48.160   -2.973     42.178 -2.972   0
l.SS2    -0.168 17.778    -35.072   -0.169     34.706 -0.168   0
l.SS3     0.491 20.965    -40.671    0.490     41.618  0.491   0
l.SS4    -0.906 13.821    -28.042   -0.907     26.206 -0.906   0
l.SS5     0.145  8.065    -15.689    0.145     15.967  0.145   0
l.SS6    -0.312 15.627    -30.993   -0.312     30.343 -0.312   0
l.SS7     0.280 20.612    -40.189    0.280     40.715  0.280   0
l.SS8     0.453 15.628    -30.230    0.453     31.111  0.453   0
l.SS9    -0.607  8.085    -16.480   -0.607     15.253 -0.607   0
l.SS10    1.244 13.832    -25.912    1.244     28.378  1.244   0
l.SS11   -0.413 20.968    -41.579   -0.413     40.720 -0.413   0
l.SS12   -0.150 17.783    -35.064   -0.150     34.735 -0.150   0
l.SS13    2.929 23.016    -42.259    2.928     48.080  2.929   0
z.SS1    -0.810 23.024    -46.015   -0.811     44.356 -0.810   0
z.SS2    -0.232 17.783    -35.146   -0.232     34.653 -0.232   0
z.SS3     0.254 20.967    -40.910    0.253     41.384  0.254   0
z.SS4     0.118 13.831    -27.037    0.117     27.250  0.118   0
z.SS5    -0.236  8.093    -16.126   -0.236     15.641 -0.236   0
z.SS6    -0.688 15.635    -31.385   -0.688     29.984 -0.687   0
z.SS7    -0.137 20.615    -40.612   -0.137     40.305 -0.137   0
z.SS8     0.197 15.636    -30.503    0.196     30.870  0.197   0
z.SS9    -0.397  8.101    -16.302   -0.397     15.495 -0.397   0
z.SS10    0.350 13.838    -26.819    0.350     27.497  0.350   0
z.SS11    0.001 20.969    -41.169    0.000     41.136  0.001   0
z.SS12    0.007 17.788    -34.918    0.006     34.902  0.007   0
z.SS13    0.794 23.024    -44.411    0.793     45.961  0.794   0

Random effects:
  Name	  Model
U11 IID2D model
U21 Copy
U12 Copy
U22 Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant   mode
alpha parameter for weibullsurv[3]  1.010 0.042      0.935    1.007      1.101  0.997
Precision for U11 (component 1)     0.897 0.136      0.659    0.887      1.192  0.867
Precision for U11 (component 2)     1.221 0.255      0.808    1.190      1.808  1.128
Rho1:2 for U11                      0.527 0.107      0.317    0.527      0.730  0.518
Beta for U12                        0.983 0.187      0.640    0.973      1.373  0.936
Beta for U22                       -0.984 0.230     -1.466   -0.970     -0.567 -0.920

Expected number of effective parameters(stdev): 271.95(10.99)
Number of equivalent replicates : 12.52 

Marginal log-Likelihood:  -9933.71 
Posterior marginals for the linear predictor and
the fitted values are computed
