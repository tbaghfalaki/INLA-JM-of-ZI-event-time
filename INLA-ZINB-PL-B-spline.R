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
#################################### Results ##################################
Call:
  c("inla(formula = formula, family = c(\"zeroinflatednbinomial0\", ", " \"binomial\", 
   \"weibullsurv\"), data = joint.data, control.compute = list(dic = TRUE), ", " control.family = 
   list(list(hyper = list(prob = list(initial = -10, ", " fixed = TRUE))), list(), list()))") 
Time used:
  Pre = 9.05, Running = 155, Post = 7.34, Total = 172 
Fixed effects:
  mean     sd 0.025quant 0.5quant 0.975quant   mode kld
mu1       0.166  1.927     -3.618    0.166      3.946  0.166   0
mu2      -0.018  1.922     -3.792   -0.018      3.752 -0.018   0
mu3       0.592  0.118      0.366    0.590      0.827  0.586   0
l.SEX     0.952  0.133      0.695    0.951      1.219  0.947   0
l.PREVOI -1.054  0.206     -1.462   -1.052     -0.652 -1.050   0
s.SEX    -0.450  0.120     -0.691   -0.448     -0.221 -0.443   0
z.SEX    -0.950  0.108     -1.168   -0.948     -0.744 -0.944   0
l.SS1    -2.845 23.023    -48.046   -2.845     42.319 -2.845   0
l.SS2    -0.345 17.784    -35.260   -0.345     34.541 -0.345   0
l.SS3     0.842 20.967    -40.324    0.841     41.973  0.842   0
l.SS4    -0.688 13.832    -27.844   -0.688     26.445 -0.688   0
l.SS5    -1.358  8.108    -17.277   -1.359     14.547 -1.358   0
l.SS6     1.027 15.646    -29.690    1.027     31.719  1.027   0
l.SS7     0.588 20.617    -39.890    0.587     41.031  0.588   0
l.SS8    -0.507 15.647    -31.227   -0.507     30.188 -0.507   0
l.SS9     0.658  8.135    -15.313    0.658     16.616  0.658   0
l.SS10    0.155 13.845    -27.027    0.154     27.314  0.155   0
l.SS11   -0.687 20.969    -41.856   -0.687     40.448 -0.687   0
l.SS12    0.733 17.791    -34.196    0.733     35.633  0.733   0
l.SS13    2.887 23.023    -42.315    2.886     48.050  2.887   0
z.SS1    -0.579 23.024    -45.783   -0.580     44.588 -0.579   0
z.SS2     0.057 17.784    -34.859    0.056     34.943  0.057   0
z.SS3     0.133 20.967    -41.033    0.132     41.264  0.133   0
z.SS4    -0.271 13.832    -27.427   -0.272     26.862 -0.271   0
z.SS5    -0.154  8.097    -16.050   -0.154     15.729 -0.154   0
z.SS6    -0.640 15.637    -31.341   -0.640     30.036 -0.640   0
z.SS7    -0.237 20.616    -40.714   -0.238     40.206 -0.237   0
z.SS8     0.454 15.638    -30.250    0.453     31.132  0.454   0
z.SS9     0.242  8.107    -15.674    0.242     16.145  0.242   0
z.SS10   -0.629 13.839    -27.800   -0.629     26.520 -0.629   0
z.SS11   -0.007 20.969    -41.176   -0.008     41.127 -0.007   0
z.SS12    0.863 17.789    -34.062    0.862     35.758  0.863   0
z.SS13    0.715 23.024    -44.490    0.714     45.882  0.715   0

Random effects:
  Name	  Model
U11 IID2D model
U21 Copy
U12 Copy
U22 Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant   mode
size for nbinomial zero-inflated observations  0.879 0.180      0.570    0.864      1.276  0.837
alpha parameter for weibullsurv[3]             1.061 0.083      0.927    1.050      1.249  1.018
Precision for U11 (component 1)                0.727 0.125      0.509    0.718      0.999  0.702
Precision for U11 (component 2)                1.203 0.260      0.779    1.172      1.796  1.112
Rho1:2 for U11                                 0.458 0.114      0.217    0.464      0.661  0.477
Beta for U12                                   0.743 0.148      0.461    0.739      1.043  0.725
Beta for U22                                  -0.716 0.198     -1.118   -0.711     -0.339 -0.693

Expected number of effective parameters(stdev): 252.36(11.71)
Number of equivalent replicates : 12.57 

Marginal log-Likelihood:  -9923.58 
Posterior marginals for the linear predictor and
the fitted values are computed

