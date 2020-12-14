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
  Pre = 6.01, Running = 75.3, Post = 3.82, Total = 85.1 
Fixed effects:
  mean     sd 0.025quant 0.5quant 0.975quant   mode kld
mu1       0.621  0.131      0.358    0.623      0.872  0.628   0
mu2      -0.732  0.144     -1.017   -0.731     -0.452 -0.730   0
mu3       0.530  0.116      0.302    0.529      0.759  0.529   0
l.SEX     0.923  0.110      0.711    0.922      1.142  0.919   0
l.PREVOI -1.066  0.170     -1.403   -1.066     -0.735 -1.064   0
s.SEX    -0.440  0.128     -0.694   -0.439     -0.192 -0.436   0
z.SEX    -0.956  0.101     -1.160   -0.955     -0.761 -0.952   0
l.SS1     2.420  0.433      1.568    2.421      3.267  2.423   0
l.SS2    -9.667  1.899    -13.387   -9.670     -5.936 -9.674   0
l.SS3    12.080  3.971      4.291   12.077     19.876 12.072   0
l.SS4    -5.800  6.040    -17.611   -5.816      6.092 -5.849   0
l.SS5     4.646  7.912    -10.920    4.657     20.138  4.679   0
l.SS6    -7.645  9.266    -25.884   -7.629     10.488 -7.596   0
l.SS7     9.289 10.203    -10.770    9.298     29.279  9.316   0
l.SS8    -8.564 11.481    -31.222   -8.525     13.853 -8.446   0
l.SS9    -3.147 30.610    -63.245   -3.148     56.900 -3.147   0
l.SS10   -0.011 31.623    -62.098   -0.012     62.023 -0.011   0
z.SS1     3.953  1.222      1.562    3.949      6.359  3.943   0
z.SS2    -9.219  5.069    -19.195   -9.211      0.704 -9.195   0
z.SS3     5.166  9.353    -13.206    5.169     23.505  5.175   0
z.SS4    -1.805 11.440    -24.239   -1.815     20.661 -1.833   0
z.SS5     5.656 11.971    -17.802    5.640     29.183  5.609   0
z.SS6    -5.308 12.017    -28.887   -5.314     18.280 -5.324   0
z.SS7     3.606 12.104    -20.148    3.602     27.361  3.594   0
z.SS8    -2.531 12.807    -27.683   -2.529     22.585 -2.524   0
z.SS9    -1.070 30.674    -61.293   -1.071     59.102 -1.070   0
z.SS10   -0.004 31.623    -62.090   -0.005     62.030 -0.004   0

Random effects:
  Name	  Model
U11 IID2D model
U21 Copy
U12 Copy
U22 Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant   mode
alpha parameter for weibullsurv[3]  1.002 0.036      0.935    1.001      1.075  0.999
Precision for U11 (component 1)     0.911 0.136      0.669    0.903      1.200  0.890
Precision for U11 (component 2)     1.265 0.270      0.786    1.251      1.836  1.230
Rho1:2 for U11                      0.524 0.110      0.290    0.530      0.719  0.542
Beta for U12                        1.044 0.191      0.671    1.043      1.420  1.041
Beta for U22                       -1.065 0.242     -1.536   -1.067     -0.587 -1.072

Expected number of effective parameters(stdev): 268.90(10.58)
Number of equivalent replicates : 12.66 

Marginal log-Likelihood:  -9909.57 
Posterior marginals for the linear predictor and
the fitted values are computed