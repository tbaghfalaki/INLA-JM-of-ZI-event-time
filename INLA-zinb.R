  rm(list=ls())
  library(INLA)
  Data=read.table("simulated_data_ZINB.txt",header=TRUE)
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
  
  
  
  
  formula = Y ~ mu + l.TIME  + l.SEX +l.PREVOI+  s.SEX - 1 + z.TIME+z.SEX+
    f(U11 , model="iid2d",n=2*n2) + 
    f(U21,  copy="U11") + f(U12, copy="U11", fixed= FALSE, param=c(0,0.01),initial = -0.2) + 
    f(U22, copy="U11",fixed = FALSE, param = c(0,0.01), initial = -1.6)
  
  
  
  formula = Y ~ mu + l.TIME  + l.SEX +l.PREVOI+  s.SEX - 1 + z.TIME+z.SEX+
    f(U11 , model="iid2d",n=2*n2) + 
    f(U21,  copy="U11") + f(U12, copy="U11", fixed= FALSE, param=c(0.01,0.01)) + 
    f(U22, copy="U11",fixed = FALSE, param = c(0.01,0.01))
  
  formula = Y ~ mu + l.TIME  + l.SEX +l.PREVOI+  s.SEX - 1 + z.TIME+z.SEX+
    f(U11 , model="iid2d",n=2*n2) + 
    f(U21,  copy="U11") + f(U12, copy="U11", fixed= FALSE, param=c(1,1)) + 
    f(U22, copy="U11",fixed = FALSE, param = c(1,1))
  
  
  
  #time00 = Sys.time()
  #joint.inla <- inla(formula, family = c("zeroinflatednbinomial0","binomial","exponentialsurv"),
   #                  data = joint.data, control.compute=list(dic=TRUE),control.family = list(
  #                     list(hyper = list(prob = list(initial = -10,
   #                                                  fixed = TRUE))),list(),list()))
  #Finaltime = Sys.time() - time00
  #print(Finaltime)
  #round(joint.inla$summary.fixed, 4)
  # summary(joint.inla)
  
  
  time00 = Sys.time()
  joint.inla <- inla(formula, family = c("zeroinflatednbinomial0","binomial","weibull.surv"),
                     data = joint.data, control.compute=list(dic=TRUE),control.family = list(
                       list(hyper = list(prob = list(initial = -10,
                                                     fixed = TRUE))),list(),list()))
  Finaltime = Sys.time() - time00
  print(Finaltime)
  summary(joint.inla)
  
  #################################### Results ##################################
  Call:
    c("inla(formula = formula, family = c(\"zeroinflatednbinomial0\", ", " \"binomial\", \"weibull.surv\"), data = 
   joint.data, control.compute = list(dic = TRUE), ", " control.family = list(list(hyper = list(prob = list(initial = 
   -10, ", " fixed = TRUE))), list(), list()))") 
  Time used:
    Pre = 5.1, Running = 149, Post = 4.14, Total = 158 
  Fixed effects:
    mean    sd 0.025quant 0.5quant 0.975quant   mode kld
  mu1       1.938 0.133      1.674    1.938      2.198  1.940   0
  mu2       0.916 0.125      0.674    0.915      1.164  0.913   0
  mu3       0.433 0.096      0.246    0.433      0.623  0.432   0
  l.TIME   -1.056 0.139     -1.329   -1.056     -0.783 -1.055   0
  l.SEX     1.029 0.097      0.840    1.028      1.219  1.028   0
  l.PREVOI -0.975 0.143     -1.256   -0.975     -0.695 -0.975   0
  s.SEX    -0.442 0.098     -0.638   -0.441     -0.252 -0.438   0
  z.TIME    1.150 0.196      0.767    1.150      1.538  1.148   0
  z.SEX    -1.079 0.096     -1.273   -1.076     -0.895 -1.072   0
  
  Random effects:
    Name	  Model
  U11 IID2D model
  U21 Copy
  U12 Copy
  U22 Copy
  
  Model hyperparameters:
    mean    sd 0.025quant 0.5quant 0.975quant   mode
  size for nbinomial zero-inflated observations  3.428 0.461      2.582    3.409      4.389  3.384
  alpha parameter for weibullsurv[3]             0.960 0.060      0.853    0.955      1.089  0.944
  Precision for U11 (component 1)                1.166 0.173      0.865    1.152      1.543  1.125
  Precision for U11 (component 2)                1.007 0.199      0.700    0.978      1.472  0.915
  Rho1:2 for U11                                 0.448 0.105      0.240    0.449      0.647  0.444
  Beta for U12                                   0.877 0.159      0.581    0.871      1.203  0.849
  Beta for U22                                  -0.856 0.152     -1.170   -0.849     -0.573 -0.824
  
  Expected number of effective parameters(stdev): 338.84(14.03)
  Number of equivalent replicates : 14.76 
  
  Marginal log-Likelihood:  -20359.42 
  Posterior marginals for the linear predictor and
  the fitted values are computed
  