#setwd("~/Dropbox (HMS)/Documents/Research/RL/SS RL/SSRL_code/")

library(glmnet)
library(MASS)

#### Functions
expit <- function(x) exp(x)/(1+exp(x)) 

###########################
##### Data Simulations ####
###########################

gen_EHR.data <- function(size,phi1,beta1,delta1,phi2,beta2,psi2,seed,return.Vstar,gamma1,gamma2){
  set.seed(seed)
  # True Q function working models:
  S1 <- ~ O1+O2+O3+O4+O5+O6+A1*(O2+O3+O4+O5+O6); S2 <- ~ Y1+O1+O2+O3+O4+O5+O6+A1+Z.21+Z.21+Z.22+A2*(O1+O2+O3+O6+A1+Z.22+Z.21)# high dimension
  # True propensity score functions:
  PS1 <- ~ O1+O2+O3; PS2 <- ~ Y1+O1+A1+Z.21+O2
  
  O.var = diag(1,6,6)
  O <- data.frame(mvrnorm(n = size, mu=rep(0,6), Sigma=O.var, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
  colnames(O) <- gsub('X','O',colnames(O))
  H1check <- model.matrix.lm(PS1,O)
  phi1 <- c(-.1,1,-1,.1)
  prob_A1 <- expit(tcrossprod(t(phi1),H1check))
  A1 <- rbinom(size,1,prob_A1)
  
  beta1 <- c(.5,.2,-1,-1,.1,-.1,.1)
  gamma1 <- c(1,-2,-2,-.1,.1,-1.5)
  H1.Q <- model.matrix.lm(S1,cbind(O,A1=A1))
  if(return.Vstar){
    # compute d1.bar and set A1=d1bar
    interaction.terms1 <- H1.Q[,c("(Intercept)",gsub(':A1','',colnames(H1.Q)[grep(':A1',colnames(H1.Q))]))]
    A1 <- 1*(tcrossprod(interaction.terms1,t(gamma1))>0)
    H1.Q <- model.matrix.lm(S1,data.frame(O,A1))
  }
  
  mean_Y1 <- tcrossprod(t(c(beta1,gamma1)),H1.Q)
  Y1 <- rbinom(size,1,expit(mean_Y1))
  
  
  Z.21 <- as.numeric(1.25*O[,1]*A1+rnorm(n=size,mean=0,sd=1) >0)
  Z.22 <- as.numeric(-1.75*O[,2]*A1+rnorm(n=size,mean=0,sd=1) >0)
  
  H2check <- model.matrix.lm(PS2,data.frame(O1=O[,1],O2=O[,2],Y1=Y1,A1=A1,Z.21=Z.21))
  mean_A2 <- tcrossprod(t(phi2[1:ncol(H2check)]),H2check)
  if (phi24 != 0){
    # model + a misspecified term which is positive when phi24 is different than zero
    phi24.vec <- t(c(phi24,phi24)/as.numeric(sqrt(crossprod(c(phi24,phi24)))))
    mean_A2 <- mean_A2 + tcrossprod(phi24.vec,H2check[,c('Y1','Z.21')])*Y1*sin(apply(H2check[,c('Y1','Z.21')],1,crossprod)/(Y1+1))
  }
  prob_A2 <- expit(mean_A2)
  A2 <- rbinom(size,1,prob_A2)
  
  beta2 <- c(1,beta1,.25,-1,-.5)
  gamma2 <- c(1,.1,-.1,.1,-.1,.25,-1,-.5)
  H2.Q <- model.matrix.lm(S2,cbind(O,A1=A1,A2=A2,Y1,Z.21,Z.22))
  if(return.Vstar){
    # compute d2.bar and set A2=d2bar
    interaction.terms2 <- H2.Q[,c("(Intercept)",gsub(':A2','',colnames(H2.Q)[grep(':A2',colnames(H2.Q))]))]
    A2 <- 1*(tcrossprod(interaction.terms2,t(gamma2))>0)
    H2.Q <- model.matrix.lm(S2,cbind(O,A1=A1,A2=A2,Y1,Z.21,Z.22))
  }
  
  mean_Y2 <- tcrossprod(t(c(beta2,gamma2)),H2.Q)
  if (beta26 != 0){
    # model + a misspecified term which is positive when beta26 is different than zero
    beta26.vec <- t(c(beta26,beta26)/as.numeric(sqrt(crossprod(c(beta26,beta26)))))
    mean_Y2 <- mean_Y2 + tcrossprod(beta26.vec,H2.Q[,c('Z.21','Z.22')])*Y1*sin(apply(H2.Q[,c('Z.21','Z.22')],1,crossprod)/(Y1+1))
  }
  Y2 <- rbinom(size,1,expit(mean_Y2))
  # Surrogates
  feat1 <- cbind(1,Y1)
  GAMMA1 <- matrix(1,nrow=3,ncol=ncol(feat1))
  feat2 <- cbind(1,Y2)
  GAMMA2 <- matrix(1,nrow=3,ncol=ncol(feat2))
   
  W1.var <- W2.var <- diag(3)/1000
  
  W1 <- floor(feat1%*%t(GAMMA1)+mvrnorm(n = size, mu=rep(0,3), Sigma=W1.var))
  colnames(W1) <- paste('W1.',1:ncol(W1),sep='')
  #summary(lm(Y1~W1))
  
  W2 <- floor(feat2%*%t(GAMMA2)+mvrnorm(n = size, mu=rep(0,3), Sigma=W2.var))
  colnames(W2) <- paste('W2.',1:ncol(W2),sep='')
  #summary(lm(Y2~W2))
  df <- cbind(O,Z.21,Z.22,A1,A2,Y1,Y2,W1,W2)
  return(df)
}

gen_synth.data <- function(sim.setting,size,phi1,beta1,delta1,phi2,beta2,psi2,seed,return.Vstar,gamma1,gamma2){
  set.seed(seed)
  # True Q functions:
  S1 <- ~ O1 + A1*O1; S2 <- ~ Y1 + O1*A1 + O2 +I((Y1*O2^2)*sin(1/((Y1+1)*O2^2))) +  (A1+O2)*A2 
  
  # True propensity score functions:
  PS1 <- ~ O1; PS2 <- ~ O1+A1+O2+A1*O1+I(sin(O2^2))+Y1
  
  O1 <- rbinom(size,1,.5)
  H1check <- model.matrix.lm(PS1,data.frame(O1))
  prob_A1 <- expit(tcrossprod(t(phi1),H1check))
  A1 <- rbinom(size,1,prob_A1)

  H1.Q <- model.matrix.lm(S1,data.frame(O1,A1))
  # compute d1.bar and set A1=d1bar
  if(return.Vstar){
    interaction.terms1 <- H1.Q[,c("(Intercept)",gsub(':A1','',colnames(H1.Q)[grep(':A1',colnames(H1.Q))]))]
    A1 <- 1*(tcrossprod(interaction.terms1,t(gamma1))>0)
    H1.Q <- model.matrix.lm(S1,data.frame(O1,A1))
  }
  
  mean_Y1 <- t(tcrossprod(beta1,H1.Q))
  
  if(sim.setting=='binary'){
    Y1 <- rbinom(size,1,expit(mean_Y1-1))
  }else{#sim.setting is continuous
    Y1 <- mean_Y1+rnorm(n=size,mean=0,sd=1)
  }
  Y1[abs(Y1)>5] <- sign(Y1[abs(Y1)>5])*5
  
  mean_O2 <- delta1['10']+delta1['11']*O1+delta1['12']*A1+delta1['13']*O1*A1
  O2 <- rnorm(n=size,mean=mean_O2,sd=sqrt(2))
  O2[abs(O2)>5] <- sign(O2[abs(O2)>5])*5
  
  H2check <- model.matrix.lm(PS2,data.frame(O1,A1,O2,Y1))
  prob_A2 <- expit(tcrossprod(t(phi2),H2check))
  A2 <- rbinom(size,1,prob_A2)
  
  H2.Q <- model.matrix.lm(S2,data.frame(Y1,O1,A1,O2,A2))
  if(return.Vstar){
    # compute d2.bar and set A2=d2bar
    interaction.terms2 <- H2.Q[,c("(Intercept)",gsub(':A2','',colnames(H2.Q)[grep(':A2',colnames(H2.Q))]))]
    A2 <- 1*(tcrossprod(interaction.terms2,t(gamma2))>0)
    H2.Q <- model.matrix.lm(S2,data.frame(Y1,O1,A1,O2,A2))
  }
  mean_Y2 <- t(tcrossprod(c(beta2[c(1,6,2,3,5,7)],psi2[1],beta2[4],psi2[-1]),H2.Q))
  if(sim.setting=='binary'){
    Y2 <- rbinom(size,1,expit(mean_Y2-3))
  }else{#sim.setting is continuous
    Y2 <- mean_Y2+rnorm(n=size,mean=0,sd=sqrt(2))
  }
  Y2[abs(Y2)>5] <- sign(Y2[abs(Y2)>5])*5
  sigma.W1 <- sqrt(ifelse(sim.setting=='continuous',.25,ifelse(sim.setting=='binary',.00001,2.2)))
  W1 <- floor(Y1 + rnorm(n=size,mean=0,sd=sigma.W1))
  summary(lm(Y1~W1))
  sigma.W2 <- sqrt(ifelse(sim.setting=='continuous',.35,ifelse(sim.setting=='binary',.00001,5)))
  W2 <- floor(Y2 + rnorm(n=size,mean=0,sd=sigma.W2))
  summary(lm(Y2~W2))
  df <- data.frame(O1,A1,Y1,O2,A2,Y2,W1,W2)
  return(df)
}

gen_data <- function(sim.setting,size,phi1,beta1,delta1,phi2,beta2,psi2,seed,return.Vstar=F,gamma1=NULL,gamma2=NULL){
  if (sim.setting %in% c('binary','continuous')){
    dat <- gen_synth.data(sim.setting,size,phi1,beta1,delta1,phi2,beta2,psi2,seed,return.Vstar,gamma1,gamma2) # low dimension
  }else if(sim.setting=='EHR'){
    dat <- gen_EHR.data(size,phi1,beta1,delta1,phi2,beta2,psi2,seed,return.Vstar,gamma1,gamma2)# high dimension
  }
  return(dat)
}

###################################################################
################# Estimate Q-functions parameters #################
###################################################################

MCparams <- function(sim.setting,S1,S2,PS1,PS2,phi1,beta1,delta1,phi2,beta2,psi2){
  n <- 1e6; N <- 10
  dat <- gen_data(sim.setting,size=n+N,phi1,beta1,delta1,phi2,beta2,psi2,seed=116687) 
  df_ls <- transform_data(S1,S2,PS1,PS2,dat,n,N)
  U_ls <- df_ls$U_ls; L_ls <- df_ls$L_ls
  results.sup <- multiroot(f = supQl.EE, start = rep(1,ncol(cbind(L_ls$H10,L_ls$H11,L_ls$Ys$Y1,L_ls$H20,L_ls$H21))),L_ls=L_ls)
  theta <- results.sup$root
  names(theta) <- colnames(cbind(L_ls$H10,L_ls$H11.A1,L_ls$Ys[,'Y1'],L_ls$H20,L_ls$H21.A2)); names(theta)[names(theta) == ""] <- 'Y1'
  the1.length <- ncol(cbind(L_ls$H10,L_ls$H11.A1))
  theta1 <- theta[1:the1.length]; theta2 <- theta[(the1.length+1):length(theta)]
  
  d1 <- crossprod(t(L_ls$H11),theta1[colnames(L_ls$H11.A1)])>0; gamma1 <- theta1[colnames(L_ls$H11.A1)]
  d2 <- crossprod(t(L_ls$H21),theta2[colnames(L_ls$H21.A2)])>0; gamma2 <- theta2[colnames(L_ls$H21.A2)]
  cat('\nOptimal treatment D*: (',mean(d1),mean(d2),')\n')

  dat <- gen_data(sim.setting,size=n+N,phi1,beta1,delta1,phi2,beta2,psi2,seed=116687,return.Vstar=T,gamma1,gamma2) 
  return(list(theta1=theta1,theta2=theta2,V.star=mean(dat$Y1+dat$Y2)))
}

################################################################
######################### Transform Data #######################
################################################################

transform_data <- function(S1,S2,PS1,PS2,dat,n,N,surrogates=NULL){
  S1 <- as.formula(S1); S2 <- as.formula(S2)
  dat.S1 <- model.matrix.lm(S1,dat)
  dat.S2 <- model.matrix.lm(S2,dat)
  dat.PS1 <- model.matrix.lm(PS1,dat)
  dat.PS2 <- model.matrix.lm(PS2,dat, na.action=na.pass)
  # Features for Q functions
  H21.A2 <- dat.S2[,grep('A2',colnames(dat.S2))]
  H21 <- dat.S2[,c('(Intercept)',gsub(':','',gsub('A2','',colnames(dat.S2)[grep('A2',colnames(dat.S2))[-1]])))]
  H20 <- dat.S2[,-grep('A2|Y',colnames(dat.S2))]
  
  H11.A1 <- dat.S1[,grep('A1',colnames(dat.S1))]
  H11 <- dat.S1[,c('(Intercept)',gsub(':','',gsub('A1','',colnames(dat.S1)[grep('A1',colnames(dat.S1))[-1]])))]
  H10 <- dat.S1[,-grep('A1',colnames(dat.S1))]
  # Outcomes/labels
  Ys <- dat[(N+1):(N+n),c('Y1','Y2')]
  # Vector of possible predictors for labeled data
  if(is.null(surrogates)){
    surrogates <- cbind(dat[,startsWith(colnames(dat),'W')])
    U_vec <- cbind(dat.S2,dat.S1,surrogates); U_vec <- U_vec[,unique(colnames(U_vec))]
    U_vec <- U_vec[,!grepl('(Intercept)|Y', colnames(U_vec))]; colnames(U_vec) <- gsub(':','x',colnames(U_vec))
  }else{
    U_vec <- dat[,surrogates]
  }
  
  # Uvecs:
  U_ls <- list(H10 = H10[1:N,], H11 = H11[1:N,], H11.A1 = H11.A1[1:N,], 
               H20 = H20[1:N,], H21 = H21[1:N,], H21.A2 = H21.A2[1:N,], 
               H1check = dat.PS1[1:N,], H2check = dat.PS2[1:N,],
               Uvec=U_vec[1:N,])
  L_ls <- list(Ys = dat[(N+1):(N+n),c('Y1','Y2')],
               H10 = H10[(N+1):(N+n),], H11 = H11[(N+1):(N+n),], H11.A1 = H11.A1[(N+1):(N+n),], 
               H20 = H20[(N+1):(N+n),], H21 = H21[(N+1):(N+n),], H21.A2 = H21.A2[(N+1):(N+n),],
               H1check = dat.PS1[(N+1):(N+n),], H2check = dat.PS2[(N+1):(N+n),],
               Uvec=U_vec[(N+1):(N+n),])
  return(list(U_ls=U_ls,L_ls=L_ls))
}

################################################################
################# Treamtent propensity score ####################
################################################################

PS.treat <- function(df_ls,PS1,PS2,label=T,small.n.correct=F){
  # A1
  if(small.n.correct){
    lambdas_vec <- list(exp(seq(0,1,by=.01))-1)[[1]]/n^.25
    compute_aic <- function(lambda){
      pi1.fit <- glmnet(x=as.matrix(df_ls$H1check[,-1]),y=df_ls$H11.A1[,'A1'],family="binomial",lambda=lambda)
      neg.2logLik <- deviance(pi1.fit)# for logistic reg. dev = 12logLik since the likelihood of the saturated model is equal to one
      k <- pi1.fit$df+1
      AIC <- neg.2logLik+2*k
      return(AIC)
    }
    AIC_vec <- unlist(lapply(lambdas_vec,compute_aic))
    lambda.star <- lambdas_vec[which.min(AIC_vec)]
    pi1.fit <- glmnet(x=as.matrix(df_ls$H1check[,-1]),y=df_ls$H11.A1[,'A1'],family="binomial",lambda=lambda.star)
    
    df_ls$pi <- data.frame(pi1=as.numeric(predict(pi1.fit,  newx = as.matrix(df_ls$H1check[,-1]),s = lambda.star,type = "response")))
  }else{
    pi1.fit <- glm(PS1,data=data.frame(A1=df_ls$H11.A1[,'A1'],df_ls$H1check),family='binomial')
    df_ls$pi <- data.frame(pi1=pi1.fit$fitted.values)
  }
  if(label){
    # A2
    if(small.n.correct){
      compute_aic <- function(lambda){
        pi2.fit <- glmnet(x=as.matrix(df_ls$H2check[,-1]),y=df_ls$H21.A2[,'A2'],family="binomial",lambda=lambda)
        neg.2logLik <- deviance(pi2.fit)# for logistic reg. dev = 12logLik since the likelihood of the saturated model is equal to one
        k <- pi2.fit$df+1
        AIC <- neg.2logLik+2*k
        return(AIC)
      }
      AIC_vec <- unlist(lapply(lambdas_vec,compute_aic))
      lambda.star <- lambdas_vec[which.min(AIC_vec)]
      pi2.fit <- glmnet(x=as.matrix(df_ls$H2check[,-1]),y=df_ls$H21.A2[,'A2'],family="binomial",lambda=lambda.star)
      df_ls$pi <- data.frame(df_ls$pi,pi2 = as.numeric(predict(pi2.fit,  newx = as.matrix(df_ls$H2check[,-1]),s = lambda.star,type = "response")))
      ###
      
    }else{
      pi2.fit <- glm(PS2,data=data.frame(A2=df_ls$H21.A2[,'A2'],df_ls$H2check),family='binomial')
      df_ls$pi <- data.frame(df_ls$pi,pi2 = pi2.fit$fitted.values)
    }
  }
  return(df_ls)
}

omegas <- function(df_ls,theta1,theta2,label=T){
  # Propensity score weights (omegas in paper)
  # optimal treatments
  d1 <- 1*(crossprod(t(df_ls$H11),theta1[colnames(df_ls$H11.A1)])>0)
  d2 <- 1*(crossprod(t(df_ls$H21),theta2[colnames(df_ls$H21.A2)])>0)
  if(label){
    df_ls$D.hat <- data.frame(d1=d1,d2=d2)
  }else{
    df_ls$D.hat.SSL <- data.frame(d1=d1,d2=d2)
  }
  # observed treatments
  A1 <- df_ls$H11.A1[,'A1']; A2 <- df_ls$H21.A2[,'A2']
  # propensity scores
  pi1 <- df_ls$pi[,'pi1']
  omega1 <- d1*A1/pi1+(1-d1)*(1-A1)/(1-pi1)
  df_ls$omegas <- data.frame(omega1=omega1,omega2=rep(NA,length(omega1)))
  if(dim(df_ls$pi)[2]>1){
    # propensity scores
    pi2 <- df_ls$pi[,'pi2']
    omega2 <- omega1*(d2*A2/pi2+(1-d2)*(1-A2)/(1-pi2))
    df_ls$omegas$omega2 <- omega2
  }
  return(df_ls)
}

VsupDR <- function(df_ls,theta1,theta2,label=T){
  # Qopt-functions & Vsup estimator
  Q1opt <- crossprod(t(df_ls$H10),theta1[colnames(df_ls$H10)]) + 
           crossprod(t(df_ls$H11),theta1[colnames(df_ls$H11.A1)])*(crossprod(t(df_ls$H11),theta1[colnames(df_ls$H11.A1)])>0)
  if(label){
    Q2opt <- crossprod(t(cbind(df_ls$Ys$Y1,df_ls$H20)),theta2[colnames(cbind(Y1=df_ls$Ys$Y1,df_ls$H20))]) + 
             crossprod(t(df_ls$H21),theta2[colnames(df_ls$H21.A2)])*(crossprod(t(df_ls$H21),theta2[colnames(df_ls$H21.A2)])>0)
    VsupDR.hat <- Q1opt + df_ls$omegas$omega1*(df_ls$Ys$Y1-Q1opt+Q2opt) + df_ls$omegas$omega2*(df_ls$Ys$Y2-Q2opt)
    df_ls$opt.SUP <- data.frame(Q1opt,Q2opt,VsupDR.hat)
  }else{
    # Compute imputed values of labeled data for influence function calculation (to derive SEs)
    Q2opt_ <- crossprod(t(df_ls$H20),theta2[colnames(df_ls$H20)]) + 
              crossprod(t(df_ls$H21),theta2[colnames(df_ls$H21.A2)])*(crossprod(t(df_ls$H21),theta2[colnames(df_ls$H21.A2)])>0)
    df_ls$opt.SSL <- data.frame(Q1opt,Q2opt_)
  }
  return(df_ls)
}

VsslDR <- function(U_ls,theta1,theta2){
  Q1opt <- U_ls$opt.SSL$Q1opt
  Q2opt_ <- U_ls$opt.SSL$Q2opt_
  VsslDR.hat <- Q1opt + U_ls$omegas$omega1*((1+theta2['Y1'])*U_ls$Ys$mu1.V-Q1opt+Q2opt_) + 
                U_ls$omegas$mu2_w2-theta2['Y1']*U_ls$omegas$mu1_w2-Q2opt_*U_ls$omegas$mu_w2
  U_ls$opt.SSL$VsslDR.hat <- VsslDR.hat#[which(U_ls$pi$pi1>.1 & U_ls$pi$pi1<.9 )]####
  
  return(U_ls)
}

################################################################
################## Random Forests & Splines ####################
################################################################

library(randomForest)

random_forests <- function(L_ls,U_ls,Q_fun){
  # Check if outcome is binary and prepare data accordingly
  if(length(unique(L_ls$Ys[,'Y1']))<=3) {# check if outcome is binary (0,1 & NA values)
    L_tmp <- list(); L_tmp$Y1.sq <- as.factor(L_ls$Ys[,'Y1']^2); L_tmp$Y1xY2 <- as.factor(L_ls$Ys[,'Y1']*L_ls$Ys[,'Y2'])
    cont.Y <- F; L_ls$Ys[,'Y1'] <- as.factor(L_ls$Ys[,'Y1']); L_ls$Ys[,'Y2'] <- as.factor(L_ls$Ys[,'Y2'])
  }else{
    cont.Y <- T
    L_tmp <- list(); L_tmp$Y1.sq <- L_ls$Ys[,'Y1']^2; L_tmp$Y1xY2 <- L_ls$Ys[,'Y1']*L_ls$Ys[,'Y2']
    }
  if(Q_fun){ # imputing for Q function
    # Y1
    RF.fit1 <- randomForest(Y1 ~ . , data = cbind(L_ls$Uvec,Y1=L_ls$Ys[,'Y1']) , subset = !is.na(L_ls$Ys[,'Y1']),ntree=500)
    if(cont.Y){
      RF.preds1 <- predict(RF.fit1,U_ls$Uvec)
    }else{
      RF.preds1 <- predict(RF.fit1,U_ls$Uvec, type = "prob")[,'1']
    }
    U_ls$Ys <- data.frame(m1=RF.preds1)

    #Y2
    RF.fit2 <- randomForest(Y2 ~ . , data = cbind(L_ls$Uvec,Y2=L_ls$Ys[,'Y2']) , subset = !is.na(L_ls$Ys[,'Y2']),ntree=500)
    if(cont.Y){
      RF.preds2 <- predict(RF.fit2,U_ls$Uvec)
    }else{
      RF.preds2 <- predict(RF.fit2,U_ls$Uvec, type = "prob")[,'1']
    }
    #mse <- with(unlabeled_data[,c(available_covs,'Y2')], mean( (Y2 - RF.preds2)^2)) 
    #mse
    U_ls$Ys['m2'] <- RF.preds2
    
    # Y1^2
    RF.fitY1.sq <- randomForest(Y1.sq ~ . , data = cbind(L_ls$Uvec,Y1.sq=L_tmp$Y1.sq) , subset = !is.na(L_ls$Ys[,'Y1']),ntree=500)
    if(cont.Y){
      RF.predsY1.sq <- predict(RF.fitY1.sq,U_ls$Uvec)
    }else{
      RF.predsY1.sq <- predict(RF.fitY1.sq,U_ls$Uvec,type = "prob")[,'1']
    }
    U_ls$Ys['m11'] <- RF.predsY1.sq
    
    # Y1xY2
    RF.fitY1xY2 <- randomForest(Y1xY2 ~ . , data = cbind(L_ls$Uvec,Y1xY2=L_tmp$Y1xY2) , subset = !is.na(L_tmp$Y1xY2),ntree=500)
    if(cont.Y){
      RF.predsY1xY2 <- predict(RF.fitY1xY2,U_ls$Uvec)
    }else{
      RF.predsY1xY2 <- predict(RF.fitY1xY2,U_ls$Uvec,type = "prob")[,'1']
    }
    U_ls$Ys['m12'] <- RF.predsY1xY2
    
  }else{ # imputing for value function
    # w2, Y1w2, Y2w2
    RF.fit_A2 <- randomForest(pi2 ~ . , data = cbind(L_ls$Uvec,pi2=L_ls$pi$pi2), ntree=500)
    mu.pi2 <- predict(RF.fit_A2, U_ls$Uvec)
    A2 <- U_ls$H21.A2[,'A2']; d2 <- U_ls$D.hat.SSL$d2
    U_ls$omegas$m_w2 <- d2*A2/mu.pi2+(1-d2)*(1-A2)/(1-mu.pi2)
    U_ls$omegas$m1_w2<- U_ls$Ys$m1*U_ls$omegas$m_w2
    U_ls$omegas$m2_w2<- U_ls$Ys$m2*U_ls$omegas$m_w2
  }
  return(U_ls)
}

library(splines)

basis_expansion <- function(L_ls,U_ls,Q_fun){
  # Check if outcome is binary and prepare data accordingly
  if(length(unique(L_ls$Ys[,'Y1']))<=3) {# check if outcome is binary (0,1 & NA values)
    L_tmp <- list(); L_tmp$Y1.sq <- as.factor(L_ls$Ys[,'Y1']^2); L_tmp$Y1xY2 <- as.factor(L_ls$Ys[,'Y1']*L_ls$Ys[,'Y2'])
    cont.Y <- F; L_ls$Ys[,'Y1'] <- as.factor(L_ls$Ys[,'Y1']); L_ls$Ys[,'Y2'] <- as.factor(L_ls$Ys[,'Y2'])
  }else{
    cont.Y <- T
    L_tmp <- list(); L_tmp$Y1.sq <- L_ls$Ys[,'Y1']^2; L_tmp$Y1xY2 <- L_ls$Ys[,'Y1']*L_ls$Ys[,'Y2']
  }
  # natural cubic splines knots number
  knots <- 2
  L.main.effs <- L_ls$Uvec[,-grep('x',colnames(L_ls$Uvec))]
  U.main.effs <- U_ls$Uvec[,-grep('x',colnames(U_ls$Uvec))]
  # Make splines matrix or return categorial column
  L.splines_ls <- apply(L.main.effs,2,function(x){tryCatch(ns(x,df = knots+1),error=function(e){as.numeric(x)})})
  # evaluate the basis at the unlabelede data
  splns.nms <- names(L.splines_ls)
  U.splines_ls <- lapply(splns.nms,function(nm){tryCatch(predict(L.splines_ls[[nm]],U.main.effs[,nm]),error=function(e){as.numeric(U.main.effs[,nm])})})
  # turn splines list into matrix
  L.splines_mat <- do.call("cbind", L.splines_ls)
  U.splines_mat <- do.call("cbind", U.splines_ls)
  colnames(U.splines_mat) <- colnames(L.splines_mat)
  if(Q_fun){# imputing for Q function
    # Y1
    if(cont.Y){
      #BE.fit1 <- lm(Y1~.,cbind(L.splines_mat,Y1=L_ls$Ys[,'Y1']))
      BE.fit1 <- cv.glmnet(x=L.splines_mat,y=L_ls$Ys[,'Y1'],family = 'gaussian',alpha=0)
    }else{
    #  BE.fit1 <- glm(Y1~.,cbind(L.splines_mat,Y1=L_ls$Ys[,'Y1']),family = 'binomial')
      BE.fit1 <- cv.glmnet(x=L.splines_mat,y=L_ls$Ys[,'Y1'],family = 'binomial',alpha=0)
    }
    BE.preds1 <- predict(BE.fit1, U.splines_mat,type='response',s = "lambda.1se")
    U_ls$Ys <- data.frame(m1=as.numeric(BE.preds1))
    
    #Y2
    if(cont.Y){
      BE.fit2 <- cv.glmnet(x=L.splines_mat,y=L_ls$Ys[,'Y2'],family = 'gaussian',alpha=0)
    }else{
      BE.fit2 <- cv.glmnet(x=L.splines_mat,y=L_ls$Ys[,'Y2'],family = 'binomial',alpha=0)
    }
    BE.preds2 <- predict(BE.fit2, U.splines_mat,type='response',s = "lambda.1se")
    U_ls$Ys['m2'] <- as.numeric(BE.preds2)
    
    
    #  unlab_dat$imputed_Y1.sq <- as.numeric(unlab_dat$imputed_Y1^2 + 
     #                                         exp(t(tcrossprod(t(alpha_1.imp),as.matrix(cbind(1,unlab_dat[,covs]))))))
    # Y1^2
    RF.fitY1.sq <- randomForest(Y1.sq ~ . , data = cbind(L_ls$Uvec,Y1.sq=L_tmp$Y1.sq) , subset = !is.na(L_ls$Ys[,'Y1']),ntree=500)
    if(cont.Y){
      BE.fit11 <- cv.glmnet(x=L.splines_mat,y=L_tmp$Y1.sq,family = 'gaussian',alpha=0)
    }else{
      BE.fit11 <- cv.glmnet(x=L.splines_mat,y=L_tmp$Y1.sq,family = 'binomial',alpha=0)
    }
    BE.preds11 <- predict(BE.fit11, U.splines_mat,type='response',s = "lambda.1se")
    U_ls$Ys['m11'] <- as.numeric(BE.preds11)
    
    # Y1xY2
    if(cont.Y){
      BE.fit12 <- cv.glmnet(x=L.splines_mat,y=L_tmp$Y1xY2,family = 'gaussian',alpha=0)
    }else{
      BE.fit12 <- cv.glmnet(x=L.splines_mat,y=L_tmp$Y1xY2,family = 'binomial',alpha=0)
    }
    BE.preds12 <- predict(BE.fit12, U.splines_mat,type='response',s = "lambda.1se")
    U_ls$Ys['m12'] <- as.numeric(BE.preds12)
    
  }else{ # imputing for value function
    # w2, Y1w2, Y2w2
    #if(sim.setting == 'continuous'){
      BE.fit_A2 <- cv.glmnet(y=L_ls$H21.A2[,'A2'],x=L.splines_mat,family = 'binomial',alpha=.1)
    #}else{
    #  BE.fit_A2 <- cv.glmnet(x=L.splines_mat,y=L_ls$pi$pi2,family = 'gaussian',alpha=0.1)
    #}
    mu.pi2 <- predict(BE.fit_A2, U.splines_mat,type='response',s = "lambda.1se")
    A2 <- U_ls$H21.A2[,'A2']; d2 <- U_ls$D.hat.SSL$d2
    U_ls$omegas$m_w2 <- d2*A2/mu.pi2+(1-d2)*(1-A2)/(1-mu.pi2)
    U_ls$omegas$m1_w2<- U_ls$Ys$m1*U_ls$omegas$m_w2
    U_ls$omegas$m2_w2<- U_ls$Ys$m2*U_ls$omegas$m_w2
  }
  return(U_ls)
}

###################################################################
######################### Imputing function #######################
###################################################################

impute <- function(L_ls,U_ls,model,CV_grps,Q_fun){
  
  k <- length(CV_grps)
  imputed_sets <- data.frame(matrix(data=NA,nrow=0,ncol = 4))
  # Sample split imputation
  for (kk in c(1:k)){
    # Define the labeled index sets L and L-k to impute with:
    L._k <- if(length(CV_grps)>1){do.call(c, CV_grps[c(1:k)[-kk]])}else{1:nrow(L_ls$Ys)}#else{CV_grps[[1]]}#
    # Impute using RF, NN or cubic splines
    if(model == 'RF'){
      tmp <- random_forests(L_ls=lapply(L_ls, function(x) x[L._k,]),U_ls,Q_fun)
    }else if(model == 'NN'){
      imputed_sets <- rbind(imputed_sets,cbind(group=kk,id=c(1:nrow(unlab_dat)),neural_network(lab_dat[L._k,],unlab_dat,available_covs,samp_forY1,samp_forY2,impute_4Ys,Q_fun)))
    }else if(model == 'BE'){
      tmp <- basis_expansion(L_ls=lapply(L_ls, function(x) x[L._k,]),U_ls,Q_fun)
    }
    if(Q_fun){
      imputed_sets <- rbind(imputed_sets,cbind(group=kk,id=1:nrow(tmp$Ys),tmp$Ys))
    }else{
      imputed_sets <- rbind(imputed_sets,cbind(group=kk,id=1:nrow(tmp$Ys),tmp$omegas))
    }
  }
  
  # Compute the mean of the imputations 
  
  if(Q_fun){ # If not debiasing for value function
    imps <- imputed_sets %>% group_by(id) %>% dplyr::summarise(m1 = mean(m1),
                                                               m11 = mean(m11),
                                                               m2 = mean(m2),
                                                               m12 = mean(m12), 
                                                               .groups = 'drop') %>% data.frame() 
    U_ls$Ys <- imps[,-1]
  }else{
    imps <- imputed_sets %>% group_by(id) %>% dplyr::summarise(m_w2 = mean(m_w2),
                                                               m2_w2 = mean(m2_w2),
                                                               m1_w2 = mean(m1_w2), 
                                                               .groups = 'drop') %>% data.frame()
    U_ls$omegas <- cbind(U_ls$omegas,imps[,-1])
  }
  # remove grouping variable
  return(U_ls)
}

###################################################################
########################## Bias correction ########################
###################################################################

refitting <- function(L_ls,U_ls,model,CV_grps,Q_fun){
  # Variables for regression on X
  X.vars <- unique(colnames(cbind(L_ls$H10,L_ls$H11.A1,L_ls$H20,L_ls$H21.A2))); X.vars <- X.vars[-grep('(Intercept)',X.vars)]
  impute.L_data <- data.frame(matrix(data=NA,nrow=0,ncol = 2*length(L_ls$Ys)))
  # Impute labeled data with labeled data using CV to estimate the errors: Y-m(x)
  for (kk in c(1:length(CV_grps))){
    # Define the label index sets L and L-k to impute with:
    L.k <- CV_grps[[kk]]
    L._k <- if(length(CV_grps)>1){do.call(c, CV_grps[c(1:length(CV_grps))[-kk]])}else{L.k}
    # Impute the k^th set using the rest of the labeled data
    tmp <- impute(L_ls=lapply(L_ls, function(x) x[L._k,]),U_ls=lapply(L_ls, function(x) x[L.k,]),
                  model=impute_with,CV_grps=list(`1`=c(1:length(L._k))),Q_fun)
    if(Q_fun){
      impute.L_data <- rbind(impute.L_data,cbind(group=kk,id=c(1:length(L.k)),tmp$Ys,cbind(L_ls$Ys,L_ls$H10,L_ls$H11.A1,L_ls$H20,L_ls$H21.A2)[L.k,c(colnames(L_ls$Ys),X.vars)]))
    }else{
      impute.L_data <- rbind(impute.L_data,cbind(group=kk,id=c(1:length(L.k)),tmp$Ys,tmp$opt.SSL,tmp$pi,tmp$D.hat.SSL,tmp$omegas,cbind(L_ls$Ys,L_ls$H10,L_ls$H11.A1,L_ls$H20,L_ls$H21.A2)[L.k,c(colnames(L_ls$Ys),X.vars)]))
    }
  }
  
  # Calculate the residuals for the k^thg group
  if(Q_fun){
    impute.L_data <- impute.L_data %>% mutate(m1_residuals=Y1-m1,
                                              m2_residuals=Y2-m2,
                                              m11_residuals=Y1^2-m11,
                                              m12_residuals=Y1*Y2-m12)
  }else{ # value function residuals:
    impute.L_data <- impute.L_data[,unique(colnames(impute.L_data))]
    impute.L_data <- impute.L_data %>% mutate(ipw2 = d2*A2/pi2+(1-d2)*(1-A2)/(1-pi2),
                                              Y1_residuals.V=omega1*(Y1-m1)/mean(omega1),
                                              w2_residuals=Q2opt_*omega1*(ipw2-m_w2)/mean(Q2opt_),
                                              w2.Y1_residuals=omega1*(Y1*ipw2-m1_w2),
                                              w2.Y2_residuals=omega1*(Y2*ipw2-m2_w2))
  }
  
  
  # Bias correction for unlabeled data (using all errors from labeled data to train eta):
    X.L <- data.frame(impute.L_data[,X.vars]); colnames(X.L) <- X.vars
    X.U <- data.frame(cbind(U_ls$Ys,U_ls$H10,U_ls$H11.A1,U_ls$H20,U_ls$H21.A2)[,X.vars]); colnames(X.U) <- X.vars
    
  # Estimate eta by regressing the residuals of the imputation result to X=(x1,x2)
    if(Q_fun){# Q function debias:
      eta <- list(m1_residuals=lm(impute.L_data$m1_residuals ~ .,data=X.L), #regression on X
                  m12_residuals=lm(m12_residuals ~ 1,data=impute.L_data), 
                  m11_residuals=lm(m11_residuals ~1,data=impute.L_data), 
                  m2_residuals=lm(impute.L_data$m2_residuals ~ .,data=X.L)) #regression on X
      # Compute mu(x)=m(x)+x*eta for unlabeled data:
      # correct names of coefficients
      U_ls$Ys$mu1 <- U_ls$Ys$m1 + predict(eta$m1_residuals,X.U)
      U_ls$Ys$mu12 <- U_ls$Ys$m12 + predict(eta$m12_residuals,X.U)
      U_ls$Ys$mu11 <- U_ls$Ys$m11 + predict(eta$m11_residuals,X.U)  
      U_ls$Ys$mu2 <- U_ls$Ys$m2 + predict(eta$m2_residuals,X.U)
      
    }else{# V function debias
      
      eta <- list(Y1_residuals.V=lm(impute.L_data$Y1_residuals.V ~ 1,data=X.L),
                  w2_residuals=lm(impute.L_data$w2_residuals ~ A2,data=X.L),
                  w2.Y1_residuals=lm(impute.L_data$w2.Y1_residuals ~ A2,data=X.L),
                  w2.Y2_residuals=lm(impute.L_data$w2.Y2_residuals ~ A2,data=X.L))
      
      # Compute mu(x)=m(x)+x*eta for unlabeled data:
      U_ls$Ys$mu1.V <- U_ls$Ys$m1 + predict(eta$Y1_residuals.V,X.U)
      U_ls$omegas$mu_w2 <- U_ls$omegas$m_w2 + predict(eta$w2_residuals,X.U)
      U_ls$omegas$mu1_w2 <- U_ls$omegas$m1_w2 + predict(eta$w2.Y1_residuals,X.U)
      U_ls$omegas$mu2_w2 <- U_ls$omegas$m2_w2 + predict(eta$w2.Y2_residuals,X.U)
    }
    return(U_ls)
}

# impute labeled set (to get errors for influence function):
reffiting4L <- function(L_ls,model,CV_grps,Q_fun){
  L.set <- Lset.4V <- list()
  tmp <- c()
  for (kk in c(1:length(CV_grps))){
    # subset data
    L.k <- CV_grps[[kk]]; 
    L._k <- if(length(CV_grps)>1){do.call(c, CV_grps[c(1:length(CV_grps))[-kk]])}else{L.k}
    L_ls._k <- lapply(L_ls, function(x) x[L._k,]) # all but indices in I_k
    L_ls.k <- lapply(L_ls, function(x) x[L.k,]) # indices in I_k
    # Step 1) Imputation
    L_ls.k <- impute(L_ls._k,L_ls.k,model,CV_grps=list(1:length(L._k)),Q_fun)
    # Step 2) Refitting step
    L_ls.k <- refitting(L_ls=L_ls._k,U_ls=L_ls.k,model,CV_grps=list(1:length(L._k)),Q_fun)
    L.set[[kk]] <- L_ls.k$Ys
    if(!Q_fun) Lset.4V[[kk]] <- L_ls.k$omegas 
  }
  # reorder the labeled set to be congruent with the imputations
  L_ls <- lapply(L_ls, function(x) x[do.call(c, CV_grps),])
  L_ls$Ys <- cbind(L_ls$Ys,do.call(rbind, L.set))
  if(!Q_fun) L_ls$omegas <- cbind(L_ls$omegas,do.call(rbind, Lset.4V))
  return(L_ls)
}

###################################################################
####################### Estimate Sigma (IF)  ######################
###################################################################

unlab_SE <- function(theta,l.th1,L_ls){
  
#  df_n <- df_n[!is.na(df_n$Y1) & !is.na(df_n$Y2),]
#  df_N <- dat.ls$dat_H2
  n <- nrow(L_ls$Ys)
  IF_theta <- matrix(NA,length(theta),n)
  ### Theta 2 ###
  theta_2 <- theta[(l.th1+1):length(theta)];l.th2 <- length(theta_2)
  # extract main and interaction term coefficients
  beta_21 <- theta_2[1]; beta_22 <- theta_2[-grep('A2',names(theta_2))][-1]
  gamma_2 <- theta_2[grep('A2',names(theta_2))]
  # Compute bar q_2 portion in the IF:
#  df_n[,"(Intercept)"] <- 1; df_N[,"(Intercept)"] <- 1
  q2 <- mean(U_ls$H20%*%beta_22+U_ls$H21.A2%*%gamma_2)
  # first row of IF2:
  psi_21i <- with(L_ls$Ys,Y1*Y2-mu12-beta_21*(Y1^2-mu11)-q2*(Y1-mu1))
  # Compute the lower part of the IF2: 
#  H2A <- as.matrix(df_n[,names(c(beta_22,gamma_2))]) # (H20,A2*H21)
  psi2_err_Y2 <- with(L_ls$Ys,Y2-mu2)
  psi2_err_Y1 <- with(L_ls$Ys,Y1-mu1)
  H2Axerr_Y2 <- t(t(cbind(L_ls$H20,L_ls$H21.A2))%*%diag(psi2_err_Y2)) # (H20,A2*H21)(Y2-mu2)
  H2Axerr_Y1 <- t(t(cbind(L_ls$H20,L_ls$H21.A2))%*%diag(psi2_err_Y1)) # (H20,A2*H21)(Y1-mu1)
  psi_22i <- H2Axerr_Y2 - beta_21*H2Axerr_Y1
  # stack IF2
  IF_2i <- cbind(Y1=psi_21i,psi_22i)
  # Compute E[X_2X_2t]:
  X2 <- cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21.A2)
  E_xxT.inv2 <- solve(t(X2)%*%X2/n) ##****
  #E_xxT.inv2 <- solve(t(X2)%*%X2/N) ##****
  rownames(E_xxT.inv2)[1] <- colnames(E_xxT.inv2)[1] <- 'Y1'
  # Varcov2 is 1/n sum IF2%*%IF2t
  Sigma2 <- matrix(0,l.th2,l.th2)
  for (i in c(1:n)){
    IF_theta[1:l.th2,i] <- E_xxT.inv2%*%IF_2i[i,]
    Sigma2 = Sigma2 + E_xxT.inv2%*%IF_2i[i,]%*%t(E_xxT.inv2%*%IF_2i[i,])
  }
  Sigma2 <- Sigma2/n ##****
  SE_theta2 <- sqrt(diag(Sigma2)/n)  ##****
  #plot(sqrt(diag(results.OLS@vcov))[1:length(SE_theta2)],SE_theta2)
  #abline(c(0,1))
  
  ###############
  ### Theta 1 ###
  ###############
  theta_1 <- theta[1:l.th1]
  # Compute E[X_1X_1t]:
  X1 <- cbind(L_ls$H10,L_ls$H11.A1)
  #X1 <- as.matrix(df_N[,names(theta_1)])
  E_xxT.inv1 <- solve(t(X1)%*%X1/n)
  #E_xxT.inv1 <- solve(t(X1)%*%X1/N)
  # Compute E[X_1(Y1,H20t)]:
  E_X1_Y1H20t <- matrix(0,ncol(X1),length(c('Y1',beta_22)))
  for (i in c(1:n)){
    E_X1_Y1H20t <- E_X1_Y1H20t + as.matrix(X1[i,])%*%cbind(L_ls$Ys$Y1,L_ls$H20)[i,]
    
  }
  E_X1_Y1H20t <- E_X1_Y1H20t/n
  # Compute E[X_1H21t|H21tgamma_2>0]:
  E_X1_H21t.cond <- matrix(0,ncol(X1),length(gamma_2))
  trt2.indx <- L_ls$H21.A2%*%gamma_2>0
  for (i in which(trt2.indx)){
    E_X1_H21t.cond <- E_X1_H21t.cond + as.matrix(X1[i,])%*%L_ls$H21.A2[i,]
  }
  E_X1_H21t.cond <- E_X1_H21t.cond/sum(trt2.indx)
  P.trt2 <- sum(trt2.indx)/n ##--
  cumm_Sigma1 <- matrix(0,length(theta_1),length(theta_1))
  for (i in c(1:n)){
    currIF_2i <- E_xxT.inv2%*%IF_2i[i,]
    IF_1i <- 
      X1[i,]*(1+beta_21)*(L_ls$Ys$Y1[i]-L_ls$Ys$mu1[i])+
      E_X1_Y1H20t%*%currIF_2i[c(names(beta_21),names(beta_22)),]+
      P.trt2*E_X1_H21t.cond%*%currIF_2i[names(gamma_2),]
    cumm_Sigma1 <- cumm_Sigma1 + E_xxT.inv1%*%IF_1i%*%t(E_xxT.inv1%*%IF_1i)
    IF_theta[(l.th2+1):length(theta),i] <- E_xxT.inv1%*%IF_1i
  }
  Sigma1 <- cumm_Sigma1/n
  SE_theta1 <- sqrt(diag(Sigma1)/n)
  #plot(sqrt(diag(results.OLS@vcov))[(length(SE_theta2)+1):length(theta)],SE_theta1)
  #abline(c(0,1))
  return(list(SE_theta=c(SE_theta1,SE_theta2),IF_theta=IF_theta))
}

lab_SE <- function(theta,l.th1,L_ls){
#  lab.indx <- with(dat_ls$all_cols,!is.na(Y1) & !is.na(Y2))
#  df <- dat_ls$all_cols[lab.indx,]
  n <- nrow(L_ls$Ys)
  IF_theta <- matrix(NA,length(theta),n)
  l.th2 <- length(theta) - l.th1
  ### Theta 2 ###
  theta_2 <- theta[(l.th1+1):length(theta)]
  # extract main and interaction term coefficients
  beta_21 <- theta_2[1]; beta_22 <- theta_2[-grep('A2',names(theta_2))][-1]
  gamma_2 <- theta_2[grep('A2',names(theta_2))]
  # Compute bar q_2 portion in the IF:
#  df[,"(Intercept)"] <- 1
  # IF2:
  X2 <- cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21.A2); colnames(X2)[1] <- 'Y1'
  IF2_err <- L_ls$Ys$Y2-X2%*%theta_2
  # Compute E[X_2X_2t]:
  E_xxT.inv2 <- solve(t(X2)%*%X2/n) ##****
  # Varcov2 is 1/n sum IF2%*%IF2t
  Sigma2 <- matrix(0,l.th2,l.th2)
  for (i in 1:n){
    IF_2i <- E_xxT.inv2%*%X2[i,]*IF2_err[i]
    IF_theta[1:l.th2,i] <- IF_2i
    Sigma2 <- Sigma2 + IF_2i%*%t(IF_2i)
  }
  Sigma2 <- Sigma2/n 
  SE_theta2 <- sqrt(diag(Sigma2)/n) 
  #plot(sqrt(diag(results.OLS@vcov))[1:length(SE_theta2)],SE_theta2)
  #abline(c(0,1))
  #plot(bb$coefficients[,2],SE_theta2)
  #plot(bb$coefficients[,2],sqrt(diag(results.OLS@vcov))[1:46])
  
  ###############
  ### Theta 1 ###
  ###############
  theta_1 <- theta[1:l.th1]
  # Compute E[X_1X_1t]:
  X1 <- cbind(L_ls$H10,L_ls$H11.A1)
  E_xxT.inv1 <- solve(t(X1)%*%X1/n)
  # Compute E[X_1(Y1,H20t)]:
  E_X1_Y1H20t <- matrix(0,ncol(X1),length(c('Y1',beta_22)))
  for (i in c(1:n)){
    E_X1_Y1H20t <- E_X1_Y1H20t + as.matrix(X1[i,])%*%cbind(L_ls$Ys$Y1,L_ls$H20)[i,]
  }
  E_X1_Y1H20t <- E_X1_Y1H20t/n
  # Compute E[X_1H21t|H21tgamma_2>0]:
  E_X1_H21t.cond <- matrix(0,ncol(X1),length(gamma_2))
  trt2.indx <- L_ls$H21.A2%*%gamma_2>0
  for (i in which(trt2.indx)){
    E_X1_H21t.cond <- E_X1_H21t.cond + as.matrix(X1[i,])%*%L_ls$H21.A2[i,]
  }
  E_X1_H21t.cond <- E_X1_H21t.cond/sum(trt2.indx)
  P.trt2 <- sum(trt2.indx)/n
  
  maxY2 <- apply(rbind(tcrossprod(theta_2,cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21)),
                       tcrossprod(theta_2,cbind(L_ls$Ys$Y1,L_ls$H20,0*L_ls$H21))),2,max)
  
  Sigma1 <- matrix(0,length(theta_1),length(theta_1))
  for (i in c(1:n)){
    Y_1i.hat <- as.numeric(theta_1%*%cbind(L_ls$H10,L_ls$H11.A1)[i,])
    psY1.hat <- L_ls$Ys$Y1[i]+maxY2[i]
    IF_2i <- E_xxT.inv2%*%X2[i,]*IF2_err[i]
    IF_1i <- X1[i,]*(psY1.hat-Y_1i.hat)+
      E_X1_Y1H20t%*%IF_2i[c(names(beta_21),names(beta_22)),]+
      P.trt2*E_X1_H21t.cond%*%IF_2i[names(gamma_2),]
    IF_theta[(l.th2+1):length(theta),i] <- E_xxT.inv1%*%IF_1i
    Sigma1 <- Sigma1 + E_xxT.inv1%*%IF_1i%*%t(E_xxT.inv1%*%IF_1i)
  }
  Sigma1 <- Sigma1/n
  SE_theta1 <- sqrt(diag(Sigma1)/n)
  #plot(sqrt(diag(results.OLS@vcov))[(length(SE_theta2)+1):length(theta)],SE_theta1)
  #abline(c(0,1))
  return(list(SE_theta=c(SE_theta1,SE_theta2),IF_theta=IF_theta))
}

###################################################################
#################### Value function estimation ####################
###################################################################

IF.derivs <- function(L_ls,theta1,theta2,label=T){  
  A1 <- L_ls$H11.A1[,'A1']; A2 <- L_ls$H21.A2[,'A2']
  pi1 <- L_ls$pi$pi1; pi2 <- L_ls$pi$pi2
  omega1 <- L_ls$omegas$omega1; omega2 <- L_ls$omegas$omega2
  if(label){
    d1 <- L_ls$D.hat$d1; d2 <- L_ls$D.hat$d2
    Q1opt <- L_ls$opt.SUP$Q1opt; Q2opt <- L_ls$opt.SUP$Q2opt
  }else{
    Y1 <- L_ls$Ys$Y1; Y2 <- L_ls$Ys$Y2
    mu1.V <- L_ls$Ys$mu1.V; mu2_w2 <- L_ls$omegas$mu2_w2; mu1_w2 <- L_ls$omegas$mu1_w2; mu_w2 <- L_ls$omegas$mu_w2 
    d1 <- L_ls$D.hat.SSL$d1; d2 <- L_ls$D.hat.SSL$d2
    Q1opt <- L_ls$opt.SSL$Q1opt; Q2opt <- L_ls$opt.SSL$Q2opt
  }
  # Derivative of Vsup with respect to theta = (theta1,theta2)
  dQ1opt.dtheta <- cbind(L_ls$H10,d1*L_ls$H11,matrix(0,nrow=n,ncol=ncol(cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21))))
  dQ2opt.dtheta <- cbind(matrix(0,nrow=n,ncol=ncol(cbind(L_ls$H10,L_ls$H11))),L_ls$Ys$Y1,L_ls$H20,d2*L_ls$H21)
  
  # Derivative of the smooth versions of omega1, omega2:
  d1.smooth <- expit(crossprod(t(L_ls$H11),theta1[colnames(L_ls$H11.A1)]))
  d2.smooth <- expit(crossprod(t(L_ls$H21),theta2[colnames(L_ls$H21.A2)]))
  
  d1.smooth.deriv <- d1.smooth*(1-d1.smooth)
  d2.smooth.deriv <- d2.smooth*(1-d2.smooth)
  
  A1.deriv.prod <- (A1/pi1-(1-A1)/(1-pi1))*d1.smooth.deriv
  A2.deriv.prod.omega1 <- (A2/pi2-(1-A2)/(1-pi2))*d2.smooth.deriv*omega1
  
  domega1.dtheta <- cbind(matrix(0,nrow=n,ncol=ncol(L_ls$H10)),as.data.frame(L_ls$H11)*A1.deriv.prod,matrix(0,nrow=n,ncol=ncol(cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21))))
  domega2.dtheta <- domega1.dtheta*(d2.smooth*A2/pi2+(1-d2.smooth)*(1-A2)/(1-pi2)) + #domega1.dtheta*(d2*A2/pi2+(1-d2)*(1-A2)/(1-pi2)) + 
    cbind(matrix(0,nrow=n,ncol=ncol(cbind(L_ls$H10,L_ls$H11,L_ls$Ys$Y1,L_ls$H20))),as.data.frame(L_ls$H21)*A2.deriv.prod.omega1)
  
  Vsup.dtheta <- dQ1opt.dtheta + domega1.dtheta*(L_ls$Ys$Y1-Q1opt+Q2opt)+L_ls$omegas$omega1*(-dQ1opt.dtheta+dQ2opt.dtheta) +
    domega2.dtheta*(L_ls$Ys$Y2-Q2opt)-L_ls$omegas$omega2*data.frame(dQ2opt.dtheta)
  
  Vsup.dtheta <- apply(Vsup.dtheta,2,mean)
  
  # Derivative of Vsup with respect to xi = (xi1, xi2)
  
  d1.deriv.prod <- -d1*A1*(1-pi1)/pi1+(1-d1)*(1-A1)*pi1/(1-pi1)
  d2.deriv.prod <- (-d2*A2*(1-pi2)/pi2+(1-d2)*(1-A2)*pi2/(1-pi2))*L_ls$omegas$omega1
  
  domega1.dxi <- cbind(L_ls$H1check*d1.deriv.prod,matrix(0,nrow=n,ncol=ncol(L_ls$H2check)))
  
  domega2.dxi <- cbind(matrix(0,nrow=n,ncol=ncol(L_ls$H1check)),L_ls$H2check*d2.deriv.prod) 
  Vsup.dxi <- as.data.frame(domega1.dxi)*(L_ls$Ys$Y1-Q1opt+Q2opt)+
              as.data.frame(domega2.dxi)*(L_ls$Ys$Y2-Q2opt)
  
  Vsup.dxi <- as.matrix(apply(Vsup.dxi,2,mean))
  # Influence function for the propensity score fitting:
  ## xi1
  IF_xi1 <- mean(crossprod(crossprod(t(L_ls$H1check)),as.matrix(pi1*(1-pi1))))^(-1)*L_ls$H1check*(A1-pi1)
  ## xi2
  IF_xi2 <- mean(crossprod(crossprod(t(L_ls$H2check)),as.matrix(pi2*(1-pi2))))^(-1)*L_ls$H2check*(A2-pi2)
  
  IF_xi <- as.matrix(cbind(IF_xi1,IF_xi2))
  
  if(label){
    L_ls$opt.SUP$V.psi <- L_ls$opt.SUP$VsupDR.hat-mean(L_ls$opt.SUP$VsupDR.hat) + t(tcrossprod(Vsup.dtheta,IF_theta.SUP)) + IF_xi%*%Vsup.dxi
  }else{
    L_ls$opt.SSL$V.psi <- omega1*(1+theta2['Y1'])*(Y1-mu1.V)+omega2*Y2-mu2_w2-theta2['Y1']*(omega2*Y1-mu1_w2)-L_ls$opt.SSL$Q2opt_*(omega2-mu_w2) +
                          t(tcrossprod(Vsup.dtheta,IF_theta.SSL)) + IF_xi%*%Vsup.dxi
    
    #V.psi_ssl <- with(labeled_data,omega1*(1+Q2.est_np[H20.nms[1]])*(Y1-Y1_V)+omega2*Y2-Y2.w2_V-Q2.est_np[H20.nms[1]]*(omega2*Y1-Y1.w2_V)-Q2opt_*(omega2-w2.inv_V)) +
    #  t(tcrossprod(Vsup.dtheta,IF_theta.SS)) + IF_xi%*%Vsup.dxi
  }
  return(L_ls)
}

###################################################################
###################### Estimating Equations #######################
###################################################################

library(rootSolve)

supQl.EE <- function(theta,L_ls){
  
  lth1 <- ncol(cbind(L_ls$H10,L_ls$H11))
  beta2 <- t(as.matrix(theta[(lth1+1):(lth1+ncol(L_ls$H20)+1)]))
  gamma2 <- t(as.matrix(theta[(lth1+2+ncol(L_ls$H20)):length(theta)]))
  theta1 <- t(as.matrix(theta[1:lth1])); theta2 <- t(as.matrix(c(beta2,gamma2)))
  # second stage EE
  X2 <- cbind(L_ls$Ys$Y1,L_ls$H20,L_ls$H21.A2)
  Y2_Xth <- L_ls$Ys$Y2-tcrossprod(X2,theta2)
  ssEE <- t(X2)%*%Y2_Xth/n
  # first stage EE
  X1 <- cbind(L_ls$H10,L_ls$H11.A1)
  H20.b2 <- tcrossprod(cbind(L_ls$Ys$Y1,L_ls$H20),beta2)
  H21.g2 <- tcrossprod(L_ls$H21,gamma2)*(tcrossprod(L_ls$H21,gamma2)>0)
  fsEE <- t(X1)%*%(L_ls$Ys$Y1 + H20.b2 + H21.g2 - tcrossprod(X1,theta1))/n
  
  c(ssEE,fsEE)
}

#theta = 1:13
SSLQl.EE <- function(theta){
  
  lth1 <- ncol(cbind(U_ls$H10,U_ls$H11))
  beta2 <- t(as.matrix(theta[(lth1+1):(lth1+ncol(U_ls$H20)+1)]))
  gamma2 <- t(as.matrix(theta[(lth1+2+ncol(U_ls$H20)):length(theta)]))
  theta1 <- t(as.matrix(theta[1:lth1])); theta2 <- t(as.matrix(c(beta2,gamma2)))
  # second stage EE
  X2 <- cbind(U_ls$H20,U_ls$H21.A2)
  augmX2a <- cbind(U_ls$Ys$mu11,U_ls$Ys$mu1*X2)
  augmX2b <- cbind(U_ls$Ys$mu1,X2)
  row1 <- mean(U_ls$Ys$mu12-tcrossprod(augmX2a,theta2))
  rwos2 <- t(X2)%*%(U_ls$Ys$mu2-tcrossprod(augmX2b,theta2))/N
  ssEE <- c(row1,rwos2)
  # first stage EE
  X1 <- cbind(U_ls$H10,U_ls$H11.A1)
  if (length(beta2[-1])>0){
    H20.b2 <- beta2[1]*U_ls$Ys$mu1 + U_ls$H20%*%beta2[-1]
  }else{
    H20.b2 <- beta2[1]*U_ls$Ys$mu1
  }
  H21.g2 <- tcrossprod(U_ls$H21,gamma2)*(tcrossprod(U_ls$H21,gamma2)>0)
  fsEE <- t(X1)%*%(U_ls$Ys$mu1 + H20.b2 + H21.g2 - tcrossprod(X1,theta1))/N
  
  c(ssEE,fsEE)
} 


###################################################################
######################## Wrapper function #########################
###################################################################

  
  
# i <- 0
# for (setting in c('binary', 'continuous' , 'EHR')){
#   for(samples in list(c(10000,500),c(1272,135))){#)){#
#     for (impute_with in c('RF','BE')){
#       for (beta26_phi24 in list(c(-1,0),c(0,-1),c(0,0),c(0,1),c(1,0))){
#         i=i+1
#         if(samples[2]==135){
#           cat(paste('sbatch run_sims.sh ',impute_with,setting, beta26_phi24[1], beta26_phi24[2], samples[1], samples[2] ,'\n\n\n',sep=' '))
#         }else{
#           cat(paste('sbatch run_sims_medium.sh ',impute_with,setting, beta26_phi24[1], beta26_phi24[2], samples[1], samples[2] ,'\n\n\n',sep=' '))
#           }
#       }
#     }
#   }
# }
# i
