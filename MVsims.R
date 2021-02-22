#setwd("~/Dropbox (HMS)/Documents/Research/RL/SS RL/SSRL_code/")
rm(list=ls())
library(dplyr)

write('args parameters', stderr())
args <- commandArgs(trailingOnly = TRUE)
impute_with <- as.character(args[2]); write(impute_with, stderr())
sim.setting <- as.character(args[4]); write(c('simulation setting: ',sim.setting), stderr())
beta26 <- as.numeric(args[6]); write(c('beta26: ',beta26), stderr())
phi24 <- as.numeric(args[8]); write(c('phi24: ',phi24), stderr())
N <- as.numeric(args[10]); write(c('N: ',N), stderr())
n <- as.numeric(args[12]); write(c('n: ',n), stderr())
if(is.na(impute_with)){
  print('local params')
  impute_with <- 'BE'#BE,RF
  sim.setting <- 'EHR'#'binary' 'continuous' or 'EHR'
  beta26 <- 1 #Q func. missp.
  phi24 <- 0 #prop. score missp. 
  # Labeled dataset size
  n <- 135#500#
  # Unlabeled dataset size
  N <-1272# 10000#
}

print_every <- 1#00

library(methods)
source('utils.R') 

#### Algorithm
# 1) Impute Y1 and Y2 using KS with labeled data set n
# 2) Run regular Q-learning method with imputed outcomes

if (sim.setting %in% c('binary','continuous')){
  # Q function working models:
  S1 <- ~ O1 + A1*O1; S2 <- Y2 ~ Y1 + O1*A1 + O2 +  (A1+O2)*A2# low dimension
  # Propensity score function working models:
  PS1 <- as.factor(A1) ~ O1; PS2 <- as.factor(A2) ~ O1+A1+O2+A1*O1+I(O2^2)+Y1# low dimension
}else if(sim.setting=='EHR'){
  S1 <- ~ O1+O2+O3+O4+O5+O6+A1*(O2+O3+O4+O5+O6); S2 <- Y2 ~ Y1+O1+O2+O3+O4+O5+O6+A1+Z.21+Z.21+Z.22+A2*(O1+O2+O3+O6+A1+Z.22+Z.21)# high dimension
  PS1 <- as.factor(A1) ~ O1+O2+O3; PS2 <- as.factor(A2) ~ Y1+O1+A1+Z.21+O2# high dimension
}
# parameters:
source('parameters.R')
# Theta estimation:
MC.theta <- MCparams(sim.setting,S1,S2,PS1,PS2,phi1,beta1,delta1,phi2,beta2,psi2)
V.star <- MC.theta$V.star
V.star
###
sims_No <- 1000


MSE_sim <- data.frame(theta1_SSL=rep(NA,sims_No),theta1_SUP=rep(NA,sims_No),
                      theta2_SSL=rep(NA,sims_No),theta2_SUP=rep(NA,sims_No))

theta2_SUP <- data.frame(matrix(NA,nrow = sims_No,ncol = length(MC.theta$theta2)))
colnames(theta2_SUP) <- names(MC.theta$theta2)
theta2_SE_SSL <- theta2_SE_SUP <- theta2_SSL <- theta2_SUP

theta1_SUP <- data.frame(matrix(NA,nrow = sims_No,ncol = length(MC.theta$theta1)))
colnames(theta1_SUP) <- names(MC.theta$theta1)
theta1_SE_SSL <- theta1_SE_SUP <- theta1_SSL <- theta1_SUP

V.hat <- data.frame(matrix(NA,nrow = sims_No,ncol = 2)); colnames(V.hat) <- c('SUP','SSL'); within_95CI_V <- V.hat.SE <- V.hat


within_95CI_SSL <- data.frame(matrix(NA,sims_No,length(MC.theta$theta1)+length(MC.theta$theta2))); colnames(within_95CI_SSL) <- names(c(MC.theta$theta1,MC.theta$theta2))
within_95CI_SUP <- within_95CI_SSL


for (sim in c(1:sims_No)){  
  # Generate data & get into suitable format:
  dat <- gen_data(sim.setting,size=n+N,phi1,beta1,delta1,phi2,beta2,psi2,seed=116687+sim) 
  df_ls <- transform_data(S1,S2,PS1,PS2,dat,n,N)
  
  U_ls <- df_ls$U_ls; L_ls <- df_ls$L_ls
  
  # groups for CV
  k <- 5         

  grps <- split(sample(c(1:n), replace = F), 1:k)
  
  # Step 1) Imputation
  U_ls <- impute(L_ls,U_ls,model=impute_with,CV_grps=grps,Q_fun=T)
  
  # Step 2) Refitting step
  U_ls <- refitting(L_ls,U_ls,model=impute_with,CV_grps=grps,Q_fun=T)
  
  # Step 3) Projection
  results.SSL <- multiroot(f = SSLQl.EE, start = rep(0,length(MC.theta$theta1)+length(MC.theta$theta2)))
  
  theta.SSL <- results.SSL$root
  names(theta.SSL) <- names(c(MC.theta$theta1,MC.theta$theta2))
  
  # compute standard errors:
  # impute labeled set to get estimates for the errors:
  L_ls <- reffiting4L(L_ls,model=impute_with,CV_grps=grps,Q_fun=T)
  ## Standard errors:
  SSL.res <- unlab_SE(theta=theta.SSL,l.th1=length(MC.theta$theta1),L_ls)
  SSL.SEs <- SSL.res$SE_theta
  IF_theta.SSL <- t(SSL.res$IF_theta)
  
  # Store results:
  
  Q2.est_SSL <- theta.SSL[(length(MC.theta$theta1)+1):length(theta.SSL)]
  MSE_sim[sim,'theta2_SSL'] <- crossprod(Q2.est_SSL-MC.theta$theta2)
  
  theta2_SSL[sim,] <- Q2.est_SSL
  theta2_SE_SSL[sim,] <- SSL.SEs[(length(MC.theta$theta1)+1):length(theta.SSL)]
  
  Q1.est_SSL <- theta.SSL[1:length(MC.theta$theta1)]
  MSE_sim[sim,'theta1_SSL'] <- crossprod(Q1.est_SSL-MC.theta$theta1)
  theta1_SSL[sim,] <- Q1.est_SSL
  theta1_SE_SSL[sim,] <- SSL.SEs[1:length(MC.theta$theta1)]
  
  
  # Labeled fitting:
  results.SUP <- multiroot(f = supQl.EE, start = rep(0,length(theta.SSL)),L_ls=L_ls)
  theta.SUP <- results.SUP$root
  names(theta.SUP) <- names(theta.SSL)
  
  ## Standard errors:
  
  SUP.res <- lab_SE(theta=theta.SUP,l.th1=length(MC.theta$theta1),L_ls)
  SUP.SEs <- SUP.res$SE_theta
  IF_theta.SUP <- t(SUP.res$IF_theta)
  
  # Store results
  Q2.est_SUP <- theta.SUP[(length(MC.theta$theta1)+1):length(theta.SUP)]
  MSE_sim[sim,'theta2_SUP'] <- crossprod(Q2.est_SUP-MC.theta$theta2)
  
  theta2_SUP[sim,] <- Q2.est_SUP
  theta2_SE_SUP[sim,] <- SUP.SEs[(length(MC.theta$theta1)+1):length(theta.SUP)]
  
  Q1.est_SUP <- theta.SUP[1:length(MC.theta$theta1)]
  
  MSE_sim[sim,'theta1_SUP'] <- crossprod(Q1.est_SUP-MC.theta$theta1)
  
  theta1_SUP[sim,] <- Q1.est_SUP
  
  theta1_SE_SUP[sim,] <- SUP.SEs[1:length(MC.theta$theta1)]
  
  #plot(x=c(apply(theta1_SE_SUP,2,mean,na.rm=T),apply(theta2_SE_SUP,2,mean,na.rm=T)),
  #     y=c(apply(theta1_SE_SSL,2,mean,na.rm=T),apply(theta2_SE_SSL,2,mean,na.rm=T)),
  #     main=paste('ASE sim:',sim),xlab='SE SUP',ylab='SE SSL')
  #abline(c(0,1))
  CI.95_SUP <- theta.SUP+cbind(LB=-1.96*SUP.SEs,UB=1.96*SUP.SEs)
  CI.95_SSL <- theta.SSL+cbind(LB=-1.96*SSL.SEs,UB=1.96*SSL.SEs)
  within_95CI_SUP[sim,] <- c(MC.theta$theta1,MC.theta$theta2)>=CI.95_SUP[,'LB'] & c(MC.theta$theta1,MC.theta$theta2)<=CI.95_SUP[,'UB']
  within_95CI_SSL[sim,] <- c(MC.theta$theta1,MC.theta$theta2)>=CI.95_SSL[,'LB'] & c(MC.theta$theta1,MC.theta$theta2)<=CI.95_SSL[,'UB']
  
  ####################################
  ##### Calculate Value function: ####
  ####################################
  
  # Propensity scores & omegas
  L_ls <- PS.treat(L_ls,PS1,PS2)
  L_ls <- omegas(df_ls=L_ls,theta1=Q1.est_SUP,theta2=Q2.est_SUP)
  # Compute supervised DR value estimate
  L_ls <- VsupDR(df_ls=L_ls,theta1=Q1.est_SUP,theta2=Q2.est_SUP)
  
  # Compute derivatives and influence function for Vsup:
  L_ls <- IF.derivs(L_ls,theta1=Q1.est_SUP,theta2=Q2.est_SUP)

  V.hat[sim,'SUP'] <- mean(L_ls$opt.SUP$VsupDR.hat)
  V.hat.SE[sim,'SUP'] <- sqrt(mean(L_ls$opt.SUP$V.psi^2)/n)
  
  #######################################
  ##### Calculate SS Value function: ####
  #######################################
  # Propensity scores & omegas
  U_ls <- PS.treat(U_ls,PS1,PS2,label=F)
  
  U_ls <- omegas(df_ls=U_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL,label=F)
  L_ls <- omegas(df_ls=L_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL,label=F)
  
  # Step 1) Imputation
  U_ls <- impute(L_ls,U_ls,model=impute_with,CV_grps=grps,Q_fun=F)
  # Compute Q opt for SE calculation
  L_ls <- VsupDR(df_ls=L_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL,label = F)
  
  # Step 2) Refitting step
  U_ls <- refitting(L_ls,U_ls,model=impute_with,CV_grps=grps,Q_fun=F)
  L_ls <- reffiting4L(L_ls,model=impute_with,CV_grps=grps,Q_fun=F)
  
  # Step 3) Projection step
  U_ls <- VsupDR(df_ls=U_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL,label = F)
  
  U_ls <- VsslDR(U_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL)
  
  L_ls <- IF.derivs(L_ls,theta1=Q1.est_SSL,theta2=Q2.est_SSL,label=F)
  
  V.hat[sim,'SSL'] <- mean(U_ls$opt.SSL$VsslDR.hat)
  V.hat.SE[sim,'SSL'] <- sqrt(mean(L_ls$opt.SSL$V.psi^2)/n)
  
  within_95CI_V[sim,'SSL'] <- V.star>=V.hat[sim,'SSL']-1.96*V.hat.SE[sim,'SSL'] & V.star<=V.hat[sim,'SSL']+1.96*V.hat.SE[sim,'SSL']
  within_95CI_V[sim,'SUP'] <- V.star>=V.hat[sim,'SUP']-1.96*V.hat.SE[sim,'SUP'] & V.star<=V.hat[sim,'SUP']+1.96*V.hat.SE[sim,'SUP']
  ###############################
  ######## print results ########
  ###############################
  
  
  if(sim %% print_every == 0){
    print(paste(sim,'N=',N,'n=',n,impute_with))
    cat('\nbias theta1:\n')
    print(cbind(SUP=apply(theta1_SUP[c(1:sim),],2,mean,na.rm=T)-MC.theta$theta1,
                SSL=apply(theta1_SSL[c(1:sim),],2,mean,na.rm=T)-MC.theta$theta1))
    
    # SE theta1
    cat('\nSE theta1:\n')

    print(cbind(SUP_ESE=apply(theta1_SUP[c(1:sim),],2,sd,na.rm=T),SUP_ASE=apply(theta1_SE_SUP,2,mean,na.rm=T),
                SSL_ESE=apply(theta1_SSL[c(1:sim),],2,sd,na.rm=T),SSL_ASE=apply(theta1_SE_SSL,2,mean,na.rm=T)))
    
    # bias theta2
    cat('\nbias theta2:\n')
    print(cbind(SUP=apply(theta2_SUP[c(1:sim),],2,mean,na.rm=T)-MC.theta$theta2,
                SSL=apply(theta2_SSL[c(1:sim),],2,mean,na.rm=T)-MC.theta$theta2))
    
    # SE theta2
    cat('\nESE theta2:\n')
    print(cbind(SUP_ESE=apply(theta2_SUP[c(1:sim),],2,sd,na.rm=T),SUP_ASE=apply(theta2_SE_SUP,2,mean,na.rm=T),
                SSL_ESE=apply(theta2_SSL[c(1:sim),],2,sd,na.rm=T),SSL_ASE=apply(theta2_SE_SSL,2,mean,na.rm=T)))
    
    cat('\nCovP SUP A2:\n')
    print(apply(within_95CI_SUP,2,mean,na.rm=T))
    cat('\nCovP SSL A2:\n')
    print(apply(within_95CI_SSL,2,mean,na.rm=T))
    
    cat('\nValue funcion bias:\n')
    print(apply(V.hat[c(1:sim),],2,mean,na.rm=T)-V.star)
    cat('\nValue funcion ESE:\n')
    print(apply(V.hat[c(1:sim),],2,sd,na.rm=T))
    cat('\nValue funcion ASE:\n')
    print(apply(V.hat.SE[c(1:sim),],2,mean,na.rm=T))
    cat('\nValue funcion CovP:\n')
    print(apply(within_95CI_V[c(1:sim),],2,mean,na.rm=T))
    cat('\n')
  }
  
  ###############################
  ######## save results ########
  ###############################
  
  if(sim %in% c(500,1000)){
    results <- list(MC.theta=MC.theta,MSE_sim=MSE_sim,
                    theta1_SUP=theta1_SUP,theta1_SSL=theta1_SSL,
                    theta1_SE_SUP=theta1_SE_SUP,theta1_SE_SSL=theta1_SE_SSL,
                    theta2_SUP=theta2_SUP,theta2_SSL=theta2_SSL,
                    theta2_SE_SUP=theta2_SE_SUP,theta2_SE_SSL=theta2_SE_SSL,
                    within_95CI_SUP=within_95CI_SUP,within_95CI_SSL=within_95CI_SSL,
                    V.hat=V.hat,V.hat.ASE=V.hat.SE,within_95CI_V=within_95CI_V)
    save(results,file=paste('../Results/SSRL_sims',sim.setting,'setting',
                        'sims_No',sim,'n',n,'N',N,
                        'impute_with',impute_with,'beta26',beta26,'phi24',phi24,format(Sys.time(), "%Y-%m-%d"),'.Rdata',sep='_'))
  }
  
}
