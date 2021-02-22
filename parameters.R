# parameters:
#######################
### schulte et. al. ###
#######################
params = 'schulte'

if(params=='schulte'){
  # A1:
  phi1 <- c(.3,-.5)
  names(phi1) <- c('10','11')
  # Y1:
  beta1 <- c(1,1,1,-2)#beta1 <- c(1,1,1,1)#
  names(beta1) <- c(10:13)
  # O2:
  delta1 <- c(0, .5,-.75, .25,1)
  names(delta1) <- c(10:14)
  
  #A2:
  phi2 <- c(0, .5, .1,-1,phi24,-.1,1)
  names(phi2) <- c(20:26)
  
  #Y2
  beta2 <- c(3, 0, .1,-.5,-.5,.1,beta26)
  psi2 <- c(1, .25, .5)#psi2 <- c(1, .25, 1.5)#
  names(beta2) <- c(20:26)
  names(psi2) <- c(20:22)
}

  #################
  ### all ones ###
  #################
if(params=='ones'){
  # A1:
  phi1 <- c(1,1)
  names(phi1) <- c('10','11')
  # Y1:
  beta1 <- c(1,1,1,-2)
  names(beta1) <- c(10:13)
  # S2:
  delta1 <- c(0, 1,1, 1,1)
  names(delta1) <- c(10:14)
  
  #A2:
  phi25 <- 0
  phi2 <- c(0, 1, 1,1,1,phi25,1)
  names(phi2) <- c(20:26)
  
  #Y2
  beta2 <- c(1, 1, 1,1,1,1)
  psi2 <- c(1, 1, 1)
  names(beta2) <- c(20:25)
  names(psi2) <- c(20:22)
}
#true_params_stage2 <- c(beta2['20'],beta2['21'],beta2['22'],beta2['23'],beta2['24'],beta2['25'],beta2['26'],psi2['20'],psi2['21'],psi2['22'])
#true_params_stage1 <- c(beta1['10'],beta1['11'],beta1['12'],beta1['13'],beta1['14'])

# Correctly speciffied models:
# Y1 = beta10 + beta11*O1 + beta12*A1 + beta13*O1*A1 + beta14*O1^2
# Y2 = beta20 + beta21*O1 + beta22*A1 + beta23*O1*A1 + beta24*O2 + beta25*Y1 + 
#      psi20*A2* + psi21*A1*A2 + psi22*A2*O2        
varsQ1 <- c('O1','A1','O1xA1')
varsQ2_np <- c('O1','A1','O1xA1','O2','imputed_Y1','A2','A2xA1','A2xO2')
varsQ2_ols <- c('O1','A1','O1xA1','O2','Y1','A2','A2xA1','A2xO2')

#if(misspecified) varsQ1 <- varsQ1[-3]; varsQ2_np <- varsQ2_np[-3]

model.Q1 <- as.formula(paste('pseudo_outcome', 
                             paste(varsQ1, collapse = " + "), 
                             sep = " ~ "))




model.Q2_np <- as.formula(paste('imputed_Y2', 
                                paste(varsQ2_np, collapse = " + "),
                             sep = " ~ "))
model.Q2_ols <- as.formula(paste('Y2', 
                                 paste(varsQ2_ols, collapse = " + "),  
                                sep = " ~ "))


# Covariates used for imputation:
# For imputing Y1
available_covs <- c("O1","A1","O1xA1","O2","A2","A2xA1","A2xO2",'W1','W2')
#available_cov_U1 <- c("O1","A1","O1xA1")
# For imputing Y2


# Next steps: 

# 1) Do one simulation with correctly specified model with the three imputation methods and see what happens
# Make this script flexible enough (with if functions) to have correct and misspecified models without having to edit the whole script and simulation functions
# 2) Do a misspecified simulation where the parameters of the misspecified model you estimated them via simulation of huge dataset and then a regular OLS fit, this way the mse is compared to the correct parameters