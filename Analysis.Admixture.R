#########################################################################
# Analysis.Admixture.R
# House functions for Admixture Simulation and Analysis for single locus
# Describe functions for case control, case only, and BMA
#########################################################################
library(BMA)
library(Rmpfr)

# Data Simulation Function
simulateData <- function() {  #Assumes that scenario parameters have already been defined
  N <- N.controls+N.cases
  Y <- c(rep(1, N.cases), rep(0, N.controls)) #Case status 1 = Case
  Z <- ifelse(Y==1,0,1) #Control status; 1 = Control
  mu.Q <- Q.mean + Y*beta.Q.cases + Z*beta.Q.controls
  Q <- rnorm(N, mean=mu.Q, sd=Q.sd) #Global Ancestry
  #mu.L <- 2*mu.Q + Y*beta.L.cases + Z*beta.L.controls
  mu.L <- 2*mu.Q - Y*beta.L.cases - Z*beta.L.controls
  L <- rnorm(N, mean=mu.L, sd=L.sd) #Local Ancestry
  d <- data.frame(Y=Y, Z=Z, Q=Q, L=L)
  d
}

# #TEST############################################################
# #Case-Only
# # r.case <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d))$coef[1,]
# r.case.only <- summary(lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d) )$coefficients
#
# m1 <- lm((Q-0.5*L)~-1+Y)
# summary(m1)
#
# r.case <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d))
#
#
# #Case-Control
# # r.casecontrol.lin <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d))$coef[2,] # case-control
# r.case.control <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d) )$coefficients
#
#
#
#
# simulateData.simple <- function(){
#   N <- N.controls+N.cases
#   Y <- c(rep(1, N.cases), rep(0, N.controls)) #Case status 1 = Case
#   Z <- ifelse(Y==1,0,1) #Control status; 1 = Control
#   diff <-
#
#
# }
# #TEST############################################################

run.BMA <- function(d){
  Y <- d$Y
  Q <- d$Q
  L <- d$L
  pmw=c(1,1)
  pmw <- pmw/sum(pmw)
  n<-length(Y)
  phi = 1#2.85
  models <- rbind(c(0,1),c(1,1))
  
  #Run regression models and extract values
  numModels <- 2
  reg <- as.list(rep(0, numModels))
  reg[[1]] <- lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d) # Case-only 
  reg[[2]] <- lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d) # case-control
  betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
  betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
  
  #Informed Prior
  # sigma1 <- (0.05/1.96)^2
  # sigma2 <- (0.05/1.96)^2
  
  #Covariance matrix without using any data
  cov1 <- matrix(sigma.1)
  cov2 <- matrix(c(sigma.1,0,0,sigma.2),ncol=2,byrow=T)
  Cov_0 <- list(cov1,cov2)
  
  #Invert the covariance matrix to get precision matrix
  lambda_0 <- list(solve(Cov_0[[1]]),solve(Cov_0[[2]]) )
  #Specify Design Matrix (list of design matrices)
  X <- list( matrix(Y),matrix(c(rep(1,length(Y)),Y),ncol=2) )
  #Specify Posterior Precision Matrix
  lambda_n <- list(t(X[[1]])%*%X[[1]]+lambda_0[[1]],t(X[[2]])%*%X[[2]]+lambda_0[[2]])
  #Specify prior mean vector for Betas
  mu_0 <- list(matrix(0),matrix(c(0,0)))
  #Specify posterior mean vector for betas
  mu_n <- list( solve(t(X[[1]])%*%X[[1]]+lambda_0[[1]]) * ( (t(X[[1]])%*%X[[1]]*betas[1]) + (lambda_0[[1]]*mu_0[[1]]) ),
                solve(t(X[[2]])%*%X[[2]]+lambda_0[[2]]) %*% ( (t(X[[2]])%*%X[[2]] %*% matrix(reg[[2]]$coefficients)) + (lambda_0[[2]] %*% mu_0[[2]]) ) )
  #Specify Prior hyperparameters for sigma^2
  #nu <- 2.58
  #lambda <- 0.28
  
  #Try setting hyper-parameters
  # nu <- 4.58
  # lambda <- 2*1.29*(sqrt(sigma1))/nu
  # lambda <- 2*1.29*(sqrt(sigma1))/nu
  # a0 <- nu/2
  # b0 <- nu*lambda/2
  
  #Just set a0 and b0 directly (without hyperparameters)
  #a0 = 2.5
  b0 = 1
  a0 = 1/(sigma.1) + 1
  
  #Specify Posterior hyperparameters for sigma^2
  an <- a0+(n/2)
  bn <- list( b0+(1/2)*(t(Q)%*%Q + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
              b0+(1/2)*(t(Q)%*%Q + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
  
  #Calculate large values using multiple precision package (Rmpfr)
  lterm1 <- exp(as(((-n/2)*log(2*pi)),"mpfr")) #1/(2pi)^(n/2)
  lterm2 <- list( exp(as(an*log(bn[[1]]),"mpfr")), exp(as(an*log(bn[[2]]),"mpfr")) )
  lterm3 <- gamma(as(an,"mpfr"))
  lterm4 <- gamma(as(a0,"mpfr"))
  
  #Closed-Form WAY
  #Calculate Marginal Likelihood
  # PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/gamma(a0)) 
  # PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/gamma(a0)) 
  PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/lterm4) 
  PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/lterm4) 
  PrDGivenM <- c(PrDGivenM1,PrDGivenM2)
  #Calculate Posterior Model Probabilities
  PrMGivenD.new <- PrDGivenM*pmw/sum( PrDGivenM*pmw )
  # betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
  # betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
  post.beta <- sum(betas*PrMGivenD.new)
  post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.new) - (post.beta^2))
  z.score <- post.beta/post.se
  p.value <- 2*pnorm(-abs(z.score))
  #BMA.cf <- as.numeric( c(post.beta, post.se, z.score, p.value))
  BMA.cf <- as.numeric( c(post.beta, post.se, z.score, p.value,
                          as.numeric(PrMGivenD.new[1]),as.numeric(PrMGivenD.new[2])))
  names(BMA.cf) <- c("post.beta", "post.se", "z.score", "p.value",
                     "PrDGivenM1","PrDGivenM1")
  BMA.cf
  
  #AIC WAY
  ll <- unlist(lapply(reg, AIC))
  fitness <- ll-log(pmw)
  PrMGivenD.aic <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness)))
  betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
  betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
  post.beta.aic <- sum(betas*PrMGivenD.aic)
  post.se.aic <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.aic) - (post.beta.aic^2))
  z.score.aic <- post.beta.aic/post.se.aic
  p.value.aic <- 2*pnorm(-abs(z.score.aic))
  BMA.aic <- c(post.beta.aic, post.se.aic, z.score.aic, p.value.aic,
               PrMGivenD.aic[1],PrMGivenD.aic[2])
  names(BMA.aic) <- c("post.beta.aic", "post.se.aic", "z.score.aic", "p.value.aic",
                      "PrDGivenM1","PrDGivenM1")
  BMA.aic
  result <- list(BMA.cf,BMA.aic)
  return(result)
}

runSim <- function(sim) {
  # if((sim %% 100) == 0) { print(sim) }
  d <- simulateData()
  #Case-Only
  # r.case <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d))$coef[1,]
  r.case.only <- summary(lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d) )$coefficients 
  #Case-Control
  # r.casecontrol.lin <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d))$coef[2,] # case-control
  r.case.control <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d) )$coefficients 
  #Control-Only
  r.control.only <- summary(lm(Q~1+offset(L/2), subset=Z==1, data=d))$coef[1,]
  #Logistic case-control
  r.logistic.cc <- summary(glm(Y ~ Q + L, family="binomial", data=d))$coef[3,]
  #BMA - 2 Types cf(closed-form formula) aic(aic estimation)
  BMA <- run.BMA(d)
  r.BMA.cf <- BMA[[1]]
  r.BMA.aic <- BMA[[2]]
  #r <- list(r.case=r.case, r.casecontrol=r.casecontrol, r.control=r.control, r.BMA=r.BMA)
  r <- list(r.case.only,r.case.control,r.control.only,r.logistic.cc,
            r.BMA.cf,r.BMA.aic)
  names(r) <- c("case.only","case.control","control.only",
                "logistic.cc","BMA.cf","BMA.aic")
  r
}
