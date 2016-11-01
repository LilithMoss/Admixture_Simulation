library(BMA)
library(Rmpfr)
source("AdmixtureBMA.Simulations_Annotated.R")

#set.seed(2016)
d <- simulateData()

d$Yc <- d$Y-mean(d$Y) #Mean center the case-status variable
Y <- d$Y
Q <- d$Q
L <- d$Q
Yc <- d$Yc
pmw=c(1,1)
#Closed-form Simulation
pmw <- pmw/sum(pmw)
#PrMGivenD <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness))) #This is what needs to be replaced
n<-length(Y)
#Specify Phi
phi = 1#2.85
#Specify models
models <- rbind(c(0,1),c(1,1))

#Run regression models and extract values
numModels <- 2
reg <- as.list(rep(0, numModels))
reg[[1]] <- lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d) # Case-only 
reg[[2]] <- lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d) # case-control
betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))

#Specify Prior Covariance matrices - Used the continuous definition as Y has only 2 categories - ask Dr. Conti later
# sigma1 <- (summary(reg[[1]])$sigma)^2 #Estimated Variance of the standard Error
# sigma2 <- (summary(reg[[2]])$sigma)^2 #Estimated Variance of the standard Error
# Cov_0 <- list((phi^2)*(betas.se[1])^2*matrix(var(Q)^-1,ncol=sum(models[1,]),nrow=sum(models[1,])), #Uses sample variance of Q as prior variance (For continuous variables - not sure if this applies)
#               (betas.se[2])^2*diag(c(var(Q),(phi^2)*var(Yc)^-1)) )

#Covariance matrix based on data (residual standard error)
# prior.var <- c((0.05)^2,(0.05)^2)
# cov_1 <- sigma1*matrix((phi^2)*(1/prior.var[1]))
# cov_2.1 <- var(Q-L/2) #Using the data
# cov_2.2 <- ((phi^2)*(1/prior.var[2]))
# cov_2 <- sigma2*matrix(c(cov_2.1,0,0,cov_2.2),byrow=T,ncol=2)
# Cov_0 <- list(cov_1,cov_2)

#Prior variance of betas based on seeing maximum a 0.05 difference. 
#Informed Prior
sigma1 <- (0.05/1.96)^2
sigma2 <- (0.05/1.96)^2
#Covariance matrix without using any data
cov1 <- matrix(sigma1)
cov2 <- matrix(c(sigma1,0,0,sigma2),ncol=2,byrow=T)
Cov_0 <- list(cov1,cov2)

#Invert the covariance matrix to get precision matrix
lambda_0 <- list(solve(Cov_0[[1]]),solve(Cov_0[[2]]) )
#Specify Design Matrix (list of design matrices)
# X <- list( matrix(Yc),matrix(c(rep(1,length(Yc)),Yc),ncol=2) )
X <- list( matrix(Y),matrix(c(rep(1,length(Yc)),Y),ncol=2) )
#Specify Posterior Precision Matrix
lambda_n <- list(t(X[[1]])%*%X[[1]]+lambda_0[[1]],t(X[[2]])%*%X[[2]]+lambda_0[[2]])
#Specify prior mean vector for Betas
# mu_0 <- list(matrix(0),matrix(c(reg[[2]]$coefficients["(Intercept)"],0) ))
mu_0 <- list(matrix(0),matrix(c(0,0)))

#Specify posterior mean vector for betas
mu_n <- list( solve(t(X[[1]])%*%X[[1]]+lambda_0[[1]]) * ( (t(X[[1]])%*%X[[1]]*betas[1]) + (lambda_0[[1]]*mu_0[[1]]) ),
              solve(t(X[[2]])%*%X[[2]]+lambda_0[[2]]) %*% ( (t(X[[2]])%*%X[[2]] %*% matrix(reg[[2]]$coefficients)) + (lambda_0[[2]] %*% mu_0[[2]]) ) )
#Specify Prior hyperparameters for sigma^2
nu <- 2.58
lambda <- 0.28
a0 <- nu/2
b0 <- nu*lambda/2
#Specify Posterior hyperparameters for sigma^2
an <- a0+(n/2)
bn <- list( b0+(1/2)*(t(Q)%*%Q + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
            b0+(1/2)*(t(Q)%*%Q + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )

#Calculate large values using multiple precision package (Rmpfr)
lterm1 <- exp(as(((-n/2)*log(2*pi)),"mpfr")) #1/(2pi)^(n/2)
lterm2 <- list( exp(as(an*log(bn[[1]]),"mpfr")), exp(as(an*log(bn[[2]]),"mpfr")) )
lterm3 <- gamma(as(an,"mpfr"))

#Closed-Form WAY
#Calculate Marginal Likelihood
PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/gamma(a0)) 
PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/gamma(a0)) 
PrDGivenM <- c(PrDGivenM1,PrDGivenM2)
#Calculate Posterior Model Probabilities
PrMGivenD.new <- PrDGivenM*pmw/sum( PrDGivenM*pmw )
betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
post.beta <- sum(betas*PrMGivenD.new)
post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.new) - (post.beta^2))
z.score <- post.beta/post.se
p.value <- 2*pnorm(-abs(z.score))
r.cf <- c(post.beta, post.se, z.score, p.value)
r.cf

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
r.aic <- c(post.beta.aic, post.se.aic, z.score.aic, p.value.aic)
r.aic

#Glib WAY (Case-control model only)
models <- as.data.frame(1)
X <- as.data.frame(model.matrix(as.formula(Q ~ 1 + offset((as.numeric(L)/2)) + Y))[,2])
r.glib <- glib(X,y=QL,error="gaussian",link="identity",phi=1,
               models=models[1],pmw=pmw[1],output.postvar=T,
               priormean=c(0,0),priorvar=Cov_0[[2]])




