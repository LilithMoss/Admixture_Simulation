library(BMA)
library(Rmpfr)
source("AdmixtureBMA.Simulations_Annotated.R")

set.seed(2016)
d <- simulateData()

#Simple Simulation
# L <- rnorm(100,1,0.5)
# Q <- runif(100,0,1)
# Y <- rbinom(100,1,0.2)
# dat <- cbind(L,Q,Y)
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
cases <- d[d$Y==1,]
Qcases <- cases$Q
sigma1 <- (summary(reg[[1]])$sigma)^2
sigma2 <- (summary(reg[[2]])$sigma)^2
# Cov_0 <- list((phi^2)*(betas.se[1])^2*matrix(var(Q)^-1,ncol=sum(models[1,]),nrow=sum(models[1,])), #Uses sample variance of Q as prior variance (For continuous variables - not sure if this applies)
#               (betas.se[2])^2*diag(c(var(Q),(phi^2)*var(Yc)^-1)) )
prior.var <- c((0.05)^2,(0.05)^2)
prior.var.alpha <- (0.01)^2
cov_1 <- sigma1*matrix((phi^2)*(1/prior.var[1]))
cov_2.1 <- var(Q-L/2)
cov_2.2 <- ((phi^2)*(1/prior.var[2]))
cov_2 <- sigma2*matrix(c(cov_2.1,0,0,cov_2.2),byrow=T,ncol=2)
Cov_0 <- list(cov_1,cov_2)
  
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

#This is the code that works
# mpfr(exp(750),128)
# gamma(as(3000,"mpfr"))

#Calculate large values using multiple precision package (Rmpfr)
lterm1 <- exp(as(((-n/2)*log(2*pi)),"mpfr")) #1/(2pi)^(n/2)
lterm2 <- list( exp(as(an*log(bn[[1]]),"mpfr")), exp(as(an*log(bn[[2]]),"mpfr")) )
lterm3 <- gamma(as(an,"mpfr"))

#Calculate Marginal Likelihood
PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/gamma(a0)) 
PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/gamma(a0)) 
PrDGivenM <- c(PrDGivenM1,PrDGivenM2)

# #Calculate Marginal Likelihood - OLD
# PrDGivenM1 <- 1/(2*pi)^(n/2)*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*(b0^a0)/(bn[[1]]^an)*(gamma(an)/gamma(a0))  
# PrDGivenM2 <- 1/(2*pi)^(n/2)*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*(b0^a0)/(bn[[2]]^an)*(gamma(an)/gamma(a0))  
# PrDGivenM <- c(PrDGivenM1,PrDGivenM2)

#Calculate Posterior Model Probabilities
PrMGivenD.new <- PrDGivenM*pmw/sum( PrDGivenM*pmw )

#Keep the rest of the calculation the same
betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
post.beta <- sum(betas*PrMGivenD.new)
post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.new) - (post.beta^2))
z.score <- post.beta/post.se
p.value <- 2*pnorm(-abs(z.score))
r.new <- c(post.beta, post.se, z.score, p.value)
r.new

#Old way
ll <- unlist(lapply(reg, AIC))
pmw <- pmw/sum(pmw)
fitness <- ll-log(pmw)
PrMGivenD <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness)))
betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
post.beta <- sum(betas*PrMGivenD)
post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD) - (post.beta^2))
z.score <- post.beta/post.se
p.value <- 2*pnorm(-abs(z.score))
r <- c(post.beta, post.se, z.score, p.value)
r

#Glib to check -- Don't do this - it's too much work.
r.co <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d))
r.cc <- summary(lm(Q~1+offset(L/2)+Y, data=d))
d$ex <- ifelse(d$Y==0,NA,1)
x <- cbind(d$Y,d$ex) #Design Matrix
y <- (d$Q-d$L/2)   #Outcome
models <- rbind(c(0),c(1))
n=nrow(d)
r.glib <- glib(x[,1],y,n, error="gaussian", link="identity",
          models=models, glimvar=TRUE,
          output.priorvar=TRUE, output.postvar=TRUE)
summary(glib)

x <- cbind(rnorm(100,2,0.5),rnorm(100,3,0.2))
y <- rnorm(100,10,0.5)
n=100
models <- rbind(c(1,0),c(1,1))
r.glib <- glib(x,y,n, error="gaussian", link="identity",
               models=models, glimvar=TRUE, nu=2.58, phi=2.85,
               output.priorvar=TRUE, output.postvar=TRUE)
summary(r.glib)
#calculate what the posterior mean and variances would be given the glib priors


library(forward)
data(vaso)
x<- vaso[,1:2]
y<- vaso[,3]
n<- rep(1,times=length(y))
models<- rbind(c(0,1),c(1,1))
glib <- glib (x,y,n, error="gaussian", link="identity",
                     models=models, glimvar=TRUE,
                     output.priorvar=TRUE, output.postvar=TRUE)
summary(glib)

m <- summary(lm(y~x[,1]+x[,2]))




