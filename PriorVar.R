library(forward)
data(vaso)
x<- vaso[,1:2]
y<- vaso[,3]
n<- rep(1,times=length(y))
finney.models<- rbind(
  c(1, 0),
  c(0, 1),
  c(1, 1))
finney.glib <- glib (x,y,n, error="binomial", link="logit",
                     models=finney.models, glimvar=TRUE,
                     output.priorvar=TRUE, output.postvar=TRUE)
summary(finney.glib)

#Hyper-parameters
models <- rbind(
  c(1, 0),
  c(0, 1),
  c(1, 1))
phi=c(2,2)
nmodel <- nrow(models)
nphi <- length(phi)
chi2 <- rep(0, nmodel)
npar <- rep(0, nmodel)
prior.var <- as.list(rep(0, nmodel))

model <- (1:ncol(x))[models[i, ] == 1]
npar[i] <- length(model) + 1
prior.var[[i]] <- array(rep(0, npar[i] * npar[i] * 
                              nphi), dim = c(npar[i], npar[i], nphi))

for (j in (1:nphi)) prior.var[[i]][, , j] <- glim.pvar(model, 
                                                       phi[j], psi, Aot, prior.spec)$pvar
#Simple Simulation
set.seed(2000)
L <- rnorm(100,1,0.5)
Q <- runif(100,0,1)
Y <- rbinom(100,1,0.2)
Yc <- Y-mean(Y)
dat <- cbind(L,Q,Y)

#Closed-form Simulation
pmw <- pmw/sum(pmw)
#PrMGivenD <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness))) #This is what needs to be replaced
n<-length(Y)
#Specify Phi
phi = 1
#Specify models
models <- rbind(c(0,1),c(1,1))
#Specify Prior Covariance matrices - Used the continuous definition as Y has only 2 categories - ask Dr. Conti later
Cov_0 <- list((phi^2)*betas.se[1]*matrix(var(Q),ncol=sum(models[1,]),nrow=sum(models[1,])), #Uses sample variance of Q as prior variance (For continuous variables - not sure if this applies)
              (phi^2)*betas.se[2]*diag(c(var(Q),var(Yc))) )
#Invert the covariance matrix to get precision matrix
lambda_0 <- list(solve(Cov_0[[1]]),solve(Cov_0[[2]]) )
#Specify Design Matrix (list of design matrices)
X <- list( matrix(Yc),matrix(c(rep(1,length(Yc)),Yc),ncol=2) )
#Specify Posterior Precision Matrix
lambda_n <- list(t(X[[1]])%*%X[[1]]+lambda_0[[1]],t(X[[2]])%*%X[[2]]+lambda_0[[2]])
#Specify prior mean vector for Betas
mu_0 <- list(matrix(0),matrix(c(reg[[2]]$coefficients["(Intercept)"],0) ))
#Specify posterior mean vector for betas
mu_n <- list( solve(t(X[[1]])%*%X[[1]]+lambda_0[[1]]) * ( (t(X[[1]])%*%X[[1]]*betas[1]) + (lambda_0[[1]]*mu_0[[1]]) ),
              solve(t(X[[2]])%*%X[[2]]+lambda_0[[2]]) %*% ( (t(X[[2]])%*%X[[2]] %*% matrix(reg[[2]]$coefficients)) 
                                                           + (lambda_0[[2]] %*% mu_0[[2]]) ) )
              
#Specify Prior hyperparameters for sigma^2
nu <- 2.58
lambda <- 0.28
a0 <- nu/2
b0 <- nu*lambda/2
#Specify Posterior hyperparameters for sigma^2
an <- a0+(n/2)
bn <- list( b0+(1/2)*(t(Q)%*%Q + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
            b0+(1/2)*(t(Q)%*%Q + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
#Calculate Marginal Likelihood
PrDGivenM1 <- 1/(2*pi)^(n/2)*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*(b0^a0)/(bn[[1]]^an)*(gamma(an)/gamma(a0))  
PrDGivenM2 <- 1/(2*pi)^(n/2)*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*(b0^a0)/(bn[[2]]^an)*(gamma(an)/gamma(a0))  
PrDGivenM <- c(PrDGivenM1,PrDGivenM2)
#Calculate Posterior Model Probabilities
PrMGivenD <- PrDGivenM*pmw/sum( PrDGivenM*pmw )

#Keep the rest of the calculation the same
betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
betas.se <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))
post.beta <- sum(betas*PrMGivenD)
post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD) - (post.beta^2))
z.score <- post.beta/post.se
p.value <- 2*pnorm(-abs(z.score))
r <- c(post.beta, post.se, z.score, p.value)
r







