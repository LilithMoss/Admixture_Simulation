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
L <- rnorm(100,)