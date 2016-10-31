######################################
library(BMA)
source("AdmixtureBMA.Simulations_Annotated.R")
set.seed(2016)
d <- simulateData()

d$Yc <- d$Y-mean(d$Y) #Mean center the case-status variable
Y <- d$Y
Q <- d$Q
L <- d$Q
Yc <- d$Yc
E <- rnorm(n,)

pmw=c(1,1)
#Closed-form Simulation
pmw <- pmw/sum(pmw)
#PrMGivenD <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness))) #This is what needs to be replaced
n<-length(Y)
#Specify Phi
phi = 1
#Specify models
models <- rbind(c(0,1),c(1,1))
