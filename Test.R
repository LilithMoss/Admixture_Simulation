data(trees)
attach(trees)
library(pscl)

model <- lm(Volume~Girth+Height)
names(model)
names(summary(model))
summary(model)
summary(model)$sigma
sqrt(var(resid(model)))
(sd(resid(model)))

beta.hat <- model$coef
n <- length(Volume)
p <- length(beta.hat)
s2 <- (n-p)*summary(model)$sigma^2
V.beta <- summary(model)$cov.unscaled

temp1 <- rgamma(1,shape=(n-p)/2,rate=s2/2)
temp2 <- rigamma(1,(n-p)/2,s2/2)
