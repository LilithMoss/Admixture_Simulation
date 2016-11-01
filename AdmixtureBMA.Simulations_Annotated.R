########################################
# Simulation Code -Admixture Analysis  #
#
#
########################################

####### functions
simulateData <- function() {  #Assumes that scenario parameters have already been defined
	N <- N.controls+N.cases
	Y <- c(rep(1, N.cases), rep(0, N.controls)) #Case status 1 = Case
	Z <- ifelse(Y==1,0,1) #Control status; 1 = Control
	mu.Q <- Q.mean + Y*beta.Q.cases + Z*beta.Q.controls
	Q <- rnorm(N, mean=mu.Q, sd=Q.sd) #Global Ancestry
	mu.L <- 2*mu.Q + Y*beta.L.cases + Z*beta.L.controls
	L <- rnorm(N, mean=mu.L, sd=L.sd) #Local Ancestry
	d <- data.frame(Y=Y, Z=Z, Q=Q, L=L)
	d
}

simulateData.new <- function() {  #Assumes that scenario parameters have already been defined
  N <- N.controls+N.cases
  Y <- c(rep(1, N.cases), rep(0, N.controls)) #Case status 1 = Case
  # Z <- ifelse(Y==1,0,1) #Control status; 1 = Control
  mu.Q <- Q.mean + Y*beta.Q.cases + Z*beta.Q.controls
  Q <- rnorm(N, mean=mu.Q, sd=Q.sd) #Global Ancestry
  mu.L <- 2*mu.Q + Y*beta.L.cases + Z*beta.L.controls
  L <- rnorm(N, mean=mu.L, sd=L.sd) #Local Ancestry
  d <- data.frame(Y=Y, Z=Z, Q=Q, L=L)
  d
}

#Simulate the difference
# alpha.controls <- 0.00 #mean in controls
# beta.cases <- 0.05 #mean difference between cases and controls
# mu.Q.cases <- 0.8
# mu.Q.controls <- 0.8
# mu.L.cases <- 2*(mu.Q.cases-beta.cases)
# sd.Q.cases <- 0.14
# sd.Q.controls <- 0.14
# 
# 
# Q.cases <- rnorm(N.cases,mean=mu.Q.cases,sd=sd.Q.cases) #Global ancestry in cases
# Q.controls <- rnorm(N.controls,mean=mu.Q.controls,sd=sd.Q.controls) #Global ancestry in controls



# Simulation Parameters
numSims <- 100	
# Original Simulation
# N.cases <- 3000
# N.controls <- 3000
#Simulation for 80% power
N.cases <- 3200
N.controls <- 3200
# N.cases <- 100
# N.controls <- 100
Q.mean <- .8
Q.sd <- 0.136
L.sd <- 0.6
beta.Q.cases <- 0.0
beta.Q.controls <- 0.0
beta.L.cases <- 0.05
beta.L.controls <- 0.0

runSim <- function(sim) {
	if((sim %% 100) == 0) { print(sim) }
	d <- simulateData()
	#r.case tests the mean of the difference in global and local
	#ancestry in cases only to see if this difference is significanly
	#different from zero.
	#Takes out intercept (that's all there is)
	r.case <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d))$coef[1,]
	
	#r.casecontrol tests the effect of Local ancestry on case status
	#adjusted for global ancestry
	#Takes out L (out of intercept, Q, L)
	r.casecontrol <- summary(glm(Y ~ Q + L, family="binomial", data=d))$coef[3,]

	#r.control tests the mean of the difference in global and local
	#ancestry in controls only to see if this difference is significanly
	#different from zero. Exactly like the case only except with controls
	#Takes out intercept (that's all there is) 
	r.control <- summary(lm(Q~1+offset(L/2), subset=Z==1, data=d))$coef[1,]
	
	#r.casecontrol.lin tests the effect of case status on the difference in 
	#local and global ancestry. 
	#Takes out Y (case-status)
	#r.casecontrol.lin <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d))$coef[2,] # case-control
	
	#r.BMA
	r.BMA <- runBMA.Admixture(d)
	r <- list(r.case=r.case, r.casecontrol=r.casecontrol, r.control=r.control, r.BMA=r.BMA)
	r
}

#############################################
# BMA for admixture is combining 2 models
# Model 1: 
runBMA.Admixture <- function(d, pmw=c(1,1)) {
	numModels <- 2
	reg <- as.list(rep(0, numModels))
	reg[[1]] <- lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d) # Case-only 
	reg[[2]] <- lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d) # case-control
	
	#Exploring what the models mean
	# Case only: the only difference in the 2 models is that the intercept in r.case 
	# is the Y coefficient in model 1. The standard error is slightly lower for m1.
	# m1 essentially tests whether the mean of the Q-L/2 difference in cases is 0
	m1 <- summary(lm(Q ~ -1 + offset((as.numeric(L)/2)) + Y, data=d))  #BMA function case-only
	r.case <- summary(lm(Q~1+offset(L/2), subset=Y==1, data=d)) # Original case-only
	
	#Case-control
	# m2 uses global ancestry as the output while r.casecontrol uses Y
	# They produce tremendously different results.
	m2 <- summary(lm(Q ~ 1 + offset((as.numeric(L)/2)) + Y, data=d)) # case-control
	r.casecontrol <- summary(glm(Y ~ Q + L, family="binomial", data=d))
	# m1 is testing the intercept and effect of case-status on global ancestry to see 
	# if they are zero. What is the intercept in m1? It is the mean difference of 
	# Q-L/2 in controls. And Y is the effect of being case on global ancestry.
	# diff <- d$Q-d$L/2
	# controls <- d[d$Z==1,]
	# diff <- controls$Q-controls$L/2
	cases <- d[d$Y==1,]
	diff <- cases$Q-cases$L/2
	
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
}

	



#  beta.Q.cases.list <- seq(from=0, to=0, by=.04)
#  beta.Q.controls.list <- seq(from=0, to=0, by=.04)
#  beta.L.cases.list <- seq(from=-0.04, to=0.04, by=.01)
#  beta.L.controls.list <- seq(from=-0.04, to=0.04, by=.01)
# 
# beta.Q.cases.list <- c(0.0)
# beta.Q.controls.list <- c(0.0)
# beta.L.cases.list <- c(0.0, 0.04)
# beta.L.controls.list <- c(0.0, 0.02, 0.04)
#  print(paste("Total Number of Scenarios:", length(beta.Q.cases.list)*length(beta.Q.controls.list)*length(beta.L.cases.list)*length(beta.L.controls.list)))
#  OverallResults <- {}
#  scenarioNum <- 1
#  for(beta.Q.cases in beta.Q.cases.list) {
#  	for(beta.Q.controls in beta.Q.controls.list) {
#  		for(beta.L.cases in beta.L.cases.list) {
#  			for(beta.L.controls in beta.L.controls.list) {
#  				print(paste("Scenario:", scenarioNum,
#  							"beta.Q.cases:", beta.Q.cases, 
#  							"beta.Q.controls:", beta.Q.controls,
#  							"beta.L.cases:", beta.L.cases,
#  							"beta.L.controls:", beta.L.controls))
#  				system.time(results <- lapply(1:numSims, runSim))
#  				r.case <- matrix(unlist(lapply(results, FUN=function(v) { v$r.case })), 4, numSims)
#  				r.casecontrol <- matrix(unlist(lapply(results, FUN=function(v) { v$r.casecontrol })), 4, numSims)
#  				r.BMA <- matrix(unlist(lapply(results, FUN=function(v) { v$r.BMA })), 4, numSims)
#  				r.control <- matrix(unlist(lapply(results, FUN=function(v) { v$r.control })), 4, numSims)
#  				OverallResults <- rbind(OverallResults, 
#  							c(beta.Q.cases, beta.Q.controls, beta.L.cases, beta.L.controls,
#  							mean(r.case[1,]), mean(r.case[2,]), sum(r.case[4,]<0.05)/numSims,
#  							mean(r.casecontrol[1,]), mean(r.casecontrol[2,]), sum(r.casecontrol[4,]<0.05)/numSims,
#  							mean(r.control[1,]), mean(r.control[2,]), sum(r.control[4,]<0.05)/numSims,
#  							mean(r.BMA[1,]), mean(r.BMA[2,]), sum(r.BMA[4,]<0.05)/numSims))
#  				scenarioNum <- scenarioNum + 1
#  			}	
#  		}
#  	}	
#  }
#  OverallResults <- as.data.frame(OverallResults)
#  names(OverallResults) <- c("beta.Q.cases", "beta.Q.controls", "beta.L.cases", "beta.L.controls",
#  								"r.case.beta", "r.case.se", "r.case.power",
#  								"r.casecontrol.beta", "r.casecontrol.se", "r.casecontrol.power",
#  								"r.control.beta", "r.control.se", "r.control.power",
#  								"r.BMA.beta", "r.BMA.se", "r.BMA.power")
#  write.table(OverallResults, file="BMA.Admixture.Results.txt", quote=F, row.names=F, sep="\t")
# 
#  #What is this actually doing
#  r.case <- matrix(unlist(lapply(results, FUN=function(v) { v$r.case })), 4, numSims)
#  results[[1]]$r.case
#  
#  
 
 
