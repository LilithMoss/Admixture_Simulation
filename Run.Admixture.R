##############################################
# Run.Admixture.R
# Run Admixture Simulation for single locus
# Run functions for case control, case only, and BMA
# Calculate and compare power for all methds
##############################################
#Source Functions
source("Analysis.Admixture.R")

# Set Parameters
numSims <- 100	
N.cases <- 3200
N.controls <- 3200
Q.mean <- .8
Q.sd <- 0.136
L.sd <- 0.6
beta.Q.cases <- 0.0
beta.Q.controls <- 0.0
beta.L.cases <- 0.05
beta.L.controls <- 0.0

#Run Only one Simulation and Analysis
sim=1
#Returns the results of 
# 1) r.case.only
# 2) r.case.control
# 3) r.control.only
# 4) r.logistic.cc
# 5) r.BMA.cf
# 6) r.BMA.aic
runSim(sim)

#Run multiple simulations and report power
#Set number of Simulations
numSims = 2

#Set sets of parameters
LEFT OFF HERE!!!

results <- lapply(1:numSims, runSim)
#Extract pieces of results by method
case.only <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$case.only }))
case.control <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$case.control }))
control.only <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$control.only }))
logistic.cc <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$logistic.cc }))
BMA.cf <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$BMA.cf }))
BMA.aic <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$BMA.aic }))

#Bind all results together (By column)
OverallResults <- as.data.frame( cbind(beta.Q.cases, beta.Q.controls, beta.L.cases, beta.L.controls, #Simulation Inputs
  mean(case.only[,1]), mean(case.only[,2]), sum(case.only[,4]<0.05)/numSims,
  mean(case.control[,1]), mean(case.control[,2]), sum(case.control[,4]<0.05)/numSims,
  mean(control.only[,1]), mean(control.only[,2]), sum(control.only[,4]<0.05)/numSims,
  mean(logistic.cc[,1]), mean(logistic.cc[,2]), sum(logistic.cc[,4]<0.05)/numSims,
  mean(BMA.cf[,1]), mean(BMA.cf[,2]), sum(BMA.cf[,4]<0.05)/numSims,
  mean(BMA.aic[,1]), mean(BMA.aic[2]), sum(BMA.aic[,4]<0.05)/numSims) )

names(OverallResults) <- c("beta.Q.cases", "beta.Q.controls", "beta.L.cases", "beta.L.controls",
  "case.only.beta", "case.only.se", "case.only.power",
  "case.control.beta", "case.control.se", "case.control.power",
  "control.only.beta", "control.only.se", "control.only.power",
  "logistic.cc.beta", "logistic.cc.se", "logistic.cc.power",
  "BMA.cf.beta", "BMA.cf.se", "BMA.cf.power",
  "BMA.aic.beta", "BMA.aic.se", "BMA.aic.power")

