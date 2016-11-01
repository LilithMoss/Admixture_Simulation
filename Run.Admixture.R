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

#Run Simulation and Analysis
sim=1
runSim(sim)
