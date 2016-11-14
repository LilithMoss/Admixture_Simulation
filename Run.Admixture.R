##############################################
# Run.Admixture.R
# Run Admixture Simulation for single locus
# Run functions for case control, case only, and BMA
# Calculate and compare power for all methds
##############################################
library(ggplot2)
library(reshape)
library(plyr)
library(gridExtra)
#Source Functions
source("Analysis.Admixture.R")
set.seed(50)
# Set Parameters
#numSims <- 100	
N.cases <- 3200
N.controls <- 3200
Q.mean <- .8
Q.sd <- 0.136
L.sd <- 0.6
# beta.Q.cases <- 0.0
# beta.Q.controls <- 0.0
# beta.L.cases <- 0.05
# beta.L.controls <- 0.0

#Run Only one Simulation and Analysis
#sim=1
#Returns the results of 
# 1) r.case.only
# 2) r.case.control
# 3) r.control.only
# 4) r.logistic.cc
# 5) r.BMA.cf
# 6) r.BMA.aic
#runSim(1)

#Run multiple simulations and report power
#Set number of Simulations
numSims = 100

#Set sets of parameters
beta.Q.cases.list <- c(0.0) #Q in cases
beta.Q.controls.list <- c(0.0) #Q in controls
#beta.L.cases.list <- c(0.0, 0.02, 0.04, 0.05, 0.06) #L effect in cases
beta.L.cases.list <- seq(0.03,0.03,by=0.01)
beta.L.controls.list <- seq(0.0, 0.07, by=0.01) #L effect in controls
#beta.L.controls.list <- c(0.00)

#All possible combinations
param <- expand.grid(beta.Q.cases.list,beta.Q.controls.list,
            beta.L.cases.list,beta.L.controls.list)
names(param) <- c("beta.Q.cases","beta.Q.controls",
                  "beta.L.cases","beta.L.controls")

OverallResults <- data.frame( beta.Q.cases=numeric(0), beta.Q.controls=numeric(0), 
                beta.L.cases=numeric(0), beta.L.controls=numeric(0), 
                case.only.beta=numeric(0), case.only.se=numeric(0), 
                case.only.power=numeric(0),case.control.beta=numeric(0),
                case.control.se=numeric(0), case.control.power=numeric(0),
                control.only.beta=numeric(0), control.only.se=numeric(0), 
                control.only.power=numeric(0),logistic.cc.beta=numeric(0), 
                logistic.cc.se=numeric(0), logistic.cc.power=numeric(0),
                BMA.cf.beta=numeric(0), BMA.cf.se=numeric(0), 
                BMA.cf.power=numeric(0), BMA.cf.ppco=numeric(0),
                BMA.cf.ppcc=numeric(0), BMA.aic.beta=numeric(0), 
                BMA.aic.se=numeric(0), BMA.aic.power=numeric(0) )

for(i in 1:nrow(param)){
  print(i)
  beta.Q.cases <- param[i,1]
  beta.Q.controls <- param[i,2]
  beta.L.cases <- param[i,3]
  beta.L.controls <- param[i,4]

  results <- lapply(1:numSims, runSim)
  #Extract pieces of results by method
  case.only <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$case.only }))
  case.control <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$case.control["Y",] }))
  control.only <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$control.only }))
  logistic.cc <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$logistic.cc }))
  BMA.cf <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$BMA.cf }))
  BMA.aic <- do.call(rbind, lapply(1:numSims,function(v){results[[v]]$BMA.aic }))
  
  #Bind all results together (By column)
  OverallResults <- rbind(OverallResults, as.data.frame( cbind(beta.Q.cases, beta.Q.controls, beta.L.cases, beta.L.controls, #Simulation Inputs
    mean(case.only[,1]), mean(case.only[,2]), sum(case.only[,4]<0.05)/numSims,
    mean(case.control[,1]), mean(case.control[,2]), sum(case.control[,4]<0.05)/numSims,
    mean(control.only[,1]), mean(control.only[,2]), sum(control.only[,4]<0.05)/numSims,
    mean(logistic.cc[,1]), mean(logistic.cc[,2]), sum(logistic.cc[,4]<0.05)/numSims,
    mean(BMA.cf[,1]), mean(BMA.cf[,2]), sum(BMA.cf[,4]<0.05)/numSims,
    mean(BMA.cf[,5]),mean(BMA.cf[,6]),
    mean(BMA.aic[,1]), mean(BMA.aic[2]), sum(BMA.aic[,4]<0.05)/numSims,
    mean(BMA.aic[,5]),mean(BMA.aic[,6])) ) )
  
}
names(OverallResults) <- c("beta.Q.cases", "beta.Q.controls", "beta.L.cases", "beta.L.controls",
                           "case.only.beta", "case.only.se", "case.only.power",
                           "case.control.beta", "case.control.se", "case.control.power",
                           "control.only.beta", "control.only.se", "control.only.power",
                           "logistic.cc.beta", "logistic.cc.se", "logistic.cc.power",
                           "BMA.cf.beta", "BMA.cf.se", "BMA.cf.power","BMA.cf.PrMGivenD1","BMA.cf.PrMGivenD2",
                           "BMA.aic.beta", "BMA.aic.se", "BMA.aic.power","BMA.aic.PrMGivenD1","BMA.aic.PrMGivenD2")
#Print Results
OverallResults
t <- Sys.Date()
write.csv(OverallResults,
          paste("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Admixture_Project/Simulation_Results/11.14.2016/",t,"_Results.csv"),row.names=F)

#Reshape data (Power,CF POsterior, AIC POsterior)
#Power: CF
mdat <- OverallResults[c("beta.L.controls","case.only.power","case.control.power",
                         "BMA.cf.power")]
dat.power <- melt(mdat,id="beta.L.controls")
names(dat.power) <- c("beta.L.controls","Model","Power")
dat.power$Model <- as.character(dat.power$Model)
dat.power[dat.power=="case.only.power"] <- "Case-Only"
dat.power[dat.power=="case.control.power"] <- "Case-Control"
dat.power[dat.power=="BMA.cf.power"] <- "BMA"
dat.power <- dat.power[order(dat.power$beta.L.controls),]
dat.power$Model <- as.factor(dat.power$Model)
dat.power$Model <- relevel(dat.power$Model,"Case-Only")

#Power: AIC
mdat.power.aic <- OverallResults[c("beta.L.controls","case.only.power","case.control.power",
                         "BMA.aic.power")]
dat.power.aic <- melt(mdat.power.aic,id="beta.L.controls")
names(dat.power.aic) <- c("beta.L.controls","Model","Power")
dat.power.aic$Model <- as.character(dat.power.aic$Model)
dat.power.aic[dat.power.aic=="case.only.power"] <- "Case-Only"
dat.power.aic[dat.power.aic=="case.control.power"] <- "Case-Control"
dat.power.aic[dat.power.aic=="BMA.cf.power"] <- "BMA"
dat.power.aic <- dat.power[order(dat.power.aic$beta.L.controls),]
dat.power.aic$Model <- as.factor(dat.power.aic$Model)
dat.power.aic$Model <- relevel(dat.power.aic$Model,"Case-Only")

#CF
mdat.cf <- OverallResults[c("beta.L.controls","BMA.cf.PrMGivenD1","BMA.cf.PrMGivenD2")]
dat.cf <- melt(mdat.cf,id="beta.L.controls")
names(dat.cf) <- c("beta.L.controls","Model","PrMGivenD")
dat.cf$Model <- as.character(dat.cf$Model)
dat.cf[dat.cf=="BMA.cf.PrMGivenD1"] <- "Case-Only"
dat.cf[dat.cf=="BMA.cf.PrMGivenD2"] <- "Case-Control"
#AIC
mdat.aic <- OverallResults[c("beta.L.controls","BMA.aic.PrMGivenD1","BMA.aic.PrMGivenD2")]
dat.aic <- melt(mdat.aic,id="beta.L.controls")
names(dat.aic) <- c("beta.L.controls","Model","PrMGivenD")
dat.aic$Model <- as.character(dat.aic$Model)
dat.aic[dat.aic=="BMA.aic.PrMGivenD1"] <- "Case-Only"
dat.aic[dat.aic=="BMA.aic.PrMGivenD2"] <- "Case-Control"

#Plot
pdf(paste("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Admixture_Project/Simulation_Results/11.14.2016/",t,"_CF_Results.pdf"),width=18)
#Power
p1 <- ggplot(dat.power, aes(x = beta.L.controls, y = Power,fill=Model)) +
  geom_bar(stat='identity') + 
  ggtitle("Power, Beta_L.Cases=0.03 ")

p1 <- ggplot(dat.power.aic, aes(x=beta.L.controls,y=Power,fill=Model)) +
  geom_bar(stat="identity",position="dodge") +
  ggtitle("Power, Beta_L.Cases=0.03, Closed-form ")

#BMA.cf (Posterior)
p2 <- ggplot(dat.cf, aes(x = beta.L.controls, y = PrMGivenD,fill=Model)) +
  geom_bar(stat='identity') + 
  ggtitle("BMA.cf Posterior Probability, Beta_L.Cases=0.03")
#BMA.aic (Posterior)
p3 <- ggplot(dat.aic, aes(x = beta.L.controls, y = PrMGivenD,fill=Model)) +
  geom_bar(stat='identity') + 
  ggtitle("BMA.aic Posterior Probability, Beta_L.Cases=0.03")

grid.arrange(p1, p2, p3, ncol=3)
dev.off()






#Plot Results (Old)
pdf(paste("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Admixture_Project/Simulation_Results/11.14.2016/",t,"_Results.pdf"),width=12)
plot(OverallResults$BMA.cf.power~OverallResults$beta.L.cases,
     type="o",pch=".", col="red",lwd=1.8,	
     main = expression(paste(beta,"_L(controls) = ", beta,"_Q(controls) = ",
                             beta,"_Q(cases) = 0")),
     xlab = expression(paste(beta,"_L(cases)")),
     ylab = "Empirical Power",
     xlim=c(-0.08,0.08),ylim=c(0,1),
     xaxt="n")
lines(OverallResults$case.only.power~OverallResults$beta.L.cases,
      type="o", pch=".", lty=2, col="blue",lwd=1.8)

lines(OverallResults$case.control.power~OverallResults$beta.L.cases,
      type="o", pch=".", lty=3, col="black",lwd=1.8)

lines(OverallResults$control.only.power~OverallResults$beta.L.cases,
      type="o", pch=".", lty=4, col="green",lwd=1.8)

lines(OverallResults$BMA.aic.power~OverallResults$beta.L.cases,
      type="o", pch=".", lty=5, col=118,lwd=1.8)

legend("top", c("BMA","Case-Only","Case-Control","Control-Only","BMA.aic"), 
       lty=c(1,2,3,4,5),cex=0.8,
       col=c("red","blue","black","green",118))
axis(1, at=OverallResults$beta.L.cases,labels=round(OverallResults$beta.L.cases,digits=2),
     col.axis="black", las=2, cex.axis=0.7, tck=-.01)
dev.off()
