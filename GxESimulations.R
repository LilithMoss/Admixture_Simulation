#Simulation 
simulate.FUN <- function(theta,B.DSL,B.DSLxE){
G.mat <- matrix(ncol=M.DSL,nrow=0)
e <- numeric()
y <- numeric()
k=1 #Counter 1
L=1 #Counter 2

repeat{
#1.  Simulate E:
#a.  E ~ rbinom(N, 1, pE)
E <- rbinom(k*N.gen,1,pE)

#2.  Simulate DSL(s):
#a.	For each DSL {DSL1, ., DSLM(DSL)}: Pr(DSLm=1 | E) = expit[]
#etaDSL <- log(pDSL/(1-pDSL)) + theta*(E-mean(E))
#etaDSL <- log(pDSL/(1-pDSL)) + theta*E
etaDSL <- do.call(cbind,lapply(1:M.DSL,function(x) log(pDSL[x]/(1-pDSL[x])) + theta[x]*(E-mean(E) )))
pDSL.E <- expit(etaDSL)
#b.  incorporates DSL-E association: DSL = (N.gen x M.DSL) Matrix
#DSL <- matrix( rbinom(k*N.gen*M.DSL,1,pDSL.E),nrow=N.gen,ncol=M.DSL  )
#DSL <- matrix( rbinom(k*N.gen*M.DSL,1,pDSL.E),nrow=k*N.gen,ncol=M.DSL  )
DSL.raw <- matrix( rbinom(k*N.gen*M.DSL,2,pDSL.E),nrow=k*N.gen,ncol=M.DSL  )
DSL <- ifelse(DSL.raw < 1,0,1) 

#3.  Simulate Y:
#a.	Pr(Y | E, DSL(s)) = expit[] : (N.gen x M.DSL) Matrix)
#DSL.mat <- B.DSL*( DSL - rep(colMeans(DSL),rep.int(nrow(DSL),ncol(DSL))) )
DSL.centered <- as.matrix( DSL - rep(colMeans(DSL),rep.int(nrow(DSL),ncol(DSL)))  )
E.centered <- as.matrix( E-mean(E) )
#DSL.centered <- DSL #For now, call them centered, but don't actually center them
#E.centered <- as.matrix(E)  #For now, call them centered, but don't actually center them

#etaY = (k*N.gen x 1) matrix Or do I want (1 x N.gen) matrix? This gives (k*N.gen x 1):
#etaY <- rep(log(pY/(1-pY)),k*N.gen) + B.E*E.centered + B.DSL*rowSums(DSL.centered) + B.DSLxE*E.centered*rowSums(DSL.centered)
#pY.E.DSL <- expit(etaY)
etaY <- rep(log(pY/(1-pY)),k*N.gen) + B.E*E.centered + 
        rowSums( do.call(cbind,lapply(1:M.DSL,function(x) B.DSL[x]*DSL.centered[,x])) ) + 
        rowSums(do.call(cbind,lapply(1:M.DSL,function(x) B.DSLxE[x]*E.centered*DSL.centered[,x])) ) 
pY.E.DSL <- expit(etaY)

#b.	incorporates DSL-Y main effects
Y <- matrix( rbinom(k*N.gen,1,pY.E.DSL),nrow=k*N.gen, ncol=1 )
GENOTYPE <- DSL #Or GENOTYPE <- cbind(DSL,G.C) 

#Bind to existing matrices
G.mat <- rbind(G.mat,GENOTYPE)
e <- c(e,E)
y <- c(y,Y)
k=0.1
L=L+1
if (sum(y) > N.cases)
  break
}

#6.  Case-control sampling
#a.	Sample N.cases (Y=1) and N.controls (Y=0)
prelim<-as.data.frame(cbind(y,e,G.mat))
names(prelim) <- c("Y","E",rep("DSL",1),rep("G.E.DSL",M.G))
cases <- prelim[prelim$Y==1,]
controls <- prelim[prelim$Y==0,]
cases <- cases[sample(nrow(cases), N.cases), ]
controls <- controls[sample(nrow(controls), N.controls), ]
Sample <- rbind(cases,controls)
return(Sample)
} 





