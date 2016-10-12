###################################################
# Simulate, analyze, and summarize BMA admixture   
# method.
###################################################
library(data.table)







#Structure and load data (Skip if using simulated data) 
d.map <- read.table("/export/home/mec/MEC_GWAS/AAPC1_AAPC2_typed/local_ancestry_byRFMix/AAMap/AAmap.chr8.liftnew.GeneticPos", header=T, sep=" ")
Q <- read.table("/export/home/mec/MEC_GWAS/AAPC1_AAPC2_typed/local_ancestry_byRFMix/EUR_AFR_asRef/AAPC_6806_chr1-22_EM0_TrioPhased.0.ForwardBackward.txt.AFR.scoreSum.GlobalAncestry", header=F, sep=" ")
Q <- as.numeric(Q[1,])
d.fam <- read.table("/export/home/mec/MEC_GWAS/AAPC1_AAPC2_typed/typed_with_alleles.fam", header=F, sep=" ")
d.fam <- d.fam[1:length(Q),]
Y <- d.fam[,6]-1 
Y <- ifelse(Y==-10, NA, Y)
Z <- ifelse(Y==1, 0, ifelse(Y==0,1, NA))

for(chr 1:22) {
    #chr <- 8
    system.time(d <- fread(paste("AAPC_6806_chr", chr, "_EM0_TrioPhased.0.ForwardBackward.txt.AFR", sep="")))
    
    #case-control linear analysis
    system.time(r <- apply(d, 1, FUN=function(L) { write.table(t(summary(lm(Q ~ Y + offset((as.numeric(L)/2))))$coef[2,]), paste("Chr", chr, "RegResults.txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F) }))
    
    # case-control logistic analysis
    system.time(r <- apply(d, 1, FUN=function(L) { 
    	L.adj <- ((as.numeric(L)/2)-Q)
    	write.table(t(summary(glm(Y ~ Q + L.adj, family="binomial"))$coef[3,]), paste("Chr", chr, "LogisticRegResults.txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F) }))
    
    # case-only analysis
    case.s <- ifelse(is.na(Y), F, ifelse(Y==1, T,F))
    Q.cases <- Q[case.s]
    
    system.time(r <- apply(subset(d, select=(case.s)), 1, FUN=function(L) { write.table(t(summary(lm(Q.cases ~ 1 + offset((as.numeric(L)/2))))$coef[1,]), paste("Chr", chr, "CaseOnlyRegResults.txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F) }))
    
    # case-only analysis with cases and controls
    
    system.time(r <- apply(d, 1, FUN=function(L) { 
    	reg <- lm(Q ~ -1 + offset((as.numeric(L)/2))+Y +Z)
    	s.reg <- summary(reg)$coef
    	write.table(t(c(s.reg[1,], s.reg[2,])), paste("Chr", chr, "CaseOnlyCaseControlsRegResults.txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F) }))
    
    
    # control-only analysis
    control.s <- ifelse(is.na(Y), F, ifelse(Y==0, T,F))
    Q.control <- Q[control.s]
    
    system.time(r <- apply(subset(d, select=(control.s)), 1, FUN=function(L) { write.table(t(summary(lm(Q.control ~ 1 + offset((as.numeric(L)/2))))$coef[1,]), paste("Chr", chr, "ControlOnlyRegResults.txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F) }))

}
