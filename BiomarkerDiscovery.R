# Getting into the feature selection steps.
setwd("~/data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/FeatureSelection/")
source('modifiedBagging.R')

# Generate results for biomarkers of size p = 2 through 10. With 80% proportion in training set, 20% in out of bag validation.
P10Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=10,stopT2=1000,proportion=.8)
write.table(P10Result$repStats,"~/Thesis/p10result.txt",sep='\t',col.names=T,row.names=T)
P9Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=9,stopT2=1000,proportion=.8)
write.table(P9Result$repStats,"~/Thesis/p9result.txt",sep='\t',col.names=T,row.names=T)
P8Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=8,stopT2=1000,proportion=.8)
write.table(P8Result$repStats,"~/Thesis/p8result.txt",sep='\t',col.names=T,row.names=T)
P7Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=7,stopT2=1000,proportion=.8)
write.table(P7Result$repStats,"~/Thesis/p7result.txt",sep='\t',col.names=T,row.names=T)
P6Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=6,stopT2=1000,proportion=.8)
write.table(P6Result$repStats,"~/Thesis/p6result.txt",sep='\t',col.names=T,row.names=T)
P5Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(P5Result$repStats,"~/Thesis/p5result.txt",sep='\t',col.names=T,row.names=T)
P4Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=4,stopT2=1000,proportion=.8)
write.table(P4Result$repStats,"~/Thesis/p4result.txt",sep='\t',col.names=T,row.names=T)
P3Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=3,stopT2=1000,proportion=.8)
write.table(P3Result$repStats,"~/Thesis/p3result.txt",sep='\t',col.names=T,row.names=T)
P2Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=2,stopT2=1000,proportion=.8)
write.table(P2Result$repStats,"~/Thesis/p2result.txt",sep='\t',col.names=T,row.names=T)

# Plots of accuracy, sensitivity, and specificity for various P
repAcc = data.frame(accuracy=c(P2Result$repStats[,1],P3Result$repStats[,1],P4Result$repStats[,1],P5Result$repStats[,1],P6Result$repStats[,1],P7Result$repStats[,1],P8Result$repStats[,1],P9Result$repStats[,1],P10Result$repStats[,1]),P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repAcc$P = as.factor(repAcc$P)
ggplot(repAcc,aes(x=accuracy,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()
repSens = data.frame(sensitivity=c(P2Result$repStats[,2],P3Result$repStats[,2],P4Result$repStats[,2],P5Result$repStats[,2],P6Result$repStats[,2],P7Result$repStats[,2],P8Result$repStats[,2],P9Result$repStats[,2],P10Result$repStats[,2]),P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repSens$P = as.factor(repSens$P)
ggplot(repSens,aes(x=sensitivity,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()
repSpec = data.frame(specificity=c(P2Result$repStats[,3],P3Result$repStats[,3],P4Result$repStats[,3],P5Result$repStats[,3],P6Result$repStats[,3],P7Result$repStats[,3],P8Result$repStats[,3],P9Result$repStats[,3],P10Result$repStats[,3]),P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repSpec$P = as.factor(repSpec$P)
ggplot(repSpec,aes(x=specificity,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()

# Find the informative set of genes
source('findInformative.R')
informativeSet = findInformative(x=as.matrix(readyTrainingSet),y=classes,rep=1000,proportion=.8,stopP=7,stopT2=1000)
write.table(informativeSet,"~/Thesis/informativeSetDetermine.txt",sep='\t',col.names=T,row.names=T)

# If needed to be read back in:
informativeSet <- read.table("~/Thesis/informativeSetDetermine.txt",sep='\t',header=T,row.names=1)

# Find a candidate level for T2 to cutoff at. We'll go from anything not associated with T2 > 2.0 and increment by 0.5 each run.
# Obtain accuracy estimates from modified bagging for genes not above each cutoff to help make decision.
infSetTwo <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.0,3:9]))]),classes,rep=100,stopP=7,stopT2=1000,proportion=.8)

# Plot when cutoff decided. Need to add lines to indicate cutoff.
ggplot(informativeSet,aes(x=Index,y=T2))+geom_point()+theme_classic()