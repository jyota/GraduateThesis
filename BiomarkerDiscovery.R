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
informativeSet = findInformative(x=as.matrix(readyTrainingSet),y=classes,rep=300,proportion=.8,stopP=5,stopT2=1000)
write.table(informativeSet,"~/Thesis/informativeSetDetermine.txt",sep='\t',col.names=T,row.names=T)

# If needed to be read back in:
informativeSet <- read.table("~/Thesis/informativeSetDetermine.txt",sep='\t',header=T,row.names=1)

informativeSet$T2 <- as.numeric(as.character(informativeSet$T2))
# Find a candidate level for T2 to cutoff at. We'll go from anything not associated with T2 > 2.0 and increment cutoff each run.
# Obtain accuracy estimates from modified bagging for genes not above each cutoff to help make decision.
infSetTwo <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetTwo$repStats,"~/Thesis/infSetTwo.txt",sep='\t',row.names=T,col.names=T)
infSetThree <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetThree$repStats,"~/Thesis/infSetThree.txt",sep='\t',row.names=T,col.names=T)
infSetFour <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetFour$repStats,"~/Thesis/infSetFour.txt",sep='\t',row.names=T,col.names=T)
infSetFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>5.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetFive$repStats,"~/Thesis/infSetFive.txt",sep='\t',row.names=T,col.names=T)
# Check out more granular cutoffs to see what may be appropriate.
infSetThreePointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.5,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetThreePointFive$repStats,"~/Thesis/infSetThreePointFive.txt",sep='\t',row.names=T,col.names=T)
infSetTwoPointSFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.75,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetTwoPointSFive$repStats,"~/Thesis/infSetTwoPointSFive.txt",sep='\t',row.names=T,col.names=T)
infSetTwoPointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.5,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetTwoPointFive$repStats,"~/Thesis/infSetTwoPointFive.txt",sep='\t',row.names=T,col.names=T)

# Plot informative set, with line for our cutoff
informativeSet$Index <- as.numeric(as.character(informativeSet$Index))
ggplot(informativeSet,aes(x=Index,y=T2))+geom_point()+theme_bw()+geom_hline(yintercept=4.0,colour='red',alpha=.8)


source('findInformativeBagging.R')
# Get estimates for final informative set of genes now that cutoff decided (1000 modified bagging schema iterations). Store associated variables this time.
infSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8)
write.table(infSetFinal,"~/Thesis/infSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Show 'discriminatory space' for LDA classifier from informative set
infBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,stopP=5,stopT2=1000)
# Fit an LDA model with biomarker variables
infFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(infBiomarker))]),class=classes,check.names=F))
infPredFit <- predict(infFit, data.frame(readyTrainingSet,check.names=F))
infPlotFit <- data.frame(LD1=infPredFit$x,pred_class=infPredFit$class,actual_class=classes)
infPlotFit$Class <- 'Healthy'
infPlotFit[infPlotFit$actual_class==1,]$Class <- 'Tumor'
ggplot(infPlotFit,aes(LD1,fill=Class))+geom_density(alpha=.8)+theme_bw()+xlab('Discriminant Score')

# Check these estimates for entire set of variables (1000 modified bagging schema iterations)
fullSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8)
write.table(fullSetFinal,"~/Thesis/fullSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Get modified bagging estimates for 'non-informative' set of genes
nonInfSetFinal <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8)
write.table(nonInfSetFinal$repStats,"~/Thesis/nonInfSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Show 'discriminatory space' for LDA classifier NOT from informative set
infBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,stopP=5,stopT2=1000)
# Fit an LDA model with biomarker variables
infFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(infBiomarker))]),class=classes,check.names=F))
infPredFit <- predict(infFit, data.frame(readyTrainingSet,check.names=F))
infPlotFit <- data.frame(LD1=infPredFit$x,pred_class=infPredFit$class,actual_class=classes)
infPlotFit$Class <- 'Healthy'
infPlotFit[infPlotFit$actual_class==1,]$Class <- 'Tumor'
ggplot(infPlotFit,aes(LD1,fill=Class))+geom_density(alpha=.8)+theme_bw()+xlab('Discriminant Score')


# Cluster by Pearson's correlation distance 
require(Hmisc)
clustResult <- varclus(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))],similarity='pearson',method='complete')
memb <- cutree(clustResult$hclust,k=10)
clusteringResult <- data.frame(cluster=seq(1:10),size=rep(0,10),use=rep(0,10),avg_use=rep(0,10))
for(i in 1:10){
	clusteringResult[i,]$size <- table(memb)[i]
	clusteringResult[i,]$use <- sum(table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:10]))[names(table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:10]))) %in% names(memb[memb==i])])
}
clusteringResult$avg_use <- clusteringResult$use / clusteringResult$size

# Identify frequently used genes -- those that are in at least 1% of OOB classifiers in informative set of genes modified bagging classifiers.
q <- table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:10]))/1000.0
frequentlyUsed <- q[q>=0.01]
# Frequent primaries will be those from frequentlyUsed that also are in clusters with >=15% average use
frequentPrimary <- names(frequentlyUsed)[names(frequentlyUsed) %in% names(memb[memb %in% c(clusteringResult[clusteringResult$avg_use>15.0,]$cluster)])]
# These intersect with those used in perfect OOB classifiers for full training set.
frequentPrimary[frequentPrimary %in% unlist(fullSetFinal[fullSetFinal$Accuracy==1,6:10])]

finalBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% frequentPrimary)]),classes,stopP=5,stopT2=1000)

# Fit an LDA model with biomarker variables
biomarkerFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(finalBiomarker))]),class=classes,check.names=F))
# Confusion matrix for training data (just for gut check, not relevant result really)
table(classes,predict(biomarkerFit,priors=c(.5,.5),data.frame(readyTrainingSet,check.names=F))$class)
