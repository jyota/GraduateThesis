# Getting into the feature selection steps.
setwd("~/FeatureSelection/")
source('findInformativeBagging.R')
source('modifiedBagging.R')

# Generate results for biomarkers of size p = 2 through 10. With 80% proportion in training set, 20% in out of bag validation.
P10Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=10,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P10Result,"~/Thesis/p10result.txt",sep='\t',col.names=T,row.names=T)
P9Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=9,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P9Result,"~/Thesis/p9result.txt",sep='\t',col.names=T,row.names=T)
P8Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P8Result,"~/Thesis/p8result.txt",sep='\t',col.names=T,row.names=T)
P7Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=7,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P7Result,"~/Thesis/p7result.txt",sep='\t',col.names=T,row.names=T)
P6Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=6,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P6Result,"~/Thesis/p6result.txt",sep='\t',col.names=T,row.names=T)
P5Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=5,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P5Result,"~/Thesis/p5result.txt",sep='\t',col.names=T,row.names=T)
P4Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=4,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P4Result,"~/Thesis/p4result.txt",sep='\t',col.names=T,row.names=T)
P3Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=3,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P3Result,"~/Thesis/p3result.txt",sep='\t',col.names=T,row.names=T)
P2Result = findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=500,stopP=2,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(P2Result,"~/Thesis/p2result.txt",sep='\t',col.names=T,row.names=T)

# Reload data if needed
P10Result <- read.table('~/Thesis/p10result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P9Result <- read.table('~/Thesis/p9result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P8Result <- read.table('~/Thesis/p8result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P7Result <- read.table('~/Thesis/p7result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P6Result <- read.table('~/Thesis/p6result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P5Result <- read.table('~/Thesis/p5result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P4Result <- read.table('~/Thesis/p4result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P3Result <- read.table('~/Thesis/p3result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)
P2Result <- read.table('~/Thesis/p2result.txt',sep='\t',header=T,row.names=1,stringsAsFactors=F)


# Plots of accuracy, sensitivity, and specificity for various P
repAcc <- data.frame(accuracy=c(as.numeric(as.character(P2Result$Accuracy)),
	as.numeric(as.character(P3Result$Accuracy)),
	as.numeric(as.character(P4Result$Accuracy)),
	as.numeric(as.character(P5Result$Accuracy)),
	as.numeric(as.character(P6Result$Accuracy)),
	as.numeric(as.character(P7Result$Accuracy)),
	as.numeric(as.character(P8Result$Accuracy)),
	as.numeric(as.character(P9Result$Accuracy)),
	as.numeric(as.character(P10Result$Accuracy))),
    P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repAcc$P <- as.factor(repAcc$P)
ggplot(repAcc,aes(x=accuracy,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()
repSens <- data.frame(sensitivity=c(as.numeric(as.character(P2Result$Sensitivity)),
	as.numeric(as.character(P3Result$Sensitivity)),
	as.numeric(as.character(P4Result$Sensitivity)),
	as.numeric(as.character(P5Result$Sensitivity)),
	as.numeric(as.character(P6Result$Sensitivity)),
	as.numeric(as.character(P7Result$Sensitivity)),
	as.numeric(as.character(P8Result$Sensitivity)),
	as.numeric(as.character(P9Result$Sensitivity)),
	as.numeric(as.character(P10Result$Sensitivity))),
	P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repSens$P <- as.factor(repSens$P)
ggplot(repSens,aes(x=sensitivity,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()
repSpec <- data.frame(specificity=c(as.numeric(as.character(P2Result$Specificity)),
	as.numeric(as.character(P3Result$Specificity)),
	as.numeric(as.character(P4Result$Specificity)),
	as.numeric(as.character(P5Result$Specificity)),
	as.numeric(as.character(P6Result$Specificity)),
	as.numeric(as.character(P7Result$Specificity)),
	as.numeric(as.character(P8Result$Specificity)),
	as.numeric(as.character(P9Result$Specificity)),
	as.numeric(as.character(P10Result$Specificity))),
 	P=c(rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)))
repSpec$P <- as.factor(repSpec$P)
ggplot(repSpec,aes(x=specificity,colour=P))+geom_density(size=1.1,alpha=.8)+theme_classic()

# Find the informative set of genes
source('findInformative.R')
informativeSet = findInformative(x=as.matrix(readyTrainingSet),y=classes,rep=500,proportion=.8,stopP=8,stopT2=1000,priors=c(.5,.5))
write.table(informativeSet,"~/Thesis/informativeSetDetermine.txt",sep='\t',col.names=T,row.names=T)

# If needed to be read back in:
#informativeSet <- read.table("~/Thesis/informativeSetDetermine.txt",sep='\t',header=T,row.names=1)

informativeSet$T2 <- as.numeric(as.character(informativeSet$T2))
# Find a candidate level for T2 to cutoff at. We'll go from anything not associated with T2 > 2.0 and increment cutoff each run.
# Obtain accuracy estimates from modified bagging for genes not above each cutoff to help make decision.
infSetThree <- modifiedBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.0,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetThree$repStats,"~/Thesis/infSetThree.txt",sep='\t',row.names=T,col.names=T)
infSetNotThree <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.0,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetNotThree$repStats,"~/Thesis/infSetNotThree.txt",sep='\t',row.names=T,col.names=T)
infSetFour <- modifiedBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=4.0,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetFour$repStats,"~/Thesis/infSetFour.txt",sep='\t',row.names=T,col.names=T)
infSetNotFour <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=4.0,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetNotFour$repStats,"~/Thesis/infSetNotFour.txt",sep='\t',row.names=T,col.names=T)
infSetThreePointFive <- modifiedBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.5,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetThreePointFive$repStats,"~/Thesis/infSetThreePointFive.txt",sep='\t',row.names=T,col.names=T)
infSetNotThreePointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.5,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetNotThreePointFive$repStats,"~/Thesis/infSetNotThreePointFive.txt",sep='\t',row.names=T,col.names=T)
infSetTwoPointFive <- modifiedBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=2.5,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetTwoPointFive$repStats,"~/Thesis/infSetTwoPointFive.txt",sep='\t',row.names=T,col.names=T)
infSetNotTwoPointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=2.5,3:10]))]),classes,rep=500,stopP=8,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetNotTwoPointFive$repStats,"~/Thesis/infSetNotTwoPointFive.txt",sep='\t',row.names=T,col.names=T)


# Plot informative set, with line for our cutoff
informativeSet$Index <- as.numeric(as.character(informativeSet$Index))
ggplot(informativeSet,aes(x=Index,y=T2))+geom_point()+theme_bw()+geom_hline(yintercept=3.0,colour='red',alpha=.8)+scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8))

source('findInformativeBagging.R')
# Get estimates for final informative set of genes now that cutoff decided (1000 modified bagging schema iterations). Store associated variables this time.
infSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.0,3:10]))]),classes,rep=1000,stopP=8,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(infSetFinal,"~/Thesis/infSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Check these estimates for entire set of variables (1000 modified bagging schema iterations)
fullSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=1000,stopP=8,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(fullSetFinal,"~/Thesis/fullSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Get modified bagging estimates for 'non-informative' set of genes
nonInfSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.0,3:10]))]),classes,rep=1000,stopP=8,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(nonInfSetFinal,"~/Thesis/nonInfSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Cluster by Pearson's correlation distance 
require(Hmisc)
clustResult <- varclus(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>=3.0,3:10]))],similarity='pearson',method='complete')
memb <- cutree(clustResult$hclust,k=20)
clusteringResult <- data.frame(cluster=seq(1:20),size=rep(0,20),use=rep(0,20),avg_use=rep(0,20))
for(i in 1:20){
	clusteringResult[i,]$size <- table(memb)[i]
	clusteringResult[i,]$use <- sum(table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:13]))[names(table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:13]))) %in% names(memb[memb==i])])
}
clusteringResult$avg_use <- clusteringResult$use / clusteringResult$size

# Identify frequently used genes -- those that are in at least 1% of OOB classifiers in informative set of genes modified bagging classifiers.
q <- table(unlist(infSetFinal[infSetFinal$Accuracy==1,6:13]))/1000.0
frequentlyUsed <- q[q>=0.01]
# Frequent primaries will be those from frequentlyUsed that also are in clusters with >5% average use
frequentPrimary <- names(frequentlyUsed)[names(frequentlyUsed) %in% names(memb[memb %in% c(clusteringResult[clusteringResult$avg_use>5.0,]$cluster)])]
# These intersect with those used in at least 1% perfect OOB classifiers for full training set.
q <- table(unlist(fullSetFinal[fullSetFinal$Accuracy==1,6:13]))/1000.0
frequentlyUsed <- q[q>=0.01]
mostFrequent <- frequentPrimary[frequentPrimary %in% names(frequentlyUsed)]

finalBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% mostFrequent)]),classes,stopP=8,stopT2=1000)

# Fit an LDA model with biomarker variables
biomarkerFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(finalBiomarker))]),class=classes,check.names=F),prior=c(.5,.5))
# Confusion matrix for training data (just for gut check, not relevant result really)
table(classes,predict(biomarkerFit,prior=c(.5,.5),data.frame(readyTrainingSet,check.names=F))$class)

# Retrieve biomarker genes from model in case they were lost somehow
finalBiomarker <- names(attributes(biomarkerFit$terms)$dataClasses)[2:9]

# Intersection of 1,000 modified bagging iteration full training set perfect classifiers and frequent primary genes.
table(unlist(fullSetFinal[fullSetFinal$Accuracy==1.0,6:13])[unlist(fullSetFinal[fullSetFinal$Accuracy==1.0,6:13]) %in% frequentPrimary])
# Intersection like the above, but for informative set of genes only.
table(unlist(infSetFinal[infSetFinal$Accuracy==1.0,6:13])[unlist(infSetFinal[infSetFinal$Accuracy==1.0,6:13]) %in% frequentPrimary])

# Examine distribution of expression values chosen for biomarker in both test and training set.
q = data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% c(names(finalBiomarker)))],check.names=F)
postGeneExamine <- melt(q, id=c('class'))
postGeneExamine$dataset <- 'Training'
q = data.frame(class=testClasses,readyTestSet[,which(colnames(readyTestSet) %in% c(names(finalBiomarker)))],check.names=F)
postGene2 <- melt(q, id=c('class'))
postGene2$dataset <- 'Test'
combGene <- rbind(postGeneExamine, postGene2)
combGene$Class <- 'Tumor'
combGene[combGene$class==0, ]$Class <- 'Healthy'
ggplot(combGene,aes(x=value, fill=Class)) + geom_density(alpha=.5) + facet_wrap(~ variable + dataset, ncol=2)+theme_bw()+xlab('Normalized & transformed expression value')

# Plot discriminatory space for test set classifier
ret = c(50, 894, 1221, 1672, 6022, 6605, 6929, 7276)
fuse = lda(x = readyTrainingSet[, which(colnames(readyTrainingSet) %in% colnames(readyTrainingSet)[ret])], grouping = classes, prior = c(.5, .5))

res = data.frame(LD = predict(fuse, readyTestSet[, which(colnames(readyTestSet) %in% colnames(readyTrainingSet)[ret])])$x, class = testClasses)

ggplot(res, aes(x = LD1, y = ifelse(class == 0, 0.1, -0.1), colour = ifelse(class == 0, "Healthy", "Tumor"))) + geom_point(shape = 5) + coord_cartesian(ylim = c(-10, 10)) + scale_colour_manual(name = "Sample Class", values = c("blue", "red")) + geom_vline(xintercept = 0, col = "purple", linetype = 2) + coord_cartesian(ylim=c(-0.5, 0.5)) + ylab("") + xlab("Sample linear score (given by R LDA model)") + scale_y_discrete(breaks = NULL) + theme_bw()

