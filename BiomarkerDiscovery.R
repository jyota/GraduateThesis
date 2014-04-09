# Getting into the feature selection steps.
setwd("~/data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/FeatureSelection/")
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

# Check discriminatory space of choice of various P
lda10 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P10Result[1,6:15])))],check.names=F),prior=c(.5,.5))
lda9 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P9Result[1,6:14])))],check.names=F),prior=c(.5,.5))
lda8 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P8Result[2,6:13])))],check.names=F),prior=c(.5,.5))
lda7 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P7Result[1,6:12])))],check.names=F),prior=c(.5,.5))
lda6 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P6Result[1,6:11])))],check.names=F),prior=c(.5,.5))
lda5 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P5Result[1,6:10])))],check.names=F),prior=c(.5,.5))
lda4 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P4Result[1,6:9])))],check.names=F),prior=c(.5,.5))
lda3 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P3Result[1,6:8])))],check.names=F),prior=c(.5,.5))
lda2 <- lda(class ~ .,data=data.frame(class=classes,readyTrainingSet[,which(colnames(readyTrainingSet) %in% as.character(unlist(P2Result[1,6:7])))],check.names=F),prior=c(.5,.5))
ldaProbs <- rbind(data.frame(probability=predict(lda10,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 10',452),class=classes),
	data.frame(probability=predict(lda9,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 9',452),class=classes),
	data.frame(probability=predict(lda8,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 8',452),class=classes),
	data.frame(probability=predict(lda7,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 7',452),class=classes),
	data.frame(probability=predict(lda6,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 6',452),class=classes),
	data.frame(probability=predict(lda5,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 5',452),class=classes),
	data.frame(probability=predict(lda4,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 4',452),class=classes),
	data.frame(probability=predict(lda3,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 3',452),class=classes),
	data.frame(probability=predict(lda2,as.data.frame(readyTrainingSet,check.names=F))$posterior[,1],P=rep('p = 2',452),class=classes))

ldaProbs$Class <- 'Tumor'
ldaProbs[ldaProbs$class==0,]$Class <- 'Healthy'
ggplot(ldaProbs,aes(x=probability,y=0,colour=Class))+geom_point(alpha=0.5, shape=5)+facet_wrap(~ P,ncol=2)+theme_bw()+scale_y_continuous(breaks=c(0))+xlab('Posterior Probability of Healthy Class')+ylab('')+theme(panel.margin=unit(2,'mm'))+scale_colour_manual(values=c('blue','red'))

# Find the informative set of genes
source('findInformative.R')
informativeSet = findInformative(x=as.matrix(readyTrainingSet),y=classes,rep=300,proportion=.8,stopP=5,stopT2=1000,priors=c(.5,.5))
write.table(informativeSet,"~/Thesis/informativeSetDetermine.txt",sep='\t',col.names=T,row.names=T)

# If needed to be read back in:
informativeSet <- read.table("~/Thesis/informativeSetDetermine.txt",sep='\t',header=T,row.names=1)

informativeSet$T2 <- as.numeric(as.character(informativeSet$T2))
# Find a candidate level for T2 to cutoff at. We'll go from anything not associated with T2 > 2.0 and increment cutoff each run.
# Obtain accuracy estimates from modified bagging for genes not above each cutoff to help make decision.
infSetTwo <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetTwo$repStats,"~/Thesis/infSetTwo.txt",sep='\t',row.names=T,col.names=T)
infSetThree <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetThree$repStats,"~/Thesis/infSetThree.txt",sep='\t',row.names=T,col.names=T)
infSetFour <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetFour$repStats,"~/Thesis/infSetFour.txt",sep='\t',row.names=T,col.names=T)
infSetFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>5.0,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetFive$repStats,"~/Thesis/infSetFive.txt",sep='\t',row.names=T,col.names=T)
# Check out more granular cutoffs to see what may be appropriate.
infSetThreePointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.5,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetThreePointFive$repStats,"~/Thesis/infSetThreePointFive.txt",sep='\t',row.names=T,col.names=T)
infSetTwoPointFive <- modifiedBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>2.5,3:7]))]),classes,rep=100,stopP=5,stopT2=1000,proportion=.8,progressBar=T,priors=c(.5,.5))
write.table(infSetTwoPointFive$repStats,"~/Thesis/infSetTwoPointFive.txt",sep='\t',row.names=T,col.names=T)


# Plot informative set, with line for our cutoff
informativeSet$Index <- as.numeric(as.character(informativeSet$Index))
ggplot(informativeSet,aes(x=Index,y=T2))+geom_point()+theme_bw()+geom_hline(yintercept=3.0,colour='red',alpha=.8)

source('findInformativeBagging.R')
# Get estimates for final informative set of genes now that cutoff decided (1000 modified bagging schema iterations). Store associated variables this time.
infSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(infSetFinal,"~/Thesis/infSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Show 'discriminatory space' for LDA classifier from informative set via posterior probabilities of class.
infBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,stopP=5,stopT2=1000)
# Fit an LDA model with biomarker variables
infFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(infBiomarker))]),class=classes,check.names=F),prior=c(.5,.5))
infPredFit <- predict(infFit, data.frame(readyTrainingSet,check.names=F))
infPlotFit <- data.frame(classProb=c(infPredFit$posterior[,1],infPredFit$posterior[,2]),
	shown_class=c(rep('Healthy',length(classes)),rep('Tumor',length(classes))),actual_class=c(classes,classes))
infPlotFit$Class <- 'Healthy'
infPlotFit[infPlotFit$actual_class==1,]$Class <- 'Tumor'
ggplot(infPlotFit,aes(classProb,fill=Class))+geom_histogram(binwidth=0.02)+facet_grid(shown_class ~ .) + theme_bw()+xlab('Posterior Probability of Class Membership')

# Check these estimates for entire set of variables (1000 modified bagging schema iterations)
fullSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(fullSetFinal,"~/Thesis/fullSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Get modified bagging estimates for 'non-informative' set of genes
nonInfSetFinal <- findInformativeBagging(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,rep=1000,stopP=5,stopT2=1000,proportion=.8,priors=c(.5,.5))
write.table(nonInfSetFinal$repStats,"~/Thesis/nonInfSetFinal.txt",sep='\t',row.names=T,col.names=T)

# Show 'discriminatory space' for LDA classifier NOT from informative set
infBiomarker <- hybridFeatureSelection(as.matrix(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>3.0,3:7]))]),classes,stopP=5,stopT2=1000)
# Fit an LDA model with biomarker variables
infFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(infBiomarker))]),class=classes,check.names=F),prior=c(.5,.5))
infPredFit <- predict(infFit, data.frame(readyTrainingSet,check.names=F))
infPlotFit <- data.frame(classProb=c(infPredFit$posterior[,1],infPredFit$posterior[,2]),
	shown_class=c(rep('Healthy',length(classes)),rep('Tumor',length(classes))),actual_class=c(classes,classes))
infPlotFit$Class <- 'Healthy'
infPlotFit[infPlotFit$actual_class==1,]$Class <- 'Tumor'
ggplot(infPlotFit,aes(classProb,fill=Class))+geom_histogram(binwidth=0.02)+facet_grid(shown_class ~ .) + theme_bw()+xlab('Posterior Probability of Class Membership')

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
biomarkerFit <- lda(class ~ .,data=data.frame(as.matrix(readyTrainingSet[,which(colnames(readyTrainingSet) %in% names(finalBiomarker))]),class=classes,check.names=F),prior=c(.5,.5))
# Confusion matrix for training data (just for gut check, not relevant result really)
table(classes,predict(biomarkerFit,priors=c(.5,.5),data.frame(readyTrainingSet,check.names=F))$class)

# Retrieve biomarker genes from model in case they were lost somehow
finalBiomarker <- names(attributes(biomarkerFit$terms)$dataClasses)[2:6]

# Intersection of 1,000 modified bagging iteration full training set perfect classifiers and frequent primary genes.
table(unlist(fullSetFinal[fullSetFinal$Accuracy==1.0,6:10])[unlist(fullSetFinal[fullSetFinal$Accuracy==1.0,6:10]) %in% frequentPrimary])
# Intersection like the above, but for informative set of genes only.
table(unlist(infSetFinal[infSetFinal$Accuracy==1.0,6:10])[unlist(infSetFinal[infSetFinal$Accuracy==1.0,6:10]) %in% frequentPrimary])
