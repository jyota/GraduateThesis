require(data.table)
require(ggplot2)
require(reshape)
require(scales)
require(matrixStats)
require(grid)
require(limma)
setwd("~/data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/")
# Load all the TCGA files into one matrixx and transform it for use.
listFiles = list.files(pattern="*rsem_gene.txt")
tmpFileMat = lapply(listFiles,read.table,sep='\t',header=T,check.names=F)
trainingSet = rbindlist(tmpFileMat)
castTrainingSet = trainingSet
castTrainingSet$transcript_id <- NULL
castTrainingSet$scaled_estimate <- NULL
castTrainingSet = cast(castTrainingSet, barcode ~ gene_id)
castTrainingSet$class = 1
castTrainingSet[as.numeric(substring(castTrainingSet$barcode,14,15))>10,]$class = 0
rownames(castTrainingSet) <- castTrainingSet$barcode
castTrainingSet$barcode <- NULL
# Save the data now, those steps above take a while...!
write.table(castTrainingSet,"mungedTrainingSet.txt",sep='\t',col.names=T,row.names=T)

# Now can just read from this following statement... instead of performing the above.
read.table("mungedTrainingSet.txt",sep='\t',header=T,check.names=F,row.names=T)
# Remove genes with 0 reads Aligned to them.
castTrainingSet = castTrainingSet[,colSums(castTrainingSet[,1:NCOL(castTrainingSet)])!=0]
classes = castTrainingSet$class 
castTrainingSet$class <- NULL

# Quartile normalize and log2 + counts per million transform.
transTrainingSet = voom(counts=t(castTrainingSet),normalize.method='quantile',plot=T)
readyTrainingSet = t(transTrainingSet$E)
colnames(readyTrainingSet) <- colnames(castTrainingSet)
rownames(readyTrainingSet) <- rownames(castTrainingSet)
transTrainingSet <- NULL

# Select some random genes for before/after plots.
set.seed(123)
randGeneNames = sample(colnames(readyTrainingSet),9)

# Histograms
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet),castTrainingSet[,which(colnames(castTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(preGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Count") + xlab("Expected Count of Short Reads Aligned to Gene") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet),readyTrainingSet[,which(colnames(readyTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(postGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Count") + xlab("Expected Count of Short Reads Aligned to Gene (transformed)") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 

# Boxplots - randomly choose samples.
randSampleNames = sample(rownames(readyTrainingSet),9)
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,]),castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,],check.names=F))
ggplot(preGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + geom_boxplot() + ylab("Expected Count of Short Reads Aligned to Gene") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) 
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,]),readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,],check.names=F))
ggplot(postGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + geom_boxplot() + ylab("Expected Count of Short Reads Aligned to Gene (transformed)") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) 

# Mean vs. Variance plots
options(scipen=999)
prefiltMeanVar = data.frame(Mean=colMeans(castTrainingSet),Variance=colVars(castTrainingSet))
ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(limits=c(0,30000000000),labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
postfiltMeanVar = data.frame(Mean=colMeans(readyTrainingSet),Variance=colVars(readyTrainingSet))
ggplot(postfiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()

# Clear out some of the data from memory prior to feature selection.
prefiltMeanVar <- NULL
postfilMeanVar <- NULL
castTrainingSet <- NULL
trainingSet <- NULL


# Getting into the feature selection steps.
setwd("~/data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/FeatureSelection/")
source('modifiedBagging.R')

# Generate results for biomarkers of size p = 2 through 10. With 80% proportion in training set, 20% in out of bag validation.
P10Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=10,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p10result.txt",sep='\t',col.names=T,row.names=T)
P9Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=9,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p9result.txt",sep='\t',col.names=T,row.names=T)
P8Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=8,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p8result.txt",sep='\t',col.names=T,row.names=T)
P7Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=7,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p7result.txt",sep='\t',col.names=T,row.names=T)
P6Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=6,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p6result.txt",sep='\t',col.names=T,row.names=T)
P5Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=5,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p5result.txt",sep='\t',col.names=T,row.names=T)
P4Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=4,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p4result.txt",sep='\t',col.names=T,row.names=T)
P3Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=3,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p3result.txt",sep='\t',col.names=T,row.names=T)
P2Result = modifiedBagging(as.matrix(readyTrainingSet),classes,rep=100,stopP=2,stopT2=1000,proportion=.8)
write.table(P10Result,"~/Thesis/p2result.txt",sep='\t',col.names=T,row.names=T)

