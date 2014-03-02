require(data.table)
require(ggplot2)
require(reshape)
require(scales)
require(matrixStats)
require(grid)
require(DESeq2)
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
castTrainingSet <- read.table("~/Thesis/mungedTrainingSet.txt",sep='\t',header=T,check.names=F,row.names=1)

# Revising the data preparation... with Upper quartile normalization then log2 transform. Filtering those with < 50 
classes <- castTrainingSet$class 
castTrainingSet$class <- NULL
source('~/Thesis/upperQuartileNormalize.R')
vTraining <- t(upperQuartileNormalize(t(castTrainingSet)))
checkGenes <- data.frame(genes=colnames(vTraining),tumorPresentPct=rep(0,length(colnames(vTraining))),healthyPresentPct=rep(0,length(colnames(vTraining))))
calculateGenes <- data.frame(vTraining>50,check.names=F)
calculateGenes$class <- classes
for(i in 1:NROW(checkGenes)){
	checkGenes[i,1] <- colnames(calculateGenes)[i]
	checkGenes[i,2] <- sum(calculateGenes[calculateGenes$class==1,i])/length(classes[classes==1])
	checkGenes[i,3] <- sum(calculateGenes[calculateGenes$class==0,i])/length(classes[classes==0])
}
genesToRemove <- checkGenes[checkGenes$tumorPresentPct < .25 & checkGenes$healthyPresentPct < .25,]$genes
interTrainingSet <- vTraining[,-which(colnames(vTraining) %in% genesToRemove)]
readyTrainingSet <- log2(interTrainingSet+1.0)

# Check out DESeq2's rlog transform.. not sure how to describe it so probably won't actually use for this thesis anyway.
q = DESeqDataSetFromMatrix(countData=as.matrix(t(checkTrainingSet)),colData=data.frame(class=classes),design=~ class)
rlogq = rlogTransformation(q,blind=T)

# Select some random genes for before/after plots.
set.seed(323)
randGeneNames = sample(colnames(readyTrainingSet),9)

# Histograms
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet),castTrainingSet[,which(colnames(castTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(preGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Count") + xlab("Expected Count of Short Reads Derived from Gene") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet),readyTrainingSet[,which(colnames(readyTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(postGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Count") + xlab("Expected Count of Short Reads Derived from Gene (transformed)") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 

# Boxplots - randomly choose samples.
randSampleNames = sample(rownames(readyTrainingSet),9)
# Just look at those after filtering... but before/after transformation with limma voom.
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]),castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))],check.names=F))
ggplot(preGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + geom_boxplot() + ylab("Expected Count of Short Reads Derived from Gene") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) 
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,]),readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,],check.names=F))
ggplot(postGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + stat_boxplot(geom='errorbar') + geom_boxplot() + ylab("Expected Count of Short Reads Derived from Gene (transformed)") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) 

# Mean vs. Variance plots (before and after limma voom, excluding genes filtered out for all plots)
options(scipen=999)
prefiltMeanVar = data.frame(Mean=colMeans(castTrainingSet[,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]),Variance=colVars(castTrainingSet[,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]))
ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(limits=c(0,30000000000),labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
postfiltMeanVar = data.frame(Mean=colMeans(readyTrainingSet),Variance=colVars(readyTrainingSet))
ggplot(postfiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()

# Clear out some of the data from memory prior to feature selection.
prefiltMeanVar <- NULL
postfilMeanVar <- NULL
castTrainingSet <- NULL
trainingSet <- NULL


# Interact with NCBI's SRA db to get metadata for independent test data set.
require(SRAdb)
sqlfile <- getSRAdbFile()
sraCon <- dbConnect(SQLite(),sqlfile)
# Get all accession identifiers for the independent test data set.
conversion <- sraConvert(c('ERP001058'),sra_con=sraCon)
conversion$description <- "placeholder"
for(i in 1:NROW(conversion)){
	rs <- dbGetQuery(sraCon,paste("SELECT * FROM sample WHERE sample_accession='",conversion[i,]$sample,"'",sep=''))
	conversion[i,]$description <- rs$description
}

# Remove runs with unreleased data.
conversion <- conversion[!(conversion$run %in% c('ERR058695','ERR318894','ERR318891','ERR318895','ERR318892','ERR318893')),]
# Class is 1 for tumor, 0 for healthy/normal
conversion$class <- 1
# Parse description to mark the normal biological samples.
conversion[grepl('apparently_normal',conversion$description),]$class <- 0
# These results will eventually be joined to the filenames of independent data set to get biological class indicator.
# Will need to resolve how to handle more than one run for a given sample...
