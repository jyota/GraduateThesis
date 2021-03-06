require(data.table)
require(ggplot2)
require(reshape)
require(scales)
require(matrixStats)
require(grid)
require(edgeR)
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

# Revising the data preparation... with Upper quartile normalization then log2 transform. Filtering based on cpm in <50% of size of smallest biological class.
classes <- castTrainingSet$class 
castTrainingSet$class <- NULL
source('~/Thesis/upperQuartileNormalize.R')
vTraining <- t(upperQuartileNormalize(t(castTrainingSet)))
zTraining <- vTraining[,colSums(vTraining)!=0]
checkGenes <- data.frame(genes=colnames(zTraining),tumorPresent=rep(0,length(colnames(zTraining))),healthyPresent=rep(0,length(colnames(zTraining))))
checkGenes$genes <- as.character(checkGenes$genes)
cutOff <- round((max(apply(cpm(t(zTraining)),2,mean,trim=0.02))+min(apply(cpm(t(zTraining)),2,mean,trim=0.02)))/2.0,0)
calculateGenes <- data.frame(cpm(t(zTraining))>=cutOff,check.names=F)
calculateGenes <- data.frame(t(calculateGenes),check.names=F)
calculateGenes$class <- classes
for(i in 1:NROW(checkGenes)){
	checkGenes[i,1] <- colnames(calculateGenes)[i]
	checkGenes[i,2] <- sum(calculateGenes[calculateGenes$class==1,i])
	checkGenes[i,3] <- sum(calculateGenes[calculateGenes$class==0,i])
}
genesToRemove <- checkGenes[checkGenes$tumorPresent < 0.2*length(classes[classes==1]) & checkGenes$healthyPresent < 0.2*length(classes[classes==0]),]$genes
interTrainingSet <- zTraining[,-which(colnames(zTraining) %in% genesToRemove)]
readyTrainingSet <- log2(interTrainingSet+1.0)

# Select some random genes for before/after plots.
set.seed(323)
randGeneNames = sample(colnames(readyTrainingSet),9)

# Histograms
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet),castTrainingSet[,which(colnames(castTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(preGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Number of Samples per Bin") + xlab("Gene Expression Value") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet),readyTrainingSet[,which(colnames(readyTrainingSet) %in% randGeneNames)],check.names=F))
ggplot(postGeneExamine,aes(value)) + geom_histogram() + facet_wrap(~ variable) + ylab("Number of Samples per Bin") + xlab("Gene Expression Value (transformed)") +  scale_x_continuous(labels=comma) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) 

# Boxplots - randomly choose samples.
set.seed(399)
randSampleNames = sample(rownames(readyTrainingSet),9)
# Just look at those after filtering... but before/after transformation.
preGeneExamine = melt(data.frame(barcode=rownames(castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]),castTrainingSet[rownames(castTrainingSet) %in% randSampleNames,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))],check.names=F))
ggplot(preGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + stat_boxplot(geom = 'errorbar') + geom_boxplot() + ylab("Gene Expression Value") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) + coord_cartesian(ylim=c(0,25000))
postGeneExamine = melt(data.frame(barcode=rownames(readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,]),readyTrainingSet[rownames(readyTrainingSet) %in% randSampleNames,],check.names=F))
ggplot(postGeneExamine,aes(factor(barcode),value,fill=factor(barcode))) + stat_boxplot(geom='errorbar') + geom_boxplot() + ylab("Gene Expression Value (transformed)") + xlab("Sample") + theme_bw() + theme(plot.margin=unit(c(1,1,1,1),'lines'),axis.text.x = element_text(angle=90,hjust=1),legend.position='none') + scale_y_continuous(labels=comma) 

# Mean vs. Variance plots (before and after upper quartile normalization excluding genes filtered out for all plots)
options(scipen=999)
prefiltMeanVar = data.frame(Mean=colMeans(castTrainingSet[,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]),Variance=colVars(castTrainingSet[,which(colnames(castTrainingSet) %in% colnames(readyTrainingSet))]))
#ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
ggplot(prefiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(limits=c(0,30000000000),labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()
postfiltMeanVar = data.frame(Mean=colMeans(readyTrainingSet),Variance=colVars(readyTrainingSet))
ggplot(postfiltMeanVar,aes(x=Mean,y=Variance)) + geom_point() + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + theme(axis.title.x=element_text(vjust=-.5),axis.title.y=element_text(vjust=-.05), plot.margin=unit(c(1,1,1,1),'lines')) + theme_bw()

# Clear out some of the data from memory prior to feature selection.
prefiltMeanVar <- NULL
postfilMeanVar <- NULL
preGeneExamine <- NULL
postGeneExamine <- NULL
vTraining <- NULL
zTraining <- NULL
interTrainingSet <- NULL

# Interact with NCBI's SRA db to get metadata for independent test data set.
require(SRAdb)
require(sqldf)
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
conversion <- conversion[!(conversion$experiment %in% c('ERX040280','ERX040279','ERX040278','ERX040277','ERX040276','ERX036541')),]
# Class is 1 for tumor, 0 for healthy/normal
conversion$class <- 1
# Parse description to mark the normal biological samples.
conversion[grepl('apparently_normal',conversion$description),]$class <- 0
# These results will eventually be joined to the filenames of independent data set to get biological class indicator.

# Load up test data set and prepare for use
setwd("~/rsem-1.2.4/")
listFiles = list.files(pattern="*genes.results")
tmpFileMat <- read.table(listFiles[1],sep='\t',header=T,check.names=F)
tmpFileMat <- tmpFileMat[!grepl('^uc',tmpFileMat$gene_id),]
tmpFileMat$run <- substring(listFiles[1],1,9)
finFileMat <- tmpFileMat
for(i in 2:NROW(listFiles)){
	tmpFileMat <- read.table(listFiles[i],sep='\t',header=T,check.names=F)
	tmpFileMat <- tmpFileMat[!grepl('^uc',tmpFileMat$gene_id),]
	tmpFileMat$run <- substring(listFiles[i],1,9)
	finFileMat <- rbind(finFileMat,tmpFileMat)
}
finFileMat$'transcript_id(s)' <- NULL
finFileMat$length <- NULL
finFileMat$effective_length <- NULL
finFileMat$TPM <- NULL
finFileMat$FPKM <- NULL
finFileMat <- data.frame(run=finFileMat$run,gene_id=finFileMat$gene_id,expected_count=finFileMat$expected_count,check.names=F)
castTestSet <- cast(finFileMat, run ~ gene_id)
castTestSet <- merge(castTestSet, conversion, by = 'run')
# Do some work to remove technical replicates, keeping those with largest library size.
res <- rowSums(castTestSet[,2:20532])
resTotal <- data.frame(run=castTestSet$run,libsize=res)
resMerge <- merge(resTotal, conversion, by = 'run')
keepLargestRun <- unlist(sqldf('select run from resMerge rs join (select sample, max(libsize) as maxLib from resMerge group by sample)a on a.sample = rs.sample and a.maxLib = rs.libsize'))
castTestSet <- castTestSet[castTestSet$run %in% keepLargestRun,]
# Technical replicates now removed, proceed.
testClasses <- castTestSet$class
castTestSet$classLabel <- 'Tumor'
castTestSet[castTestSet$class==0,]$classLabel <- 'Healthy'
testClassLabels <- castTestSet$classLabel
castTestSet$classLabel <- NULL
castTestSet$class <- NULL
castTestSet$study <- NULL
castTestSet$submission <- NULL
castTestSet$sample <- NULL
castTestSet$experiment <- NULL
castTestSet$description <- NULL
castTestSet <- data.frame(class=testClassLabels, castTestSet, check.names=F)
rownames(castTestSet) <- castTestSet$run
castTestSet$run <- NULL
# Need to verify the above removes technical replicates when all data finished processing.

# Write the test set file with class identifier. 
write.table(castTestSet,"~/Thesis/mungedTestSet.txt",sep='\t',col.names=T,row.names=T)
# Another version for professor.
write.table(t(castTestSet),"~/Thesis/preparedTestSet.txt",sep='\t',col.names=T,row.names=T)

# Read back in if needed.
castTestSet <- read.table("~/Thesis/mungedTestSet.txt",sep='\t',header=T,check.names=F,row.names=1)
testClasses <- rep(1, nrow(castTestSet))
testClasses[castTestSet$class=='Healthy'] <- 0
testClassLabels <- castTestSet$class 
castTestSet$class <- NULL
source('~/Thesis/upperQuartileNormalize.R')
# Mean of upper quartiles from training data set
s <- 1874.198
# Scale & transform test data for use.
testUq <- t(upperQuartileScale(t(castTestSet),s))
readyTestSet <- log2(testUq+1.0)

# Write this 'ready' version for professor.
writeTestSet <- data.frame(class=testClassLabels,readyTestSet,check.names=F)
write.table(t(writeTestSet),"~/Thesis/transformedTestSet.txt",sep='\t',col.names=T,row.names=T)

# Check how RefSeq genes in NCBI GEO accession GSE40419 may differ from genes available in training & test data sets.
gse40419 <- read.table("GSE40419_LC-87_RPKM_expression.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
gseGenes <- as.character(unique(gse40419$gene))
trainingSetGenes <- as.character(colnames(castTrainingSet))
trainingSetGenes <- trainingSetGenes[trainingSetGenes != "class"]
trainingSetEntrez <- unique(sub(".*\\|", "", trainingSetGenes))
trainingSetGeneFull <- data.frame(gene = toupper(sub("\\|.*", "", trainingSetGenes)), entrez_id = sub(".*\\|", "", trainingSetGenes), stringsAsFactors = FALSE)
require(annotation)
library(hgu95av2.db)
require(sqldf)

xx <- select(hgu95av2.db, keys = as.character(gse40419$accession), columns = c("SYMBOL", "GENENAME", "ALIAS", "ENTREZID", "REFSEQ"), keytype = "REFSEQ")
xx$ALIAS <- toupper(xx$ALIAS)
xx$SYMBOL <- toupper(xx$SYMBOL)
xx$GENENAME <- toupper(xx$GENENAME)
gseEntrez <- unique(xx[!is.na(xx$ENTREZID), ]$ENTREZID)
joinedToGse <- sqldf("SELECT tf.*, xx.ALIAS, xx.SYMBOL, xx.GENENAME FROM trainingSetGeneFull tf join xx on (tf.entrez_id = xx.ENTREZID)")
# Joining on entrez ID joins all but 107 from training set to GSE matrix
length(unique(joinedToGse$entrez_id))
leftOvers = trainingSetGeneFull[!trainingSetGeneFull$entrez_id %in% unique(joinedToGse$entrez_id), ]
#check the leftovers via entrez ID...
xx <- select(hgu95av2.db, keys = as.character(leftOvers$entrez_id), columns = c("SYMBOL", "GENENAME", "ALIAS", "ENTREZID", "REFSEQ"), keytype = "ENTREZID")
# do these join up to the GSE on REFSEQ?
xx <- xx[!is.na(xx$REFSEQ), ]
gseAddl <- sqldf("select xx.*, gse.gene, gse.accession from xx join gse40419 gse on gse.accession = xx.REFSEQ", drv = "SQLite")
# No, zero.
nrow(gseAddl)
# Looks like based on RefSeq -> Entrez ID conversion then join, these don't match between files.
# At least some seem due to outdated RefSeq accessions in the test matrix file... 

