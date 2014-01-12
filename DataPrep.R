require(data.table)
require(ggplot2)
require(reshape)
require(limma)
setwd("~/data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/")
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
write.table(castTrainingSet,"mungedTrainingSet.txt",sep='\t',col.names=T,row.names=T)
