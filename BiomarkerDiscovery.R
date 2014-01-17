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

