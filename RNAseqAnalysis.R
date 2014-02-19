# Check out how violation of normal-based parametric assumptions may have influenced results.

# Check variances. Do our 'informative set of genes' have a tendency for higher variance?
x <- colVars(readyTrainingSet[,which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))])
y <- colVars(readyTrainingSet[,-which(colnames(readyTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))])
cmpVars <- data.frame(variance=c(x,y),setOfVars=c(rep("Informative Set",length(x)),rep("Noninformative Set",length(y))))
# Plot indicates for this data that informative set of genes selected with parametric methods found genes with somewhat higher variance important.
ggplot(cmpVars,aes(y=variance,x=setOfVars,fill=setOfVars))+stat_boxplot(geom='errorbar')+geom_boxplot()+theme_bw()+guides(fill=FALSE)+ylab("Variance")+xlab("")+coord_cartesian(ylim=c(0,10))

# Let's check out the distributions of these genes' variances prior to the normalization transform.
x <- colVars(interTrainingSet[,which(colnames(interTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))])
y <- colVars(interTrainingSet[,-which(colnames(interTrainingSet) %in% unlist(informativeSet[informativeSet$T2>4.0,3:7]))])
cmpVars <- data.frame(variance=c(x,y),setOfVars=c(rep("Informative Set",length(x)),rep("Noninformative Set",length(y))))
# Now the 'non-informative' gene set box plot shifted up. This indicates the transformation variance may have influenced the parametric based gene selection.
ggplot(cmpVars,aes(y=variance,x=setOfVars,fill=setOfVars))+stat_boxplot(geom='errorbar')+geom_boxplot()+theme_bw()+guides(fill=FALSE)+ylab("Variance")+xlab("")+coord_cartesian(ylim=c(0,80000000))+scale_y_continuous(labels=comma)
