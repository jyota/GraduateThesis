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

# Plot maximum variance genes class separation prior to transform non informative set
nonInfVars <- cmpVars[cmpVars$setOfVars=='Noninformative Set',]
nonInfVars <- nonInfVars[order(nonInfVars$variance,decreasing=T),]
niGenes <- rownames(nonInfVars)[1:4]
maxVarPlot <- data.frame(gene=c(rep(niGenes[1],length(classes)),rep(niGenes[2],length(classes)),
	rep(niGenes[3],length(classes)),rep(niGenes[4],length(classes))),
    class=c(rep(classes,4)),value=c(castTrainingSet[,which(colnames(castTrainingSet) %in% niGenes[1])],
    	castTrainingSet[,which(colnames(castTrainingSet) %in% niGenes[2])],castTrainingSet[,which(colnames(castTrainingSet) %in% niGenes[3])],
    	castTrainingSet[,which(colnames(castTrainingSet) %in% niGenes[4])]))

maxVarPlot$Class <- 'Tumor'
maxVarPlot[maxVarPlot$class==0,]$Class <- 'Healthy'
ggplot(maxVarPlot,aes(value,fill=Class))+geom_density(alpha=.8)+theme_bw()+scale_x_continuous(labels=comma)+facet_wrap(~ gene)+theme(axis.text.x=element_text(angle=90))+labs(fill='Class')+xlab('')

# Check out post transform plot of these
nonInfVars <- cmpVars[cmpVars$setOfVars=='Noninformative Set',]
nonInfVars <- nonInfVars[order(nonInfVars$variance,decreasing=T),]
niGenes <- rownames(nonInfVars)[1:4]
maxVarPlot <- data.frame(gene=c(rep(niGenes[1],length(classes)),rep(niGenes[2],length(classes)),
	rep(niGenes[3],length(classes)),rep(niGenes[4],length(classes))),
    class=c(rep(classes,4)),value=c(readyTrainingSet[,which(colnames(readyTrainingSet) %in% niGenes[1])],
    	readyTrainingSet[,which(colnames(readyTrainingSet) %in% niGenes[2])],readyTrainingSet[,which(colnames(readyTrainingSet) %in% niGenes[3])],
    	readyTrainingSet[,which(colnames(readyTrainingSet) %in% niGenes[4])]))

maxVarPlot$Class <- 'Tumor'
maxVarPlot[maxVarPlot$class==0,]$Class <- 'Healthy'
ggplot(maxVarPlot,aes(value,fill=Class))+geom_density(alpha=.8)+theme_bw()+scale_x_continuous(labels=comma)+facet_wrap(~ gene)+theme(axis.text.x=element_text(angle=90))+labs(fill='Class')+xlab('')
