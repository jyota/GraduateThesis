fischerLda<-function(X, Y){	
require(psych)
require(geigen)
# X is matrix of variables that may vary per class (for example, gene expression levels)
# Y is vector or matrix of class variable(s) (for example, a vector of binary classes specifying cancerous vs. non-cancerous)
# Uses framework of matrix regression for LDA related analyses.
# Here, a different approach is attempted than the MASS lda function.
# This solves the eigenproblem for Hv = wEv, where H is the hypothesis matrix and E is the error matrix.
# Then, projects the data onto the eigenvector associated with the largest eigenvalue for the two class problem.
# This does not rescale the data like MASS lda.
Y_ = matrix(ncol=NCOL(Y)+1,nrow=NROW(Y))
Y_[,2:(NCOL(Y)+1)]=Y
Y_[,1]=1
if(rcond(crossprod(Y_)) < .Machine$double.eps){
cat('Matrix singular, cancelling.\n')
return()
}else{
BETA = solve(crossprod(Y_)) %*% t(Y_) %*% X
}

X_ = Y_ %*% BETA
ERROR_ = X - X_

MEAN_X=matrix(ncol=NCOL(X),nrow=NROW(X))
if (NCOL(X)==1){
MEAN_X[,1]=mean(X)
}else {
for(k in 1:NCOL(X)){
MEAN_X[,k]=mean(X[,k])
}
}

SSCP_regression = crossprod(X_) - crossprod(MEAN_X)
SSCP_residual = crossprod(ERROR_)
SSCP_total = crossprod(X) - crossprod(MEAN_X)

if(rcond(SSCP_residual) < .Machine$double.eps){
cat('SSCP residual matrix singular, cancelling.\n')
return()
}else{
T2 = tr(SSCP_regression %*% solve(SSCP_residual))
# Assume NOT a symmetric matrix.
    eigs <- geigen(A=SSCP_regression,B=SSCP_residual,symmetric=F)
    # We're only doing two class classification here so all we care about to get the scores from
    # X is the eigenvector associated with the largest eigenvalue.
    # Not re-scaling scores around zero like MASS lda
    #return(eigs)
    scores <- rowSums(X %*% diag(eigs$vectors[,1]))
    subDat <- data.frame(X,class=Y,check.names=F)
    scoreClass <- data.frame(score=scores,class=Y,class0Prob=rep(NA,length(scores)),class1Prob=rep(NA,length(scores)),predictedClass=rep(NA,length(scores)))
    N = nrow(X)
    J = length(unique(Y))
    t = 1 # Since this is a two class problem, t is always going to be 1
    nj0 = length(Y[Y==0])
    nj1 = length(Y[Y==1])
    wbarj0 = mean(scoreClass[scoreClass$class==0,]$score)
    wbarj1 = mean(scoreClass[scoreClass$class==1,]$score)
    print(c(wbarj0,wbarj1))
    qj0 = length(Y[Y==0])/length(Y)
    qj1 = length(Y[Y==1])/length(Y)
    for(i in 1:NROW(scoreClass)){
     scoreClass[i,]$class0Prob = (scores[i] - wbarj0)^2 - 0.72*log(qj0)
     scoreClass[i,]$class1Prob = (scores[i] - wbarj1)^2 - 0.72*log(qj1)
if(scoreClass[i,]$class0Prob < scoreClass[i,]$class1Prob){
scoreClass[i,]$predictedClass = 0
}else{
scoreClass[i,]$predictedClass = 1
}
}
return(scoreClass)
}
}