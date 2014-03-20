pooledCov <- function(x,grouping){
	dat <- data.frame(x, class=grouping)
	g0  <- var(dat[dat$class==0,1:NCOL(x)])
	g1  <- var(dat[dat$class==1,1:NCOL(x)])
	Sp  <- (((length(grouping[grouping==0])-1) * g0) + ((length(grouping[grouping==1])-1) * g1))/(length(grouping)-2)
	return(Sp)
}

linscore <- function(x, grouping){
	dat <- data.frame(x, class=grouping)
	m0 <- colMeans(dat[dat$class==0,1:NCOL(x)])
	m1 <- colMeans(dat[dat$class==1,1:NCOL(x)])
	p0 <- length(grouping[grouping==0])/length(grouping)
	p1 <- length(grouping[grouping==1])/length(grouping)	
	Sp <- pooledCov(x,grouping)
	#s0 <- -0.5 * rep((crossprod(m0,solve(Sp)) %*% m0),NROW(x)) + (crossprod(m0,solve(Sp)) %*% t(as.matrix(dat[,1:NCOL(x)]))) + log(p0)
	#s1 <- -0.5 * rep((crossprod(m1,solve(Sp)) %*% m1),NROW(x)) + (crossprod(m1,solve(Sp)) %*% t(as.matrix(dat[,1:NCOL(x)]))) + log(p1)
	s0 <- (as.matrix(dat[,1:NCOL(x)]) %*% solve(Sp) %*% as.matrix(m0)) - rep(0.5 * (m0 %*% solve(Sp) %*% as.matrix(m0)),nrow(x)) + log(p0)
	s1 <- (as.matrix(dat[,1:NCOL(x)]) %*% solve(Sp) %*% as.matrix(m1)) - rep(0.5 * (m1 %*% solve(Sp) %*% as.matrix(m1)),nrow(x)) + log(p1)
	result <- matrix(ncol=2,nrow=nrow(x))
	result[,1] <- s0
	result[,2] <- s1
	colnames(result) <- c('class0','class1')
	return(data.frame(result))
}

