
makepd <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)


	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		x3 <- scaled.sumstat[,3][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2+x3,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2],x3=target.s[3]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2],
	#	scaled.sumstat[,3][wt1]=target.s[3]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1]),predmean = predmean, fv = fv)
	}
	l1
}


makepd2 <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)


	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1]),
		predmean = predmean, fv = fv)
	}
	l1
}


makepd5 <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	scaled.sumstat[,4] <- normalise(sumstat[,4],sumstat[,4][gwt])
	scaled.sumstat[,5] <- normalise(sumstat[,5],sumstat[,5][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	target.s[4] <- normalise(target[4],sumstat[,4][gwt])
	target.s[5] <- normalise(target[5],sumstat[,5][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2 +
	(scaled.sumstat[,4]-target.s[4])^2 +
	(scaled.sumstat[,5]-target.s[5])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)

	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		x3 <- scaled.sumstat[,3][wt1]
		x4 <- scaled.sumstat[,4][wt1]
		x5 <- scaled.sumstat[,5][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2+x3+x4+x5,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2],x3=target.s[3],x4=target.s[4],
			x5=target.s[5]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2],
	#	scaled.sumstat[,3][wt1]=target.s[3],
	#	scaled.sumstat[,4][wt1]=target.s[4],
	#	scaled.sumstat[,5][wt1]=target.s[5]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1],sumstat[,4][wt1],sumstat[,5][wt1]),predmean = predmean, fv = fv)
	}
	l1
}

makepd4 <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	scaled.sumstat[,4] <- normalise(sumstat[,4],sumstat[,4][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	target.s[4] <- normalise(target[4],sumstat[,4][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2 +
	(scaled.sumstat[,4]-target.s[4])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)

	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		x3 <- scaled.sumstat[,3][wt1]
		x4 <- scaled.sumstat[,4][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2+x3+x4,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2],x3=target.s[3],x4=target.s[4]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2],
	#	scaled.sumstat[,3][wt1]=target.s[3],
	#	scaled.sumstat[,4][wt1]=target.s[4]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1],sumstat[,4][wt1]),predmean = predmean, fv = fv)
	}
	l1
}

makepd6 <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	scaled.sumstat[,4] <- normalise(sumstat[,4],sumstat[,4][gwt])
	scaled.sumstat[,5] <- normalise(sumstat[,5],sumstat[,5][gwt])
	scaled.sumstat[,6] <- normalise(sumstat[,6],sumstat[,6][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	target.s[4] <- normalise(target[4],sumstat[,4][gwt])
	target.s[5] <- normalise(target[5],sumstat[,5][gwt])
	target.s[6] <- normalise(target[6],sumstat[,6][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2 +
	(scaled.sumstat[,4]-target.s[4])^2 +
	(scaled.sumstat[,5]-target.s[5])^2 +
	(scaled.sumstat[,6]-target.s[6])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)

	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		x3 <- scaled.sumstat[,3][wt1]
		x4 <- scaled.sumstat[,4][wt1]
		x5 <- scaled.sumstat[,5][wt1]
		x6 <- scaled.sumstat[,6][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2+x3+x4+x5+x6,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2],x3=target.s[3],x4=target.s[4],x5=target.s[5],x6=target.s[6]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2],
	#	scaled.sumstat[,3][wt1]=target.s[3],
	#	scaled.sumstat[,4][wt1]=target.s[4],
	#	scaled.sumstat[,5][wt1]=target.s[5],
	#	scaled.sumstat[,6][wt1]=target.s[6]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1],sumstat[,4][wt1],sumstat[,5][wt1],sumstat[,6][wt1]),predmean = predmean, fv = fv)
	}
	l1
}


makepd8 <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{

	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	scaled.sumstat[,4] <- normalise(sumstat[,4],sumstat[,4][gwt])
	scaled.sumstat[,5] <- normalise(sumstat[,5],sumstat[,5][gwt])
	scaled.sumstat[,6] <- normalise(sumstat[,6],sumstat[,6][gwt])
	scaled.sumstat[,7] <- normalise(sumstat[,7],sumstat[,7][gwt])
	scaled.sumstat[,8] <- normalise(sumstat[,8],sumstat[,8][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	target.s[4] <- normalise(target[4],sumstat[,4][gwt])
	target.s[5] <- normalise(target[5],sumstat[,5][gwt])
	target.s[6] <- normalise(target[6],sumstat[,6][gwt])
	target.s[7] <- normalise(target[7],sumstat[,7][gwt])
	target.s[8] <- normalise(target[8],sumstat[,8][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2 +
	(scaled.sumstat[,4]-target.s[4])^2 +
	(scaled.sumstat[,5]-target.s[5])^2 +
	(scaled.sumstat[,6]-target.s[6])^2 +
	(scaled.sumstat[,7]-target.s[7])^2 +
	(scaled.sumstat[,8]-target.s[8])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)

	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	
	if(rejmethod){
		l1 <- list(x=x[wt1],wt=0)
	}
	else{
		regwt <- 1-dist[wt1]^2/abstol^2
		x1 <- scaled.sumstat[,1][wt1]
		x2 <- scaled.sumstat[,2][wt1]
		x3 <- scaled.sumstat[,3][wt1]
		x4 <- scaled.sumstat[,4][wt1]
		x5 <- scaled.sumstat[,5][wt1]
		x6 <- scaled.sumstat[,6][wt1]
		x7 <- scaled.sumstat[,7][wt1]
		x8 <- scaled.sumstat[,8][wt1]
		fit1 <- lm(x[wt1] ~ x1+x2+x3+x4+x5+x6+x7+x8,weight=regwt)
		predmean <- predict.lm(fit1,data.frame(
			x1=target.s[1],x2=target.s[2],x3=target.s[3],x4=target.s[4],x5=target.s[5],x6=target.s[6],x7=target.s[7],x8=target.s[8]))
	#	predmean <- predict.lm(fit1,data.frame(scaled.sumstat[,1][wt1]=target.s[1],
	#	scaled.sumstat[,2][wt1]=target.s[2],
	#	scaled.sumstat[,3][wt1]=target.s[3],
	#	scaled.sumstat[,4][wt1]=target.s[4],
	#	scaled.sumstat[,5][wt1]=target.s[5],
	#	scaled.sumstat[,6][wt1]=target.s[6]
	#	scaled.sumstat[,7][wt1]=target.s[7]
	#	scaled.sumstat[,8][wt1]=target.s[8]))
		fv <- predict.lm(fit1)
	
		l1 <- list(x=x[wt1]+predmean-fv,vals = x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1],sumstat[,4][wt1],sumstat[,5][wt1],sumstat[,6][wt1],sumstat[,7][wt1],sumstat[,8][wt1]),predmean = predmean, fv = fv)
	}
	l1
}



weightfun <- function(target,x,sumstat,tol,gwt)
{
	scaled.sumstat <- sumstat
	scaled.sumstat[,1] <- normalise(sumstat[,1],sumstat[,1][gwt])
	scaled.sumstat[,2] <- normalise(sumstat[,2],sumstat[,2][gwt])
	scaled.sumstat[,3] <- normalise(sumstat[,3],sumstat[,3][gwt])
	target.s <- target
	target.s[1] <- normalise(target[1],sumstat[,1][gwt])
	target.s[2] <- normalise(target[2],sumstat[,2][gwt])
	target.s[3] <- normalise(target[3],sumstat[,3][gwt])
	dist <- sqrt((scaled.sumstat[,1]-target.s[1])^2 + 
	(scaled.sumstat[,2]-target.s[2])^2 + 
	(scaled.sumstat[,3]-target.s[3])^2)
	dist[!gwt] <- floor(max(dist[gwt])+10)


	abstol <- quantile(dist,tol)
	wt1 <- dist < abstol
	regwt <- 1-dist[wt1]^2/abstol^2
	list(x=x[wt1],wt=regwt,ss = cbind(sumstat[,1][wt1],sumstat[,2][wt1],
	sumstat[,3][wt1]))
}

normalise <- function(x,y){
(x-(mean(y)))/sqrt(var(y))
}
