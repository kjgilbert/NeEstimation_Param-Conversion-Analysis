# Acquired from Stuart Barker, as used in his 2011 Molecular Ecology paper "Effective population size of natural populations of Drosophila buzzatii, with a comparative evaluation of nine methods of estimation", to analyze outputs from the program TMVP


#TMVP1pout.r

#R function to analyse one "out" file from TMVP1p  checking
#convergence and modes for each variable

#To run in the R console, need to enter the following three lines,
#where Bul3 is the identifier for the population and replicate being 
#analysed (use your own identifier here!!):

#Bul3=read.table('Bul3.out')
#source('TMVP1pout.r')
#mc(Bul3)

#prob=0.1 is set to get 90% HPD limits on univariate modes

# out file columns:
# Column 2 = loglike
# Column 3 = Ne

mc=function(object,cut1=0,cut2=0.5,prob=0.1,alpha=0.5,maxk=100)
{
#convert to mcmc object
	#start and end can be set to trim off ends.
	#thin sets intervals between consecutive observations
	a1=round(cut1*dim(object)[1])+1 
	library(coda)
	start=round(dim(object)[1]/10)+1
	end=dim(object)[1]
	object.mcmc=mcmc(data=object,start=round(dim(object)[1]/10)+1,end=dim(object)[1],thin=1)
	if(a1==1){object.1=object[start:end,]
		#plot individual columns (ie variables)
		#dev.new() ; plot(object.mcmc[,2:3]) # UNCOMMENT THIS IF WANT TO SEE TRACES PLOTTED ALONE
		plot(object.mcmc[,2:3])
	}else{object.1=object[a1:end,]
		object1.mcmc=mcmc(data=object,start=a1,end=dim(object)[1],thin=1)
		#plot individual columns (ie variables)
		# dev.new() ; plot(object1.mcmc[,2:3]) #  dev.new() ; par(mfrow=c(3,1)) ; plot(fit2) ; plot(fit3) # UNCOMMENT THIS IF WANT TO SEE TRACES PLOTTED ALONE
		plot(object1.mcmc[,2:3])
	}
#plot individual columns (ie variables)
	cat("Gelman and Rubin's convergence diagnostic","\n")
	#To run this you need two parallel chains
	#I have just broken the original columns 2 - 3 in two.  
	a=round(cut1*dim(object.mcmc)[1])+1 
	b=round(cut2*dim(object.mcmc)[1])
	object.list=mcmc.list(mcmc(object.mcmc[a:b,2]),mcmc(object.mcmc[(b+1):dim(object.mcmc)[1],2]))
	cat("Variable:log(loglike)","\n")
	print(gelman.diag(object.list, confidence = 0.95, transform=FALSE))
	object.list=mcmc.list(mcmc(object.mcmc[a:b,3]),mcmc(object.mcmc[(b+1):dim(object.mcmc)[1],3]))
	cat("Variable:mean_Ne","\n")
	print(gelman.diag(object.list, confidence = 0.95, transform=FALSE))
	library(locfit)
	fit2=locfit(~object.1[,2],alpha=alpha)
	fit3=locfit(~object.1[,3],alpha=alpha)
	cat("*********For the following, prob is set at = ", prob,"************", "\n") ; cat("\n")
	stats=loc1stats(object.1[,2],prob=prob)
	mode.x=stats[1]
	cat("\n") ; cat("******** mode for log_like = ", mode.x, "********","\n")
	cat("\n") ; cat("** Limits will occur in pairs if there is more than one mode.**","\n")
	cat("** Choose the appropriate set that encloses the reported mode.********","\n")
	print(stats[-1])
	cat("\n") ; cat("______________________________________________________________________","\n")
	cat("\n") ; cat("*********For the following, prob is set at = ", prob,"************", "\n")
	cat("\n")
	stats=loc1stats(object.1[,3],prob=prob)
	mode.x=stats[1]
	cat("\n") ; cat("******** mode for mean Ne = ", mode.x, "********","\n")
	cat("\n") ; cat("** Limits will occur in pairs if there is more than one mode.**","\n")
	cat("** Choose the appropriate set that encloses the reported mode.********","\n")
	print(stats[-1])
	cat("\n") ; cat("______________________________________________________________________","\n")
	#  dev.new() ; par(mfrow=c(3,1)) ; plot(fit2) ; plot(fit3) # UNCOMMENT THIS IF WANT TO SEE TRACES PLOTTED ALONE
	par(mfrow=c(3,1)) ; plot(fit2) ; plot(fit3)
	answers <- c(mode.x,stats[-1]) # these are the values I want to put in my file, mode for mean Ne, lower CI, upper CI
	return(answers)
}

loc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

loc2plotw <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,wt,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100,weight=wt)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100,weight=wt)
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}
gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
#	d1 <- (x-px)^2+(y-py)^2
#	best <- d1 == min(d1)
#	lev <- mean(fitted(fit)[best])
	lev <- predict.locfit(fit,list(px,py))
	slev <- sort(fitted(fit))
	indic <- slev <= lev
	sum(indic)/length(x)
}
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}
loc1stats <- function(x,prob,alpha=0.5,xlim,...)
{
	if(missing(xlim)) fit <- locfit(~x,alpha=alpha)
	else fit <- locfit(~x,alpha=alpha,xlim=xlim)
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}	
	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
	print("log difference from max is ")
	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}			
linestats <- function(x,y,prob)
{
	n <- length(x)
	minval <- x[1]
	maxval <- x[n]
	# This bit just guarantees that irrespective of what is 
	# in x and y, we have 10000 *evenly* spaced (interpolated) points between
	# max and min
	xx <- seq(minval,maxval,len=10000)
	yy <- numeric(10000)
	sum1 <- 0
	for(j in 1:10000){
		yy[j] <- approxfun(x,y)(xx[j])
		sum1 <- sum1 + yy[j]
	}
	x.modef <- max(yy)
	x.mode <- xx[yy == x.modef]
	if(length(x.mode)>1)x.mode <- x.mode[1]
	
	yy2 <- sort(yy)
	pval <- 0
	for(j in 1:10000){
		pval <- pval+yy2[j]/sum1
		if(pval > prob)break
	}
	lev <- yy2[j]
	print("log difference from max is ")
	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	ii <- 2
	flip <- TRUE
	for(j in 2:length(xx)){
		if(flip && yy[j] > lev){
			l1[[ii]] <- xx[j-1]
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && yy[j] < lev){
			l1[[ii]] <- xx[j]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip &&  j == length(xx)){
			l1[[ii]] <- xx[j]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}
tloc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha)
	else fit <- locfit(~x+y,alpha=alpha,
	xlim=xlim)
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

