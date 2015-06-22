setwd("  ")

m1 <- t(matrix(scan("priors.txt"),nrow=9))
lambda <- (-0.2)
ne <- m1[,1]
ne1 <- (ne^lambda - 1)/lambda
m <- m1[,2]
ld <- m1[,3]
lnb <- m1[,4]
hetx <- m1[,5] 
xhet <- m1[,6] 
nals <- m1[,7]
mho <- m1[,8]
vho <- m1[,9]
target.s <- scan("output.txt")
source("make_pd.r")
result1 <- makepd8(target.s[c(2,3,4,5,6,7,8,9)],ne1,cbind(m,ld,lnb,hetx,xhet,nals,mho,vho),0.01,seq(1,len=50000),F)
result1$x <- (lambda*result1$x+1)^(1/lambda)
source("loc2plot.r")
library(locfit)
#mode <- loc1stats(result1$x,prob=0.05)[1]
mean <- (lambda*result1$predmean+1)^(1/lambda)
median <- median(result1$x)
vari <- var(result1$x)
qntlci <- quantile(result1$x,c(0.025,0.975))
#hpdlu <- loc1stats(result1$x,prob=0.05)[2:3]
write(c("Below are the mean, median, and 95% credible limits for the posterior distribution of Ne from OneSamp"),file=paste("final.txt"),ncol=1,append=T)
write(c("mean  median    lower95%CL   upper95%CL"),file=paste("final.txt"),ncol=1,append=T)
write(c(mean,   median,  qntlci[1],   qntlci[2]),file=paste("final.txt"),ncol=4,append=T)
#q()
