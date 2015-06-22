
#af1 <- "~/Desktop/IsolatedSource/mig50_meta0.01_gen115_t1_1_001.dat"
#cf1 <- "~/Desktop/IsolatedSource/mig50_meta0.01_gen115_t2_1_001.dat"
#sf1 <- "~/Desktop/IsolatedSource/mig50_meta0.01_gen501_1_1.dat"


#currentdata <- cf1
#ancientdata <- af1
#sourcedata <- sf1
#ninds= 250
#nsub=250
#gensbtwn="1.00"
#mlne.gens="1"
#maxNe=1000
#nsamps=2
#maxall=256
#IBD=FALSE
#na.s = c("0", "00", "000", "0000", "00000", "000000", "NA")


read.nemo.migr <- function (currentdata, ancientdata, sourcedata, ninds=250, nsub, maxall=256, nsamps=2, gensbtwn="4.00", mlne.gens="4", maxNe=1000, IBD=FALSE, na.s = c("0", "00", "000", "0000", "00000", "000000", "NA")) 
{
    x <- scan(currentdata, n = 4) #read first line of file, gives number of loci, alleles, etc
    #first line tells me:  13 patches, number of loci + 4, number of alleles, digits per genotype
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames <- scan(currentdata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames <- c("Pop", lnames)  #add first column heading name to be "pop"
    
    dat <- scan(currentdata, skip = nloc + 5, what = character(), na.strings = na.s, nlines=ninds*3)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    	#only read in first 3 patches when dealing with this isolated pop stuff
    dat <- data.frame(matrix(dat, ncol = nloc + 5, byrow = TRUE))  #genotype data in matrix format now
    dat <- dat[,-((length(dat[1,])-3):length(dat[1,]))]  #remove last 4 columns of nemo data
    names(dat) <- lnames  #add column names
  

#split genotypes into alleles of 3 digits
    for(i in 2:(nloc+1)){  #go from first column of genotype data to the last column
    	if(i %% 100==0) print(i)	
    	ans <- matrix(NA, ncol=2, nrow=length(dat[,i]))
		for(j in 1:length(dat[,i])){
			test <- strsplit(as.character(dat[,i]), "")[[j]]
			ans[j,1] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[1]
			ans[j,2] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[2]
			}  #end for loop going through individual genotype column, splitting genotypes
		dat[,i] <- ans	
    	}  #end for loop going through genotype columns
#remove loci that don't meet criteria
    remove.loci <- NULL		
    	for(k in 2:(nloc+1)){
    		if(length(table(dat[,k]))==1){ #if an allele at a locus is fixed, remove that locus
    			r.l <- k
    			remove.loci <- c(remove.loci, r.l)
    		}

    	}
    
    
### REPEAT FOR SECOND DATASET
x <- scan(ancientdata, n = 4) #read first line of file, gives number of loci, alleles, etc
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames2 <- scan(ancientdata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames2 <- c("Pop", lnames2)  #add first column heading name to be "pop"
    dat2 <- scan(ancientdata, skip = nloc + 5, what = character(), na.strings = na.s, nlines=ninds*3)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    	#only read in first 3 patches when dealing with this isolated pop stuff
    dat2 <- data.frame(matrix(dat2, ncol = nloc + 5, byrow = TRUE))  #genotype data in matrix format now
    dat2 <- dat2[,-((length(dat2[1,])-3):length(dat2[1,]))]  #remove last 4 columns of nemo data
    names(dat2) <- lnames2  #add column names
    
    
#split genotypes into alleles of 3 digits
    for(i in 2:(nloc+1)){  #go from first column of genotype data to the last column
    	if(i %% 100==0) print(i)	
    	ans <- matrix(NA, ncol=2, nrow=length(dat2[,i]))
		for(j in 1:length(dat2[,i])){
			test <- strsplit(as.character(dat2[,i]), "")[[j]]
			ans[j,1] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[1]
			ans[j,2] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[2]
			}  #end for loop going through individual genotype column, splitting genotypes
		dat2[,i] <- ans
    	}  #end for loop going through genotype columns
    	
#remove loci that don't meet criteria
    remove.loci2 <- NULL		
    	for(k in 2:(nloc+1)){
    		if(length(table(dat2[,k]))==1){ #if an allele at a locus is fixed, remove that locus
    			r.l2 <- k
    			remove.loci2 <- c(remove.loci2, r.l2)
    		}
    	}

 
    
### REPEAT FOR ISOLATED DATASET
x <- scan(sourcedata, n = 4) #read first line of file, gives number of loci, alleles, etc
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames3 <- scan(sourcedata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames3 <- c("Pop", lnames3)  #add first column heading name to be "pop"
    dat3 <- scan(sourcedata, skip = nloc + 5 + ninds*14, what = character(), na.strings = na.s)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    	#only read in patch 15 when dealing with this isolated pop stuff
    dat3 <- data.frame(matrix(dat3, ncol = nloc + 5, byrow = TRUE))  #genotype data in matrix format now
    dat3 <- dat3[,-((length(dat3[1,])-3):length(dat3[1,]))]  #remove last 4 columns of nemo data
    names(dat3) <- lnames3  #add column names
    
    
#split genotypes into alleles of 3 digits
    for(i in 2:(nloc+1)){  #go from first column of genotype data to the last column
    	if(i %% 100==0) print(i)	
    	ans <- matrix(NA, ncol=2, nrow=length(dat3[,i]))
		for(j in 1:length(dat3[,i])){
			test <- strsplit(as.character(dat3[,i]), "")[[j]]
			ans[j,1] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[1]
			ans[j,2] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[2]
			}  #end for loop going through individual genotype column, splitting genotypes
		dat3[,i] <- ans
    	}  #end for loop going through genotype columns
    	
#remove loci that don't meet criteria
    remove.loci3 <- NULL		
    	for(k in 2:(nloc+1)){
    		if(length(table(dat3[,k]))==1){ #if an allele at a locus is fixed, remove that locus
    			r.l3 <- k
    			remove.loci3 <- c(remove.loci3, r.l3)
    		}
    	}
    	
### now take only the loci that meet criteria for BOTH datasets, i.e. where they overlap

	loci.to.remove <- unique(c(remove.loci, remove.loci2, remove.loci3))    
        if(length(loci.to.remove)>0) dat <- dat[,-loci.to.remove]
        if(length(loci.to.remove)>0) dat2 <- dat2[,-loci.to.remove]
        if(length(loci.to.remove)>0) dat3 <- dat3[,-loci.to.remove]
        # note, these numbers indicate the column number, not the locus number, so they will be one higher than the locus number beingremoved

### now randomly take 40 loci out of these

final40 <- sample(lnames[-c(1,loci.to.remove)], size=40)  #take 40 loci out of the remaining loci (get rid of pop column and removed loci column)
if(length(loci.to.remove)==0){
	final40 <- sample(lnames[-1], size=40)
}

finalcurrentdat_presubset <- dat[c("Pop", final40)]
finalancientdat_presubset <- dat2[c("Pop", final40)]
finalsourcedat_presubset <- dat3[c("Pop", final40)]


dat <- NULL
dat2 <- NULL
dat3 <- NULL

#NOW SEPARATE INTO PATCHES
finalcurrentdat_presubset1 <- finalcurrentdat_presubset[finalcurrentdat_presubset$Pop=="1",]
finalancientdat_presubset1 <- finalancientdat_presubset[finalancientdat_presubset$Pop=="1",]

finalcurrentdat_presubset2 <- finalcurrentdat_presubset[finalcurrentdat_presubset$Pop=="2",]
finalancientdat_presubset2 <- finalancientdat_presubset[finalancientdat_presubset$Pop=="2",]

finalcurrentdat_presubset3 <- finalcurrentdat_presubset[finalcurrentdat_presubset$Pop=="3",]
finalancientdat_presubset3 <- finalancientdat_presubset[finalancientdat_presubset$Pop=="3",]


finalsourcedat_presubset15 <- finalsourcedat_presubset[finalsourcedat_presubset$Pop=="15",]



who.to.sample <- sample(1:ninds, size=nsub)


#1
finalcurrentdat1 <- finalcurrentdat_presubset1[who.to.sample,]
finalancientdat1 <- finalancientdat_presubset1[who.to.sample,]
#make more space in memory:
finalcurrentdat_presubset1 <- NULL
finalancientdat_presubset1 <- NULL

#2
finalcurrentdat2 <- finalcurrentdat_presubset2[who.to.sample,]
finalancientdat2 <- finalancientdat_presubset2[who.to.sample,]
#make more space in memory:
finalcurrentdat_presubset2 <- NULL
finalancientdat_presubset2 <- NULL

#3
finalcurrentdat3 <- finalcurrentdat_presubset3[who.to.sample,]
finalancientdat3 <- finalancientdat_presubset3[who.to.sample,]
#make more space in memory:
finalcurrentdat_presubset3 <- NULL
finalancientdat_presubset3 <- NULL


#15
finalsourcedat15 <- finalsourcedat_presubset15[who.to.sample,]
#make more space in memory:
finalsourcedat_presubset15 <- NULL

number.of.loci <- 40
finalnames <-  final40
   
   
   #island model migration
#Nf = patch 1 or 2 or 3 alone
#Nf = 1 or 2 or 3 and source = 15  (misID source)
#Nf = 1 or 2 or 3 and source = 4 or 7 or 11 (correct source)
#Nf = 1 or 2 or 3 and source = 6 or 9 or 13  (distantly related source)
#Nf = 5 and source = 6  (closely related source in the metapop)
#Nf = 1 or 2 or 3 and source = 4-14 (as one source)


	# IBD:
# Nf = 1 & source = 2
# Nf = 6 & source = 7
# Nf = 1 & source = 15
# Nf = 5 & source = 9
# Nf = 1-15 all as one source   
   
   
   #focal pops = 1, 2, 3, 5, 6, all  -> so these become separate input files for every program but MLNe which will take immigrant pops separately
   #source pops = 2, 4, 6, 7, 9, 11, 13, 15, 4-14
   
#### NOW GO INTO EACH INPUT FORMAT and write those data files from this script!

 ## ##  Nf = 1        source = 15  (misID source)	"1.15"

mlne.currdat <- finalcurrentdat1
mlne.ancdat <- finalancientdat1
mlne.srcdat <- finalsourcedat15

#get allele counts for current dataset
	mfulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqscurr <- as.data.frame(table(mlne.currdat[,i]))
		moutputcurr <- cbind(i, names(mlne.currdat[i]), mfreqscurr)
		mfulloutcurr <- rbind(mfulloutcurr, moutputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutcurr) <- c("Number", "Locus", "Allele", "Count")
## REPEAT FOR SECOND TIME POINT
	mfulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqsanc <- as.data.frame(table(mlne.ancdat[,i]))
		moutputanc <- cbind(i, names(mlne.ancdat[i]), mfreqsanc)
		mfulloutanc <- rbind(mfulloutanc, moutputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutanc) <- c("Number", "Locus", "Allele", "Count")	
## REPEAT FOR SOURCE POP
	mfulloutsrc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqssrc <- as.data.frame(table(mlne.srcdat[,i]))
		moutputsrc <- cbind(i, names(mlne.srcdat[i]), mfreqssrc)
		mfulloutsrc <- rbind(mfulloutsrc, moutputsrc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutsrc) <- c("Number", "Locus", "Allele", "Count")
#now need to put zeros where some alleles are not there
mlne.fulllist <- NULL
allele.list <- NULL
allelescurr <- list()
allelesanc <- list()
allelessrc <- list()
for(m in 2:(number.of.loci+1)){
	mtempallelescurr <- subset(mfulloutcurr, mfulloutcurr$Number==m)
	mtempallelesanc <- subset(mfulloutanc, mfulloutanc$Number==m)
	mtempallelessrc <- subset(mfulloutsrc, mfulloutsrc$Number==m)
	mpotalleles <- as.matrix(cbind(seq(0, maxall), rep(0, maxall+1)))
	mfpotallelescurr <- as.matrix(cbind(as.factor(mpotalleles[,1]), mpotalleles[,2]))
	mfpotallelesanc <- mfpotallelescurr
	mfpotallelessrc <- mfpotallelescurr
	#time point 1, put in allele counts
	for(k in 1:length(mfpotallelescurr[,1])){
		for(i in 1:length(mtempallelescurr$Allele)){
			if(mfpotallelescurr[k,1]==as.numeric(as.character(mtempallelescurr$Allele[i]))){mfpotallelescurr[k,2] <- mtempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(mfpotallelesanc[,1])){
		for(i in 1:length(mtempallelesanc$Allele)){
			if(mfpotallelesanc[k,1]==as.numeric(as.character(mtempallelesanc$Allele[i]))){mfpotallelesanc[k,2] <- mtempallelesanc$Count[i]}
	}}
	#time point 2 source pop, put in allele counts
	for(k in 1:length(mfpotallelessrc[,1])){
		for(i in 1:length(mtempallelessrc$Allele)){
			if(mfpotallelessrc[k,1]==as.numeric(as.character(mtempallelessrc$Allele[i]))){mfpotallelessrc[k,2] <- mtempallelessrc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(mfpotallelescurr[,1])){
		if(mfpotallelescurr[zz,2]==0 && mfpotallelesanc[zz,2]==0 && mfpotallelessrc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	if(is.null(bb)==FALSE){
		mfinalalleles1 <- mfpotallelescurr[-bb,]
		mfinalalleles2 <- mfpotallelesanc[-bb,]
		mfinalalleles3 <- mfpotallelessrc[-bb,]
	}else{mfinalalleles1 <- mfpotallelescurr ; mfinalalleles2 <- mfpotallelesanc ; mfinalalleles3 <- mfpotallelessrc}
tempalls1 <- paste(if(length(mfinalalleles1)>2){paste(mfinalalleles1[,2])}else{paste(mfinalalleles1[2])}, collapse=" ")
tempalls2 <- paste(if(length(mfinalalleles2)>2){paste(mfinalalleles2[,2])}else{paste(mfinalalleles2[2])}, collapse=" ")
tempalls3 <- paste(if(length(mfinalalleles3)>2){paste(mfinalalleles3[,2])}else{paste(mfinalalleles3[2])}, collapse=" ")
allelescurr <- c(allelescurr, list(tempalls1))
allelesanc <- c(allelesanc, list(tempalls2))
allelessrc <- c(allelessrc, list(tempalls3))
allcounts <- paste(if(length(mfinalalleles1)>2){paste(length(mfinalalleles1[,2]))}else{paste(length(mfinalalleles1[2]))}, collapse=" ")
allele.list <- c(allele.list, allcounts)
}  #end for loop
#now there are 4 lists, allele.list with each locus's allele count, allelescurr with each element of the list being an allele count for a locus at the current time sample, allelesanc for the ancient sample allele counts per locus, and allelessrc for the ancient source pop sample allele counts per locus
#so loop through and paste out each element of the lists together
for(f in 1:length(allelescurr)){
		mtemplist <- paste(c(allelescurr[f], ""))
		mlne.fulllist <- c(mlne.fulllist, mtemplist)
	}	
mlne.fulllist <- paste(c(unlist(allelescurr),"", unlist(allelesanc),"","", unlist(allelessrc)))
between <- paste(c("0", mlne.gens), collapse=",")

mlne_input <- paste(c(paste("1"), paste("0"), paste(maxNe), paste("2"), paste("0"), paste(number.of.loci), paste(allele.list, collapse=","), paste(nsamps), between, "", mlne.fulllist, "", paste("1"), sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(mlne_input, paste("Outputs/MLNeIN_1.15_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))
#	writeLines(mlne_input, paste("MLNeIN_1.15_", tail(unlist(strsplit(cf1, "/")),n=1), sep=""))
 
 
  ## ##  Nf = 2        source = 15  (misID source)	"1.15"

mlne.currdat <- finalcurrentdat2
mlne.ancdat <- finalancientdat2
mlne.srcdat <- finalsourcedat15

#get allele counts for current dataset
	mfulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqscurr <- as.data.frame(table(mlne.currdat[,i]))
		moutputcurr <- cbind(i, names(mlne.currdat[i]), mfreqscurr)
		mfulloutcurr <- rbind(mfulloutcurr, moutputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutcurr) <- c("Number", "Locus", "Allele", "Count")
## REPEAT FOR SECOND TIME POINT
	mfulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqsanc <- as.data.frame(table(mlne.ancdat[,i]))
		moutputanc <- cbind(i, names(mlne.ancdat[i]), mfreqsanc)
		mfulloutanc <- rbind(mfulloutanc, moutputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutanc) <- c("Number", "Locus", "Allele", "Count")	
## REPEAT FOR SOURCE POP
	mfulloutsrc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqssrc <- as.data.frame(table(mlne.srcdat[,i]))
		moutputsrc <- cbind(i, names(mlne.srcdat[i]), mfreqssrc)
		mfulloutsrc <- rbind(mfulloutsrc, moutputsrc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutsrc) <- c("Number", "Locus", "Allele", "Count")
#now need to put zeros where some alleles are not there
mlne.fulllist <- NULL
allele.list <- NULL
allelescurr <- list()
allelesanc <- list()
allelessrc <- list()
for(m in 2:(number.of.loci+1)){
	mtempallelescurr <- subset(mfulloutcurr, mfulloutcurr$Number==m)
	mtempallelesanc <- subset(mfulloutanc, mfulloutanc$Number==m)
	mtempallelessrc <- subset(mfulloutsrc, mfulloutsrc$Number==m)
	mpotalleles <- as.matrix(cbind(seq("0", maxall), rep(0, maxall+1)))
	mfpotallelescurr <- as.matrix(cbind(as.factor(mpotalleles[,1]), mpotalleles[,2]))
	mfpotallelesanc <- mfpotallelescurr
	mfpotallelessrc <- mfpotallelescurr
	#time point 1, put in allele counts
	for(k in 1:length(mfpotallelescurr[,1])){
		for(i in 1:length(mtempallelescurr$Allele)){
			if(mfpotallelescurr[k,1]==as.numeric(as.character(mtempallelescurr$Allele[i]))){mfpotallelescurr[k,2] <- mtempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(mfpotallelesanc[,1])){
		for(i in 1:length(mtempallelesanc$Allele)){
			if(mfpotallelesanc[k,1]==as.numeric(as.character(mtempallelesanc$Allele[i]))){mfpotallelesanc[k,2] <- mtempallelesanc$Count[i]}
	}}
	#time point 2 source pop, put in allele counts
	for(k in 1:length(mfpotallelessrc[,1])){
		for(i in 1:length(mtempallelessrc$Allele)){
			if(mfpotallelessrc[k,1]==as.numeric(as.character(mtempallelessrc$Allele[i]))){mfpotallelessrc[k,2] <- mtempallelessrc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(mfpotallelescurr[,1])){
		if(mfpotallelescurr[zz,2]==0 && mfpotallelesanc[zz,2]==0 && mfpotallelessrc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	if(is.null(bb)==FALSE){
		mfinalalleles1 <- mfpotallelescurr[-bb,]
		mfinalalleles2 <- mfpotallelesanc[-bb,]
		mfinalalleles3 <- mfpotallelessrc[-bb,]
	}else{mfinalalleles1 <- mfpotallelescurr ; mfinalalleles2 <- mfpotallelesanc ; mfinalalleles3 <- mfpotallelessrc}
tempalls1 <- paste(if(length(mfinalalleles1)>2){paste(mfinalalleles1[,2])}else{paste(mfinalalleles1[2])}, collapse=" ")
tempalls2 <- paste(if(length(mfinalalleles2)>2){paste(mfinalalleles2[,2])}else{paste(mfinalalleles2[2])}, collapse=" ")
tempalls3 <- paste(if(length(mfinalalleles3)>2){paste(mfinalalleles3[,2])}else{paste(mfinalalleles3[2])}, collapse=" ")
allelescurr <- c(allelescurr, list(tempalls1))
allelesanc <- c(allelesanc, list(tempalls2))
allelessrc <- c(allelessrc, list(tempalls3))
allcounts <- paste(if(length(mfinalalleles1)>2){paste(length(mfinalalleles1[,2]))}else{paste(length(mfinalalleles1[2]))}, collapse=" ")
allele.list <- c(allele.list, allcounts)
}  #end for loop
#now there are 4 lists, allele.list with each locus's allele count, allelescurr with each element of the list being an allele count for a locus at the current time sample, allelesanc for the ancient sample allele counts per locus, and allelessrc for the ancient source pop sample allele counts per locus
#so loop through and paste out each element of the lists together
for(f in 1:length(allelescurr)){
		mtemplist <- paste(c(allelescurr[f], ""))
		mlne.fulllist <- c(mlne.fulllist, mtemplist)
	}	
mlne.fulllist <- paste(c(unlist(allelescurr),"", unlist(allelesanc),"","", unlist(allelessrc)))
between <- paste(c("0", mlne.gens), collapse=",")

mlne_input <- paste(c(paste("1"), paste("0"), paste(maxNe), paste("2"), paste("0"), paste(number.of.loci), paste(allele.list, collapse=","), paste(nsamps), between, "", mlne.fulllist, "", paste("1"), sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(mlne_input, paste("Outputs/MLNeIN_2.15_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))
 
 
  ## ##  Nf = 3        source = 15  (misID source)	"1.15"

mlne.currdat <- finalcurrentdat3
mlne.ancdat <- finalancientdat3
mlne.srcdat <- finalsourcedat15

#get allele counts for current dataset
	mfulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqscurr <- as.data.frame(table(mlne.currdat[,i]))
		moutputcurr <- cbind(i, names(mlne.currdat[i]), mfreqscurr)
		mfulloutcurr <- rbind(mfulloutcurr, moutputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutcurr) <- c("Number", "Locus", "Allele", "Count")
## REPEAT FOR SECOND TIME POINT
	mfulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqsanc <- as.data.frame(table(mlne.ancdat[,i]))
		moutputanc <- cbind(i, names(mlne.ancdat[i]), mfreqsanc)
		mfulloutanc <- rbind(mfulloutanc, moutputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutanc) <- c("Number", "Locus", "Allele", "Count")	
## REPEAT FOR SOURCE POP
	mfulloutsrc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqssrc <- as.data.frame(table(mlne.srcdat[,i]))
		moutputsrc <- cbind(i, names(mlne.srcdat[i]), mfreqssrc)
		mfulloutsrc <- rbind(mfulloutsrc, moutputsrc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutsrc) <- c("Number", "Locus", "Allele", "Count")
#now need to put zeros where some alleles are not there
mlne.fulllist <- NULL
allele.list <- NULL
allelescurr <- list()
allelesanc <- list()
allelessrc <- list()
for(m in 2:(number.of.loci+1)){
	mtempallelescurr <- subset(mfulloutcurr, mfulloutcurr$Number==m)
	mtempallelesanc <- subset(mfulloutanc, mfulloutanc$Number==m)
	mtempallelessrc <- subset(mfulloutsrc, mfulloutsrc$Number==m)
	mpotalleles <- as.matrix(cbind(seq("0", maxall), rep(0, maxall+1)))
	mfpotallelescurr <- as.matrix(cbind(as.factor(mpotalleles[,1]), mpotalleles[,2]))
	mfpotallelesanc <- mfpotallelescurr
	mfpotallelessrc <- mfpotallelescurr
	#time point 1, put in allele counts
	for(k in 1:length(mfpotallelescurr[,1])){
		for(i in 1:length(mtempallelescurr$Allele)){
			if(mfpotallelescurr[k,1]==as.numeric(as.character(mtempallelescurr$Allele[i]))){mfpotallelescurr[k,2] <- mtempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(mfpotallelesanc[,1])){
		for(i in 1:length(mtempallelesanc$Allele)){
			if(mfpotallelesanc[k,1]==as.numeric(as.character(mtempallelesanc$Allele[i]))){mfpotallelesanc[k,2] <- mtempallelesanc$Count[i]}
	}}
	#time point 2 source pop, put in allele counts
	for(k in 1:length(mfpotallelessrc[,1])){
		for(i in 1:length(mtempallelessrc$Allele)){
			if(mfpotallelessrc[k,1]==as.numeric(as.character(mtempallelessrc$Allele[i]))){mfpotallelessrc[k,2] <- mtempallelessrc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(mfpotallelescurr[,1])){
		if(mfpotallelescurr[zz,2]==0 && mfpotallelesanc[zz,2]==0 && mfpotallelessrc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	if(is.null(bb)==FALSE){
		mfinalalleles1 <- mfpotallelescurr[-bb,]
		mfinalalleles2 <- mfpotallelesanc[-bb,]
		mfinalalleles3 <- mfpotallelessrc[-bb,]
	}else{mfinalalleles1 <- mfpotallelescurr ; mfinalalleles2 <- mfpotallelesanc ; mfinalalleles3 <- mfpotallelessrc}
tempalls1 <- paste(if(length(mfinalalleles1)>2){paste(mfinalalleles1[,2])}else{paste(mfinalalleles1[2])}, collapse=" ")
tempalls2 <- paste(if(length(mfinalalleles2)>2){paste(mfinalalleles2[,2])}else{paste(mfinalalleles2[2])}, collapse=" ")
tempalls3 <- paste(if(length(mfinalalleles3)>2){paste(mfinalalleles3[,2])}else{paste(mfinalalleles3[2])}, collapse=" ")
allelescurr <- c(allelescurr, list(tempalls1))
allelesanc <- c(allelesanc, list(tempalls2))
allelessrc <- c(allelessrc, list(tempalls3))
allcounts <- paste(if(length(mfinalalleles1)>2){paste(length(mfinalalleles1[,2]))}else{paste(length(mfinalalleles1[2]))}, collapse=" ")
allele.list <- c(allele.list, allcounts)
}  #end for loop
#now there are 4 lists, allele.list with each locus's allele count, allelescurr with each element of the list being an allele count for a locus at the current time sample, allelesanc for the ancient sample allele counts per locus, and allelessrc for the ancient source pop sample allele counts per locus
#so loop through and paste out each element of the lists together
for(f in 1:length(allelescurr)){
		mtemplist <- paste(c(allelescurr[f], ""))
		mlne.fulllist <- c(mlne.fulllist, mtemplist)
	}	
mlne.fulllist <- paste(c(unlist(allelescurr),"", unlist(allelesanc),"","", unlist(allelessrc)))
between <- paste(c("0", mlne.gens), collapse=",")

mlne_input <- paste(c(paste("1"), paste("0"), paste(maxNe), paste("2"), paste("0"), paste(number.of.loci), paste(allele.list, collapse=","), paste(nsamps), between, "", mlne.fulllist, "", paste("1"), sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(mlne_input, paste("Outputs/MLNeIN_3.15_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))
 
 

 
   

       	
   # return(dat)
   
}  #end function
