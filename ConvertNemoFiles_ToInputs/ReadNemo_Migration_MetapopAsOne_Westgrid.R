



read.nemo.migr <- function (currentdata, ancientdata, ninds, nsub, maxall=256, nsamps=2, gensbtwn="4.00", mlne.gens="4", maxNe=1000, IBD=FALSE, na.s = c("0", "00", "000", "0000", "00000", "000000", "NA")) 
{
    x <- scan(currentdata, n = 4) #read first line of file, gives number of loci, alleles, etc
    #first line tells me:  13 patches, number of loci + 4, number of alleles, digits per genotype
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames <- scan(currentdata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames <- c("Pop", lnames)  #add first column heading name to be "pop"
    
    dat <- scan(currentdata, skip = nloc + 755, what = character(), na.strings = na.s, nlines = 2750)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    	# skip patches 1-3  (250 inds * 3 = 750)
    	# max number of lines to read is 2750 because want to skip patch 15 too
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
    dat2 <- scan(ancientdata, skip = nloc + 755, what = character(), na.strings = na.s, nlines = 2750)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
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

### now take only the loci that meet criteria for BOTH datasets, i.e. where they overlap

	loci.to.remove <- unique(c(remove.loci, remove.loci2))    
        if(length(loci.to.remove)>0) dat <- dat[,-loci.to.remove]
        if(length(loci.to.remove)>0) dat2 <- dat2[,-loci.to.remove]
        # note, these numbers indicate the column number, not the locus number, so they will be one higher than the locus number beingremoved

### now randomly take 40 loci out of these

final40 <- sample(lnames[-c(1,loci.to.remove)], size=40)  #take 40 loci out of the remaining loci (get rid of pop column and removed loci column)
if(length(loci.to.remove)==0){
	final40 <- sample(lnames[-1], size=40)
}

finalcurrentdat_presubset <- dat[c("Pop", final40)]
finalancientdat_presubset <- dat2[c("Pop", final40)]


dat <- NULL
dat2 <- NULL


#NOW SEPARATE INTO PATCHES  -- don't need to, have only read in data for patches 4-14 which is in the metapop only


#subset individuas to sample  
inds.4.14 <- 2750  # want to take a sample out of any inds in patches 4-14 file for "all patches" so don't want to use the same who.to.sample which is only out of 250 inds

who.to.sample.more.patches <- sample(1:inds.4.14, size=nsub)




# 4-14 all as one source
finalcurrentdat4.14 <- finalcurrentdat_presubset[who.to.sample.more.patches,]
finalancientdat4.14 <- finalancientdat_presubset[who.to.sample.more.patches,]
#make more space in memory:
finalcurrentdat_presubset <- NULL
finalancientdat_presubset <- NULL


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

####################################################################################################
##   NeEstimator
# format is FSTAT with multiple sample years listed as separate populations


##  ##  PATCH 1  ##  ##
NEcurrdat <- finalcurrentdat4.14	# pop 1
NEancdat <- finalancientdat4.14	# pop 2

#make locus rows for first part of file
	NElocirows <- paste(finalnames, sep="\n")	
#make the first row the # samples, loci, max alleles, digits per allele
	NErow1 <- paste(c(nsamps, number.of.loci, maxall, "3"), collapse="\t")
#now concatenate genotypes back together in each column
for(xx in 2:(number.of.loci+1)){
	if(xx %% 10==0) print(xx)
	#go through each genotype column
	for(zz in 1:length(NEcurrdat[,1])){
		NEcurrdat[zz, xx] <- paste(NEcurrdat[zz, xx], collapse="")
		}
		NEcurrdat[,xx] <- as.matrix(NEcurrdat[,xx][,1])
	}
#now make each row a character
	finalNEcurrdat <- NULL
	for(qq in 1:length(NEcurrdat[,1])){
		aa <- paste(NEcurrdat[qq,], collapse="\t") #make tab delimited for NeEstimator
		finalNEcurrdat <- c(finalNEcurrdat, aa)
	}	
	
# REPEAT FOR ANCIENT DATASET

# change pop name for ancient dataset (second sample)
	NEancdat[,1] <- factor(2)   ##### THIS LINE WILL HAVE TO CHANGE IF DIFFERENT POPS IN ONE FILE
#now concatenate genotypes back together in each column
for(xx in 2:(number.of.loci+1)){
	if(xx %% 10==0) print(xx)
	#go through each genotype column
	for(zz in 1:length(NEancdat[,1])){
		NEancdat[zz, xx] <- paste(NEancdat[zz, xx], collapse="")
		}
		NEancdat[,xx] <- as.matrix(NEancdat[,xx][,1])
	}
#now make each row a character
	finalNEancdat <- NULL
	for(qq in 1:length(NEancdat[,1])){
		aa <- paste(c(nsamps, NEancdat[qq,][-1]), collapse="\t") #make tab delimited for NeEstimator
		finalNEancdat <- c(finalNEancdat, aa)
	}
	
	NeEstimator_input <- c(NErow1, NElocirows, finalNEcurrdat, finalNEancdat)

	writeLines(NeEstimator_input, paste("Outputs/NeEstIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))

####################################################################################################



####################################################################################################
##   ONeSAMP
#format - first row=anything, following rows are locus names plus repeat length, this must be made 2, so alleles must be multiplied by 2, one row with 'pop' folowed by genotype rows as 6 digit alleles
## make an input file only for the current data sample
#will need to be made to windows line endings


##  ##  PATCH 1  ##  ##

OSdat <- finalancientdat4.14

#go through and multiply by two to fix motif repeat number
for(kk in 2:(number.of.loci+1)){
	OSdat[,kk][,1] <- (as.numeric(OSdat[,kk])*2)[1:(length(OSdat[,kk])/2)]
	OSdat[,kk][,2] <- (as.numeric(OSdat[,kk])*2)[((length(OSdat[,kk])/2)+1):length(OSdat[,kk])]

	OSdat[,kk][,1] <- as.character(OSdat[,kk][,1])
	OSdat[,kk][,2] <- as.character(OSdat[,kk][,2])
	
	for(jj in 1:length(OSdat[,kk])){
		if(nchar(OSdat[,kk][jj])<3){OSdat[,kk][jj] <- paste("0", OSdat[,kk][jj], sep="")}
	}
		for(jjj in 1:length(OSdat[,kk])){
		if(nchar(OSdat[,kk][jjj])<3){OSdat[,kk][jjj] <- paste("0", OSdat[,kk][jjj], sep="")}
	}
	for(jj in 1:length(finalcurrentdat4.14[,kk])){
		if(nchar(OSdat[,kk][jj])<2){OSdat[,kk][jj] <- paste("00", OSdat[,kk][jj], sep="")}
        }
}

#now concatenate genotypes back together in each column
#go through each locus
for(xx in 2:(number.of.loci+1)){
	if(xx %% 10==0) print(xx)
	#go through each genotype column
	for(zz in 1:length(OSdat[,1])){
		OSdat[zz, xx] <- paste(OSdat[zz, xx], collapse="")
		}
		OSdat[,xx] <- as.matrix(OSdat[,xx][,1])
	}
	
#make locus rows for first part of file
OSrows <- paste(finalnames, sep="\n")
	#add commas after locus name
	for(jj in 1:length(OSrows)){
		OSrows[jj] <- paste(OSrows[jj], ", 2", sep="")
	}
OSrowsfinal <- paste(c(OSrows,"Pop"), sep="\n")
	# now have the locus rows and the pop row
	
#make the first row the nemo file name
	OSrow1 <- paste(tail(unlist(strsplit(currentdata, "/")),n=1))

#add commas after pop name
	col1 <- as.character(OSdat[,1])
	for(j in 1:length(col1)){
		aaa <- paste(c(paste(col1[j]), ","), collapse="")
		col1[j] <- aaa
	}	
OSdat[,1] <- col1	
		 	
#now make each row a character
	finalOSdat <- NULL
	for(qq in 1:length(OSdat[,1])){
		aa <- paste(OSdat[qq,], collapse=" ") #make space delimited for ONeSAMP
		finalOSdat <- c(finalOSdat, aa)
	}
	
	onesamp_input <- c(OSrow1, OSrowsfinal, finalOSdat)

	writeLines(onesamp_input, paste("Outputs/ONeSampIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))

   
  
 
   
   
    


####################################################################################################
##   Estim
# .gen format, first line can be anything, make it the file name
## WILL NEED TO CONVERT TO WINDOWS LINE ENDINGS


##  ##  PATCH 1  ##  ##

ESTdat <- finalancientdat4.14

#now concatenate genotypes back together in each column
#go through each locus
for(xx in 2:(number.of.loci+1)){
	if(xx %% 10==0) print(xx)
	#go through each genotype column
	for(zz in 1:length(ESTdat[,1])){
		ESTdat[zz, xx] <- paste(ESTdat[zz, xx], collapse="")
		}
		ESTdat[,xx] <- as.matrix(ESTdat[,xx][,1])
	}
	
#make locus rows for first part of file
ESTrows <- paste(finalnames, sep="\n")
ESTrowsfinal <- paste(c(ESTrows,"POP"), sep="\n")
	
#make the first row the nemo file name
	ESTrow1 <- paste(tail(unlist(strsplit(currentdata, "/")),n=1))
	
	#add commas after pop name
	estcol1 <- as.character(ESTdat[,1])
	for(j in 1:length(estcol1)){
		aaa <- paste(c(paste(estcol1[j]), ","), collapse="")
		estcol1[j] <- aaa
	}	
ESTdat[,1] <- estcol1	
		 	
#now make each row a character
	finalESTdat <- NULL
	for(qq in 1:length(ESTdat[,1])){
		aa <- paste(ESTdat[qq,], collapse=" ") #make space delimited for ONeSAMP
		finalESTdat <- c(finalESTdat, aa)
	}
	
	estim_input <- c(ESTrow1, ESTrowsfinal, finalESTdat)

	writeLines(estim_input, paste("Outputs/EstimIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))






####################################################################################################
##   TMVP
# format is allele counts from each time sample


##  ##  PATCH 1  ##  ##

tmvp.currdat <- finalcurrentdat4.14
tmvp.ancdat <- finalancientdat4.14

#get allele counts for current dataset
	fulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		freqscurr <- as.data.frame(table(tmvp.currdat[,i]))
		outputcurr <- cbind(i, names(tmvp.currdat[i]), freqscurr)
		fulloutcurr <- rbind(fulloutcurr, outputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(fulloutcurr) <- c("Number", "Locus", "Allele", "Count")

## REPEAT FOR SECOND TIME POINT
	fulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		freqsanc <- as.data.frame(table(tmvp.ancdat[,i]))
		outputanc <- cbind(i, names(tmvp.ancdat[i]), freqsanc)
		fulloutanc <- rbind(fulloutanc, outputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(fulloutanc) <- c("Number", "Locus", "Allele", "Count")
	
#now need to put zeros where some alleles are not there
fulllist <- NULL
for(m in 2:(number.of.loci+1)){
	tempallelescurr <- subset(fulloutcurr, fulloutcurr$Number==m)
	tempallelesanc <- subset(fulloutanc, fulloutanc$Number==m)
	
	potalleles <- as.matrix(cbind(seq("0", maxall), rep(0, maxall+1)))
	fpotallelescurr <- as.matrix(cbind(as.factor(potalleles[,1]), potalleles[,2]))
	fpotallelesanc <- fpotallelescurr
	
	#time point 1, put in allele counts
	for(k in 1:length(fpotallelescurr[,1])){
		for(i in 1:length(tempallelescurr$Allele)){
			if(fpotallelescurr[k,1]==as.numeric(as.character(tempallelescurr$Allele[i]))){fpotallelescurr[k,2] <- tempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(fpotallelesanc[,1])){
		for(i in 1:length(tempallelesanc$Allele)){
			if(fpotallelesanc[k,1]==as.numeric(as.character(tempallelesanc$Allele[i]))){fpotallelesanc[k,2] <- tempallelesanc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(fpotallelescurr[,1])){
		if(fpotallelescurr[zz,2]==0 && fpotallelesanc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	finalallelescurr <- fpotallelescurr[-bb,]
	finalallelesanc <- fpotallelesanc[-bb,]

	#if statements there for when there is only one allele at a locus -> it doesn't read it as a matrix
	templist <- paste(c(
	paste(c(nsamples=nsamps, 
	if(length(finalallelescurr)>2){paste(length(finalallelescurr[,2]))}else{paste(length(finalallelescurr[2]))}), collapse=" "), 
		#number of samples, number of alleles
	paste(c("0.00", if(length(finalallelescurr)>2){paste(finalallelescurr[,2])}else{paste(finalallelescurr[2])}), collapse=" "), 
		#time point, allele counts
	paste(c(generationsbetweensamples=gensbtwn, if(length(finalallelesanc)>2){paste(finalallelesanc[,2])}else{paste(finalallelesanc[2])}), collapse=" "), 
		#time point, allele counts
	 ""))
	fulllist <- paste(c(fulllist, templist))
}

tmvp_input <- paste(c(paste(number.of.loci), fulllist, sep="\n"))  #add number of loci at start of file

	writeLines(tmvp_input, paste("Outputs/TMVPIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))




####################################################################################################
##   CoNe
#similar format to tmvp except:
### FILE 1 MUST BE OLDEST TIME SAMPLE!!!, i.e. line two is the most recent sample for each locus's data


##  ##  PATCH 1  ##  ##

cone.currdat <- finalcurrentdat4.14
cone.ancdat <- finalancientdat4.14

#get allele counts for current dataset
	conefulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		freqscurr <- as.data.frame(table(cone.currdat[,i]))
		outputcurr <- cbind(i, names(cone.currdat[i]), freqscurr)
		conefulloutcurr <- rbind(conefulloutcurr, outputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(conefulloutcurr) <- c("Number", "Locus", "Allele", "Count")

## REPEAT FOR SECOND TIME POINT
	conefulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		freqsanc <- as.data.frame(table(cone.ancdat[,i]))
		outputanc <- cbind(i, names(cone.ancdat[i]), freqsanc)
		conefulloutanc <- rbind(conefulloutanc, outputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(conefulloutanc) <- c("Number", "Locus", "Allele", "Count")
	
#now need to put zeros where some alleles are not there
cone.fulllist <- NULL
for(m in 2:(number.of.loci+1)){
	tempallelescurr <- subset(conefulloutcurr, conefulloutcurr$Number==m)
	tempallelesanc <- subset(conefulloutanc, conefulloutanc $Number==m)
	
	potalleles <- as.matrix(cbind(seq("0", maxall), rep(0, maxall+1)))
	fpotallelescurr <- as.matrix(cbind(as.factor(potalleles[,1]), potalleles[,2]))
	fpotallelesanc <- fpotallelescurr
	
	#time point 1, put in allele counts
	for(k in 1:length(fpotallelescurr[,1])){
		for(i in 1:length(tempallelescurr$Allele)){
			if(fpotallelescurr[k,1]==as.numeric(as.character(tempallelescurr$Allele[i]))){fpotallelescurr[k,2] <- tempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(fpotallelesanc[,1])){
		for(i in 1:length(tempallelesanc$Allele)){
			if(fpotallelesanc[k,1]==as.numeric(as.character(tempallelesanc$Allele[i]))){fpotallelesanc[k,2] <- tempallelesanc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(fpotallelescurr[,1])){
		if(fpotallelescurr[zz,2]==0 && fpotallelesanc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	conefinalallelescurr <- fpotallelescurr[-bb,]
	conefinalallelesanc <- fpotallelesanc[-bb,]

	#if statements there for when there is only one allele at a locus -> it doesn't read it as a matrix
    ctemplist <- paste(c(
	paste(if(length(conefinalallelesanc)>2){paste(length(conefinalallelesanc[,2]))}else{paste(length(conefinalallelesanc[2]))}, collapse=" "), 
		#number of alleles
	paste(if(length(conefinalallelesanc)>2){paste(conefinalallelesanc[,2])}else{paste(conefinalallelesanc[2])}, collapse=" "),  
		## allele counts most ancient sample
	paste(if(length(conefinalallelescurr)>2){paste(conefinalallelescurr[,2])}else{paste(conefinalallelescurr[2])}, collapse=" "),  
		## allele counts most recent sample 
	 ""))
	cone.fulllist <- paste(c(cone.fulllist, ctemplist))
}	
       	
       	cone_input <- paste(c(paste("0"), paste(nsamples=nsamps), paste(number.of.loci), "", cone.fulllist, sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(cone_input, paste("Outputs/CoNeIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))



################################################################################################
####  COLONY

#complex input file, not inclding any of the parent-offspring-sibship relationship parts





## ## PATCH 1 ## ##

colony.currdat <- finalancientdat4.14

#get allele counts for current dataset
	fulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		freqscurr <- as.data.frame(table(colony.currdat[,i]))
		outputcurr <- cbind(i, names(colony.currdat[i]), freqscurr)
		fulloutcurr <- rbind(fulloutcurr, outputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(fulloutcurr) <- c("Number", "Locus", "Allele", "Count")

tab <- rowSums(table(fulloutcurr$Locus, fulloutcurr$Allele))
tab2 <- sort(tab)
nums.alls.per.loc <- unname(unlist(tab2))

unsortedallelecounts <- unname(unlist(tab))
j <- 0
for(i in fulloutcurr$Number){
	j <- j+1
	fulloutcurr$NumAllsAtLocus[j] <- unsortedallelecounts[i-1]
}  
fulloutcurr$Freq <- fulloutcurr$Count/fulloutcurr$NumAllsAtLocus

#sorted locus names
Locus <- names(tab2)
ord <- 1:40
ordered.loci <- cbind(Locus, as.integer(ord))

for.sorting <- merge(fulloutcurr, ordered.loci, by = "Locus")
for.sorting$order <- as.integer(as.character(for.sorting$V2))
sorted.matrix <- for.sorting[order(for.sorting$order),]
# now the matrix is sorted to match the order of loci from fewest to most alleles per locus

by.locus <- split(sorted.matrix, f= sorted.matrix$order)

all.loci <- NULL
for(k in 1:number.of.loci){
	temp.mat <- by.locus[[k]]
	allele.numbers <- paste(as.character(1:temp.mat$NumAllsAtLocus[1]), collapse=" ")
	allele.frequencies <- paste(temp.mat$Freq, collapse=" ")
	temp.row <- c(allele.numbers, allele.frequencies)
	all.loci <- c(all.loci, temp.row)
}  #this makes the allale frequencies per allele per locus for the earlier part of the input file

#now make the offspring IDs and genotypes part of the input file

#first just match up allele IDs per locus with each possible allele
allele.nums <- NULL
for(zz in 1:length(nums.alls.per.loc)){
	temp.allele.nums <- 1:nums.alls.per.loc[zz]
	allele.nums <- c(allele.nums, temp.allele.nums)
}

sorted.matrix$all.num <- allele.nums
for.matching.genos <- sorted.matrix[,c(1,8,9,3)]  #now each allele per locus is matched to it's numerical ID

#then go through the actual genotype data per individual and paste out the genotypes by numerical identifier

ind.numbers <- 1:nsub  #plain list of numbers
ind.ids <- NULL  # this will be the final list of IDs for each individual
for(jj in 1:nsub){ #go through each number
	individual.id <- paste("O", ind.numbers[jj], sep="")  #past an O in front of it
	ind.ids <- c(ind.ids, individual.id)  #put it in the list
}  #now have a vector of individual identifiers

all.genotypes <- NULL
for(gg in 1:nsub){  #take a row of raw genotype data
	temp.genos <- colony.currdat[gg,] #this is one row, i.e. one individual's raw genotype
	genotype <- ind.ids[gg] #give that individual it's ID from the list
	for(hh in 1:number.of.loci){  #take one locus per that row (per individual)
		aaa <- which(names(temp.genos)==Locus[hh])  
			#find the genotypes that correspond to the locus I want (since the raw data isn't ordered, but I need to keep things ordered ascendingly in the input file)
		to.replace <- temp.genos[[aaa]]  #this is one locus's genotype (diploid)
		temp.to.match <- subset(for.matching.genos[for.matching.genos$Locus==Locus[hh],])  
			#get all possible alleles for that locus
		
		for(ii in 1:length(temp.to.match$order)){  #replace that allele with its identifier number
			if(to.replace[1]==temp.to.match$Allele[ii]){first.geno <- paste(temp.to.match$all.num[ii])}
			if(to.replace[2]==temp.to.match$Allele[ii]){second.geno <- paste(temp.to.match$all.num[ii])}
		}
		geno <- paste(first.geno, second.geno) #this is one locus per one individual's genotype
		genotype <- paste(genotype, geno)  #this is one individual's genotype for all loci
	}
	all.genotypes <- c(all.genotypes, genotype)
}

colony_input <- paste(c(
	paste(tail(unlist(strsplit(currentdata, "/")),n=1)), 
	paste("ColonyOUT_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""), 
	paste(ninds), 
	paste(number.of.loci), 
	paste("1234"), #Seed for random number generator
	paste("0"), # 0/1 Not updating/updating allele frequency
	paste("1"), # 2/1 Dioecious/Monoecious species
	paste("0"), # 0/1 Diploid species/HaploDiploid species
	paste("0 0"), # 0/1 Polygamy/Monogamy for males & females
	paste("0"), # 0/1	Without/with sibship size prior
	paste("1"), # 0/1 Unknown/Known population allele frequency
	paste(nums.alls.per.loc, collapse=" "), #number of alleles per locus from locus 1 in ascending order
	paste(all.loci), #allele numerical identifiers with frequencies
	paste(""),
	paste("1"), # number of runs
	paste("2"),	# Length of run 1 2 3 or 4  short medium long very long
	paste("0"),	# 0/1=Monitor method by Iterate#/Time in second
	paste("100000"),	# Monitor interval in Iterate# / in seconds
	paste("1"),	# Windows version  1 equals DOS mode      0 equals windows
	paste("1"),	# Full likelihood method   0 equals pair likelihood method
	paste("3"),	# 1/2/3=low/medium/high Precision for Full likelihood
	paste(""),
	paste(Locus, collapse=" "), #marker names
	paste(rep(0, length(Locus)), collapse=" "), #0's for codominant
	paste(rep("0.0", length(Locus)), collapse=" "), #dropout rate
	paste(rep("0.0", length(Locus)), collapse=" "), #error rates
	paste(""),
	paste(all.genotypes, sep="/n"),
	paste("0.0 0.0"), 
	paste("0 0"), 
	paste("0"), 
	paste("0"), 
	paste("0"), 
	paste("0"), 
	paste("0"), 
	paste("0"), 
	paste("0"), 
	paste("0"),      sep="\n"))  

	writeLines(colony_input, paste("Outputs/ColonyIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))
	
	
 

##    MLNe

#Nf = patch 1 or 2 or 3  no source					"1-2-3alone"
#Nf = 1-2-3     source = 15  (misID source)			"1-2-3.15"
#Nf = 1-2-3     source = 4-7-11 (correct source)	"1-2-3.4-7-11"
#Nf = 1-2-3     source = 6-9-13						"1-2-3.6-9-13"
#Nf = 5    		source = 6							"5.6"
#Nf = 1-2-3     source = 4-14 (as one source)		"1-2-3.4-14"
#NF = 1			source = 2 (IBD one end)			"ibd.1.2"
#Nf = 6         source = 7							"ibd.6.7"
#Nf = 5         source = 9							"ibd.5.9"
#Nf = 1-15      no source 							"1-15"

####################################################################################################
##   MLNe
#odd input, section from focal current pop, then ancient focal pop, if there's migration going on, set a 3rd ancient pop
#0	  #0= Ne only, 1= Ne and m
#0	 #0=not at drift-migr equil, 1=at equil
#10000  #max Ne allowed
#3   # 0,1,2,3 updates to console (all still goes to output)
#2   # number of threads/ CPUs to use 0 will use all available
#40  # num loci
#7,25,13,5,13,21,9,10,10,8,7,2,7,17,17,13,20,14,8    #num alleles at each locus
#2   #number of samples in time
#0,2  #generations of each sample time going back, 0 for current
# 1 at the end of the file for starting point in max likelihood




## ##  Nf = patch 4-14  no source		"1alone"

mlne.currdat <- finalcurrentdat4.14
mlne.ancdat <- finalancientdat4.14

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
#now need to put zeros where some alleles are not there
mlne.fulllist <- NULL
allele.list <- NULL
allelescurr <- list()
allelesanc <- list()
for(m in 2:(number.of.loci+1)){
	mtempallelescurr <- subset(mfulloutcurr, mfulloutcurr$Number==m)
	mtempallelesanc <- subset(mfulloutanc, mfulloutanc$Number==m)	
	mpotalleles <- as.matrix(cbind(seq("0", maxall), rep(0, maxall+1)))
	mfpotallelescurr <- as.matrix(cbind(as.factor(mpotalleles[,1]), mpotalleles[,2]))
	mfpotallelesanc <- mfpotallelescurr	
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
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(mfpotallelescurr[,1])){
		if(mfpotallelescurr[zz,2]==0 && mfpotallelesanc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	if(is.null(bb)==FALSE){
		mfinalalleles1 <- mfpotallelescurr[-bb,]
		mfinalalleles2 <- mfpotallelesanc[-bb,]
	}else{mfinalalleles1 <- mfpotallelescurr ; mfinalalleles2 <- mfpotallelesanc}
tempalls1 <- paste(if(length(mfinalalleles1)>2){paste(mfinalalleles1[,2])}else{paste(mfinalalleles1[2])}, collapse=" ")
tempalls2 <- paste(if(length(mfinalalleles2)>2){paste(mfinalalleles2[,2])}else{paste(mfinalalleles2[2])}, collapse=" ")
allelescurr <- c(allelescurr, list(tempalls1))
allelesanc <- c(allelesanc, list(tempalls2))
allcounts <- paste(if(length(mfinalalleles1)>2){paste(length(mfinalalleles1[,2]))}else{paste(length(mfinalalleles1[2]))}, collapse=" ")
allele.list <- c(allele.list, allcounts)
}  #end for loop
#now there are 3 lists, allele.list with each locus's allele count, allelescurr with each element of the list being an allele count for a locus at the current time sample abd allelesanc for the ancient sample allele counts per locus
#so loop through and paste out each element of the lists together
for(f in 1:length(allelescurr)){
		mtemplist <- paste(c(allelescurr[f], ""))
		mlne.fulllist <- c(mlne.fulllist, mtemplist)
	}	
mlne.fulllist <- paste(c(unlist(allelescurr),"", unlist(allelesanc)))
between <- paste(c("0", mlne.gens), collapse=",")
mlne_input <- paste(c(paste("0"), paste("0"), paste(maxNe), paste("2"), paste("0"), paste(number.of.loci), paste(allele.list, collapse=","), paste(nsamps), between, "", mlne.fulllist, "", paste("1"), sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(mlne_input, paste("Outputs/MLNeIN_P4-14_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))


 
   

       	
   # return(dat)
   
}  #end function
