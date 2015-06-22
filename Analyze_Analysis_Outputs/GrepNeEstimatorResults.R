#read file in

setwd("  ")

all.files <- list("SmallerSample_NeEstOUT_P5_mig50_meta0.01_gen115_t2_1_001.dat","SmallerSample_NeEstOUT_P5_mig50_meta0.01_gen115_t2_1_002.dat"
)


results <- NULL
for(i in all.files){
	d <- readLines(i)
	str1 <- "LINKAGE DISEQUILIBRIUM METHOD"
	str2 <- "HETEROZYGOTE EXCESS METHOD"
	str3 <- "MOLECULAR COANCESTRY METHOD"
	str4 <- "(Pollak)"
	str5 <- "(Nei/Tajima)"
	str6 <- "(Jorde/Ryman)"

	LD.line1 <- grep(str1, d)[1]
	LD.line2 <- grep(str1, d)[2]
	het.line1 <- grep(str2, d)[1]
	het.line2 <- grep(str2, d)[2]
	coanc.line1 <- grep(str3, d)[1]
	coanc.line2 <- grep(str3, d)[2]
	Pollak.line <- grep(str4, d)
	NeiTaj.line <- grep(str5, d)
	JorRy.line <- grep(str6, d)
	

## LD method
	temp.LDNe1 <- d[LD.line1 + 6]
	temp.LDNe1 <- gsub(" +", "_", temp.LDNe1)
	LDNe1 <- unlist(strsplit(temp.LDNe1, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.LDNe1.lci <- d[LD.line1 + 12]
	temp.LDNe1.lci <- gsub(" +", "_", temp.LDNe1.lci)
	LDNe1.lci <- unlist(strsplit(temp.LDNe1.lci, "_"))[5]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.LDNe1.uci <- d[LD.line1 + 13]
	temp.LDNe1.uci <- gsub(" +", "_", temp.LDNe1.uci)
	LDNe1.uci <- unlist(strsplit(temp.LDNe1.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
	temp.LDNe2 <- d[LD.line2 + 6]
	temp.LDNe2 <- gsub(" +", "_", temp.LDNe2)
	LDNe2 <- unlist(strsplit(temp.LDNe2, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.LDNe2.lci <- d[LD.line2 + 12]
	temp.LDNe2.lci <- gsub(" +", "_", temp.LDNe2.lci)
	LDNe2.lci <- unlist(strsplit(temp.LDNe2.lci, "_"))[5]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.LDNe2.uci <- d[LD.line2 + 13]
	temp.LDNe2.uci <- gsub(" +", "_", temp.LDNe2.uci)
	LDNe2.uci <- unlist(strsplit(temp.LDNe2.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
## het excess method
	temp.hetNe1 <- d[het.line1 + 5]
	temp.hetNe1 <- gsub(" +", "_", temp.hetNe1)
	hetNe1 <- unlist(strsplit(temp.hetNe1, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.hetNe1.lci <- d[het.line1 + 8]
	temp.hetNe1.lci <- gsub(" +", "_", temp.hetNe1.lci)
	hetNe1.lci <- unlist(strsplit(temp.hetNe1.lci, "_"))[3]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.hetNe1.uci <- d[het.line1 + 9]
	temp.hetNe1.uci <- gsub(" +", "_", temp.hetNe1.uci)
	hetNe1.uci <- unlist(strsplit(temp.hetNe1.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
	temp.hetNe2 <- d[het.line2 + 5]
	temp.hetNe2 <- gsub(" +", "_", temp.hetNe2)
	hetNe2 <- unlist(strsplit(temp.hetNe2, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.hetNe2.lci <- d[het.line2 + 8]
	temp.hetNe2.lci <- gsub(" +", "_", temp.hetNe2.lci)
	hetNe2.lci <- unlist(strsplit(temp.hetNe2.lci, "_"))[3]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.hetNe2.uci <- d[het.line2 + 9]
	temp.hetNe2.uci <- gsub(" +", "_", temp.hetNe2.uci)
	hetNe2.uci <- unlist(strsplit(temp.hetNe2.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
## mol coancestry method
	temp.coancNe1 <- d[coanc.line1 + 4]
	temp.coancNe1 <- gsub(" +", "_", temp.coancNe1)
	coancNe1 <- unlist(strsplit(temp.coancNe1, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.coancNe1.lci <- d[coanc.line1 + 7]
	temp.coancNe1.lci <- gsub(" +", "_", temp.coancNe1.lci)
	coancNe1.lci <- unlist(strsplit(temp.coancNe1.lci, "_"))[5]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.coancNe1.uci <- d[coanc.line1 + 8]
	temp.coancNe1.uci <- gsub(" +", "_", temp.coancNe1.uci)
	coancNe1.uci <- unlist(strsplit(temp.coancNe1.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
	temp.coancNe2 <- d[coanc.line2 + 4]
	temp.coancNe2 <- gsub(" +", "_", temp.coancNe2)
	coancNe2 <- unlist(strsplit(temp.coancNe2, "_"))[4] ## THIS IS AT MAF 0.05 CUTOFF
	temp.coancNe2.lci <- d[coanc.line2 + 7]
	temp.coancNe2.lci <- gsub(" +", "_", temp.coancNe2.lci)
	coancNe2.lci <- unlist(strsplit(temp.coancNe2.lci, "_"))[5]  ## THIS IS AT MAF 0.05 CUTOFF
	temp.coancNe2.uci <- d[coanc.line2 + 8]
	temp.coancNe2.uci <- gsub(" +", "_", temp.coancNe2.uci)
	coancNe2.uci <- unlist(strsplit(temp.coancNe2.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
## temporal methods
	# Pollak method
	temp.PNe <- d[Pollak.line + 4]
	temp.PNe <- gsub(" +", "_", temp.PNe)
	PNe <- unlist(strsplit(temp.PNe, "_"))[5] ## THIS IS AT MAF 0.05 CUTOFF
	temp.PNe.lci <- d[Pollak.line + 9]
	temp.PNe.lci <- gsub(" +", "_", temp.PNe.lci)
	PNe.lci <- unlist(strsplit(temp.PNe.lci, "_"))[6]  ## THIS IS AT MAF 0.05 CUTOFF, jackknife
	temp.PNe.uci <- d[Pollak.line + 10]
	temp.PNe.uci <- gsub(" +", "_", temp.PNe.uci)
	PNe.uci <- unlist(strsplit(temp.PNe.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF
	
	# Nei/Tajima method
	temp.NTNe <- d[NeiTaj.line[2] + 4]
	temp.NTNe <- gsub(" +", "_", temp.NTNe)
	NTNe <- unlist(strsplit(temp.NTNe, "_"))[5] ## THIS IS AT MAF 0.05 CUTOFF
	temp.NTNe.lci <- d[NeiTaj.line[2] + 9]
	temp.NTNe.lci <- gsub(" +", "_", temp.NTNe.lci)
	NTNe.lci <- unlist(strsplit(temp.NTNe.lci, "_"))[6]  ## THIS IS AT MAF 0.05 CUTOFF, jackknife
	temp.NTNe.uci <- d[NeiTaj.line[2] + 10]
	temp.NTNe.uci <- gsub(" +", "_", temp.NTNe.uci)
	NTNe.uci <- unlist(strsplit(temp.NTNe.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF

	# Jorde/Ryman method
	temp.JRNe <- d[JorRy.line[2] + 4]
	temp.JRNe <- gsub(" +", "_", temp.JRNe)
	JRNe <- unlist(strsplit(temp.JRNe, "_"))[5] ## THIS IS AT MAF 0.05 CUTOFF
	temp.JRNe.lci <- d[JorRy.line[2] + 9]
	temp.JRNe.lci <- gsub(" +", "_", temp.JRNe.lci)
	JRNe.lci <- unlist(strsplit(temp.JRNe.lci, "_"))[6]  ## THIS IS AT MAF 0.05 CUTOFF, jackknife
	temp.JRNe.uci <- d[JorRy.line[2] + 10]
	temp.JRNe.uci <- gsub(" +", "_", temp.JRNe.uci)
	JRNe.uci <- unlist(strsplit(temp.JRNe.uci, "_"))[2]  ## THIS IS AT MAF 0.05 CUTOFF

	file.name <- tail(unlist(strsplit(i, "/")), n=1)
	
	# order will be:  filename, pop1 LD & CIs, pop2 LD &CIs, pop1 het & CIs, pop2 het & CIs, pop1 coanc & CIs, pop2 coanc & CIs, temp Pollak & CIs, temp NeiTaj & CIs, temp JorRy & CIs
	answer <- paste(file.name, LDNe1, LDNe1.lci, LDNe1.uci, LDNe2, LDNe2.lci, LDNe2.uci, hetNe1, hetNe1.lci, hetNe1.uci, hetNe2, hetNe2.lci, hetNe2.uci, coancNe1, coancNe1.lci, coancNe1.uci, coancNe2, coancNe2.lci, coancNe2.uci, PNe, PNe.lci, PNe.uci, NTNe, NTNe.lci, NTNe.uci, JRNe, JRNe.lci, JRNe.uci, sep=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/NeEstimator_Estimates"))
