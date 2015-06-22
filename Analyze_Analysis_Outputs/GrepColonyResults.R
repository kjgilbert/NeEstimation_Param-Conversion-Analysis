#read file in

setwd("  ")

all.files <- list(
"ColonyOUT_P1_mig50_meta0.1_gen015_t2_1_004.dat.Ne","ColonyOUT_P1_mig50_meta0.1_gen015_t2_1_080.dat.Ne"
)



results <- NULL
for(i in all.files){
	d <- readLines(i)
	str1 <- "----------------------------------------------------------------------------------------------------"
	Ne.line.sib <- grep(str1, d)[2] #get the 2nd line which is just above the Ne estimates based on sibship data
	Ne.line.het <- grep(str1, d)[4] #get the 2nd line which is just above the Ne estimates based on sibship data

	# just pull out the last few lines I want
	data.sib <- d[(Ne.line.sib+2):(Ne.line.sib+4)]
	data.het <- d[(Ne.line.het+1):(Ne.line.het+3)]

	Ne.sib.temp <- gsub(" +", "_", data.sib)
	Ne.sib <- unlist(strsplit(Ne.sib.temp, "_"))[3]
	lci.sib <- unlist(strsplit(Ne.sib.temp, "_"))[6]
	uci.sib <- unlist(strsplit(Ne.sib.temp, "_"))[9]
	
	Ne.het.temp <- gsub(" +", "_", data.het)
	Ne.het <- unlist(strsplit(Ne.het.temp, "_"))[3]
	lci.het <- unlist(strsplit(Ne.het.temp, "_"))[6]
	uci.het <- unlist(strsplit(Ne.het.temp, "_"))[9]
	
	name <- tail(unlist(strsplit(i, "/")), n=1)
	
	answer <- paste(name, Ne.sib, lci.sib, uci.sib, Ne.het, lci.het, uci.het, sep=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/Colony_Estimates"))
