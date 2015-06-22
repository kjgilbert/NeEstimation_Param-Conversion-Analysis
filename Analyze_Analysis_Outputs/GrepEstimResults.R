#read file in

setwd("  ")

all.files <- list(
"EstimOUT_LONG_ideal50_gen3502_001","EstimOUT_LONG_ideal50_gen3502_002"
)


results <- NULL
for(i in all.files){
	d <- readLines(i)
	str1 <- "__________________________________________________"
	Ne.line <- grep(str1, d)[5] #get the 5th line which is just above the Ne estimates

	# just pull out the last few lines I want
	data <- d[Ne.line+1]

	Ne.temp <- gsub(" +", "_", data)
	Ne <- unlist(strsplit(Ne.temp, "_"))[3]
	Ne.lci <- unlist(strsplit(Ne.temp, "_"))[4]
	Ne.uci <- unlist(strsplit(Ne.temp, "_"))[5]
	m <- unlist(strsplit(Ne.temp, "_"))[6]
	m.lci <- unlist(strsplit(Ne.temp, "_"))[7]
	m.uci <- unlist(strsplit(Ne.temp, "_"))[8]
	
	if(is.na(m.uci)){ # id m.uci is NA, that means the data went on to the next line, grab it instead
		last.dat <- d[Ne.line+2]
		uci.temp <- gsub(" +", "_", last.dat)
		m.uci <- unlist(strsplit(uci.temp, "_"))[1]
	}
	
	# if there is an infinite Ne estimate, need to shift answers over a spot
	if(Ne == "->"){
		Ne <- unlist(strsplit(Ne.temp, "_"))[4]
		Ne.lci <- unlist(strsplit(Ne.temp, "_"))[5]
		Ne.uci <- unlist(strsplit(Ne.temp, "_"))[6]
		m <- unlist(strsplit(Ne.temp, "_"))[7]
		m.lci <- unlist(strsplit(Ne.temp, "_"))[8]
		m.uci <- unlist(strsplit(Ne.temp, "_"))[9]
		
		if(is.na(m.uci)){ # id m.uci is NA, that means the data went on to the next line, grab it instead
			last.dat <- d[Ne.line+2]
			uci.temp <- gsub(" +", "_", last.dat)
			m.uci <- unlist(strsplit(uci.temp, "_"))[1]
		}
	}
	
	name <- tail(unlist(strsplit(i, "/")), n=1)
	
	answer <- paste(name, Ne, Ne.lci, Ne.uci, m, m.lci, m.uci, sep=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/Estim_Estimates"))
