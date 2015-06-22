#read file in

setwd("  ")



all.files <- list(
"",""
)


results <- NULL
for(i in all.files){
	d <- readLines(i)
	str1 <- "MaxByParabolicInterpolation"
	Ne.line <- grep(str1, d)

	if(length(Ne.line)==0){
		Ne <- "fail"
		lci <- "fail"
		uci <- "fail"
	}else{   # IF NO OUTPUT PRODUCED, Ne.line is an empty vector because that text doesn't exist, so instead print this to the output
	
	# just pull out the last few lines I want
	data <- d[Ne.line:(Ne.line+2)]

	Ne <- unlist(strsplit(data[1], "  "))[2]
	lci <- unlist(strsplit(data[2], " "))[4]
	uci <- unlist(strsplit(data[3], " "))[4]
}
	name <- tail(unlist(strsplit(i, "/")), n=1)

	answer <- paste(name, Ne, lci, uci, sep=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/CoNe_Estimates"))
