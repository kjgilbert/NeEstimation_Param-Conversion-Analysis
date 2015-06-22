source('MC_Script_AnalyzeTMVP.R', chdir = TRUE)

setwd("  ")


all.files <- list(

results <- NULL
for(i in all.files){
	filename <- tail(unlist(strsplit(i, "/")), n=1)
	dat <- read.table(i)
	raw.answer <- mc(dat, maxk=100)
	answer <- paste(filename, paste(raw.answer, collapse=" "), collapse=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/TMVP_Estimates"))

