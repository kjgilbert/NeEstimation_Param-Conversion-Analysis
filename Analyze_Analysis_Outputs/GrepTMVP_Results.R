source('MC_Script_AnalyzeTMVP.R', chdir = TRUE)

setwd("  ")


all.files <- list("out_TMVP_ibd50_m0.1_P1_gen015_rep058_50000its_7500-5-0.5","out_TMVP_ibd50_m0.1_P1_gen015_rep059_50000its_7500-5-0.5")

results <- NULL
for(i in all.files){
	filename <- tail(unlist(strsplit(i, "/")), n=1)
	dat <- read.table(i)
	raw.answer <- mc(dat, maxk=100)
	answer <- paste(filename, paste(raw.answer, collapse=" "), collapse=" ")
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/TMVP_Estimates"))


