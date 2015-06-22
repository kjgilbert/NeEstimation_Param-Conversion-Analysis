#read file in

setwd("  ")


all.files <- list(
"MLNeOUT_SmallerSample_P5_mig50_meta0.01_gen115_t2_1_001.dat","MLNeOUT_SmallerSample_P5_mig50_meta0.01_gen115_t2_1_002.dat"
)


results <- NULL
for(i in all.files){
	d <- readLines(i)
	
	# first, did the file also estimate migration?
	str1 <- "Estimating"
	Ne.lik <- grep(str1, d)
	temp <- unlist(strsplit(d[Ne.lik], " "))
	is.isolated <- "only," %in% temp
	
	if(is.isolated==TRUE){ # when migration is not estimated
		str.lik <- "Likelihood estimates"
		lik.dat <- (grep(str.lik, d))+2  #data is 2 lines down from the title of the section
		Ne.lik.temp <- gsub(" +", "_", d[lik.dat])
		Ne.lik <- unlist(strsplit(Ne.lik.temp, "_"))[3]
		uci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[4]
		lci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[5]
		if(uci.lik==">"){ # if it was an infinite limit
			uci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[5]
			lci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[6]
		}
		
		str.mom <- "Moment     estimates"
		mom.dat <- (grep(str.mom, d))+2  #data is 2 lines down from the title of the section
		Ne.mom.temp <- gsub(" +", "-", d[mom.dat])
		Ne.mom.temp.temp <- which(unlist(strsplit(Ne.mom.temp, "-"))=="MT_Ne=")
		Ne.mom <- unlist(strsplit(Ne.mom.temp, "-"))[Ne.mom.temp.temp+1]
	}
	if(is.isolated==FALSE){ # when migration is estimated
		str.lik <- "Likelihood estimates"
		lik.dat <- (grep(str.lik, d))+2  #data is 2 lines down from the title of the section
		Ne.lik.temp <- gsub(" +", "_", d[lik.dat])
		Ne.lik <- unlist(strsplit(Ne.lik.temp, "_"))[3]
		uci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[4]
		lci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[5]
		if(uci.lik==">"){ # if it was an infinite limit
			uci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[5]
			lci.lik <- unlist(strsplit(Ne.lik.temp, "_"))[6]
		}
		mig.lik.temp <- gsub(" +", "_", d[lik.dat+1]) # migration data is next line down
		m.lik <- unlist(strsplit(mig.lik.temp, "_"))[3]
		uci.m.lik <- unlist(strsplit(mig.lik.temp, "_"))[4]
		lci.m.lik <- unlist(strsplit(mig.lik.temp, "_"))[5]
		if(lci.m.lik=="<"){ # if it was an infinite limit
			lci.m.lik <- unlist(strsplit(mig.lik.temp, "_"))[6]
		}
		if(uci.m.lik==">"){ # if it was an infinite limit
			uci.m.lik <- unlist(strsplit(mig.lik.temp, "_"))[5]
			lci.m.lik <- unlist(strsplit(mig.lik.temp, "_"))[6]
		}
		
		str.mom <- "Moment     estimates"
		mom.dat <- (grep(str.mom, d))+2  #data is 2 lines down from the title of the section
		Ne.mom.temp <- gsub(" +", "`", d[mom.dat])
		Ne.mom.temp.temp <- which(unlist(strsplit(Ne.mom.temp, "`"))=="MT_Ne=")
		Ne.mom <- unlist(strsplit(Ne.mom.temp, "`"))[Ne.mom.temp.temp+1]
		if(length(Ne.mom)==0){
			Ne.mom <- "inf"
		}
		m.mom.temp <- which(unlist(strsplit(Ne.mom.temp, "`"))=="MT_m=")
		m.mom <- unlist(strsplit(Ne.mom.temp, "`"))[m.mom.temp+1]
	}
	
	name <- tail(unlist(strsplit(i, "/")), n=1)
	
	## NOTE upper CI comes first!
	if(is.isolated==TRUE){ # when migration is estimated
		answer <- paste(name, Ne.lik, uci.lik, lci.lik, Ne.mom, sep=" ")
	}
	if(is.isolated==FALSE){ # when migration is estimated
		answer <- paste(name, Ne.lik, uci.lik, lci.lik, Ne.mom, m.lik, uci.m.lik, lci.m.lik, m.mom, sep=" ")
	}
	results <- paste(c(paste(results), paste(answer), sep=" "))
}

writeLines(results, paste("~/Documents/My_Documents/UBC/Research/Ne_SimulationProject/Ne_Estimations/GrepResults/MLNe_Estimates"))
