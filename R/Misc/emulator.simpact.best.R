#Read the data from the emulator best 274.
summaryparameters.best <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                            "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                            "ART.cov.wom.18.50", "median.wom.18.50.AD")
targets <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)


file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-BESTEmulator1-274-2017-01-26.csv") # param.varied

file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-BESTEmulatorRERUN1-7-2017-01-27.csv") # param.varied

##file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-Emu222-274-2017-01-19.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.completeReminder <- data.frame(read.csv(file = file.name.csv, header = TRUE))

inputANDoutput.analysis <- inputANDoutput.completeReminder

#get all the average in all the columns in the selected df
inputANDoutput.analysis <- aggregate(inputANDoutput.analysis, by = list(inputANDoutput.analysis$sim.id), FUN = "mean")

#inputANDoutput.analysis <- subset(inputANDoutput.analysis, select=summaryparameters)

#label witch rows are complete
inputANDoutput.analysis$is.complete <- complete.cases(inputANDoutput.analysis)
sum(inputANDoutput.analysis$is.complete)

#Ploting xlim and number of bars if need be
x.lim <-list(c(-0.02,0.04),c(0,0.3),c(0,0.07),c(0,0.9),
             c(0.1,0.9),c(0,1),c(0.1,0.9), c(0.1,0.45), c(0.1,0.45),c(3.5,6.5))
n.bars <-c(13,6,13,9,8,10,8,10,8,6)

#Check which of the incidence fall within the CI
inputANDoutput.analysis$inc.men.ok <- (inputANDoutput.analysis$inc.men.20.25 >0.011 &
                                         inputANDoutput.analysis$inc.men.20.25 <0.025)

inputANDoutput.analysis$inc.wom.ok <- (inputANDoutput.analysis$inc.wom.20.25 >0.033 &
                                         inputANDoutput.analysis$inc.wom.20.25 <0.056)

#Which agree on both cases men and women incidece are within the CI
inputANDoutput.analysis$inc.complete <- inputANDoutput.analysis$inc.men.ok*inputANDoutput.analysis$inc.wom.ok

inputANDoutput.analysis$met.cat <- inputANDoutput.analysis$is.complete
#Ask which of the rows are complete and meet the CI
inputANDoutput.analysis$met.cat[inputANDoutput.analysis$is.complete==TRUE &
                                  inputANDoutput.analysis$inc.complete==TRUE] <- 3

#Want to then produce pairs plots on simpact paramters that meet the criteria
inputANDoutput.analysis.p <- inputANDoutput.analysis

inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat>2] <-
  inputANDoutput.analysis.p$sim.id[inputANDoutput.analysis.p$met.cat>2]

#Just give those rows different numbers so we see them in the pairs plot
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==191] <- 3
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==264] <- 4
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==257] <- 5
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==233] <- 6
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==183] <- 7
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==16] <- 8
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==52] <- 9

for (i in seq(1, length(x.variables), 5)) {#specific for this work
pair.plot <- pairs(inputANDoutput.analysis.p[, x.variables[i:(i+4)]],
                   col = 1+inputANDoutput.analysis.p$met.cat, pch = 16, cex = 3)
}

inputANDoutput.analysis <- subset(inputANDoutput.analysis, is.complete == TRUE)

simpact.z.analysis <- dplyr::select_(inputANDoutput.analysis,.dots=z.variables)
par(mfrow=c(1,1))
multi.hist(simpact.z.analysis)

par(mfrow=c(4,3)) # distribution of statistics
j <- 0
for (i in z.variables){
  j <- j + 1
  hist(simpact.z.analysis[,i],
       main = paste("Dist of ",i, sep = " "), xlab = "simpact.statistic")
  abline(v = as.numeric(targets[j]), col="red", lwd=3, lty=2)
}

#get the sum of square difference

check.ssd.c <- as.data.frame(t(apply(simpact.z.analysis, 1, function(x)(x-t(targets))^2)))
names(check.ssd.c) <- z.variables
check.ssd.c.plot <- multi.hist(check.ssd.c)


#with limits on the x-axis and the number of hist bars
par(mfrow=c(4,3)) # distribution of statistics
j <- 0
for (i in z.variables){
  j <- j + 1
  hist(simpact.z.analysis[,i], n.bars[j], xlim = x.lim[[j]],
       main = paste("Dist of ",i, sep = " "), xlab = "simpact.statistic")
  abline(v = as.numeric(targets[j]), col="red", lwd=3, lty=2)
}



