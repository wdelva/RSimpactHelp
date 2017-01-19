#Read the data from the emulator best 221.
emulator.simpact.analysis <- inputANDoutput.chunk.df
summaryparameters.best <- summaryparameters
targets <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)

file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-Emu222-274-2017-01-19.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.completeReminder <- data.frame(read.csv(file = file.name.csv, header = TRUE))

emulator.simpact.analysis.join <- rbind(emulator.simpact.analysis, inputANDoutput.completeReminder)

#get all the average in all the columns in the selected df
inputANDoutput.analysis <- aggregate(emulator.simpact.analysis.join, by = list(emulator.simpact.analysis.join$sim.id), FUN = "mean")

#### This will generate pairs that with each parameter
inputANDoutput.analysis$is.complete <- complete.cases(inputANDoutput.analysis)
sum(inputANDoutput.analysis$is.complete)

pair.plot <- pairs(inputANDoutput.analysis[, x.variables[1:8]],
                   col = 1+inputANDoutput.analysis$is.complete, pch = 16, cex = 2)

inputANDoutput.analysis <- subset(inputANDoutput.analysis, is.complete == TRUE)


simpact.z.analysis <- dplyr::select_(inputANDoutput.analysis,.dots=z.variables)
par(mfrow=c(1,1))
multi.hist(simpact.z.analysis)

par(mfrow=c(4,3)) # distribution of statistics
j <- 0
for (i in z.variables){
  j <- j + 1
  hist(simpact.z.analysis[,i], main = paste("Dist of ",i, sep = " "), xlab = "simpact.statistic")
  abline(v = as.numeric(targets[j]), col="red", lwd=3, lty=2)
}



