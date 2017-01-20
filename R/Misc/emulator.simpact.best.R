#Read the data from the emulator best 221.
emulator.simpact.analysis <- inputANDoutput.chunk.df
summaryparameters.best <- summaryparameters
targets <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)


file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-1-274-2017-01-18.csv") # param.varied

##file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-Emu222-274-2017-01-19.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.completeReminder <- data.frame(read.csv(file = file.name.csv, header = TRUE))

emulator.simpact.analysis.join <- rbind(emulator.simpact.analysis, inputANDoutput.completeReminder)


inputANDoutput.analysis <- inputANDoutput.completeReminder

#get all the average in all the columns in the selected df
inputANDoutput.analysis <- aggregate(emulator.simpact.analysis.join, by = list(emulator.simpact.analysis.join$sim.id), FUN = "mean")

inputANDoutput.analysis <- aggregate(inputANDoutput.analysis, by = list(inputANDoutput.analysis$sim.id), FUN = "mean")

inputANDoutput.analysis <- subset(inputANDoutput.analysis, select=summaryparameters)

#### This will generate pairs that with each parameter
inputANDoutput.analysis$is.complete <- complete.cases(inputANDoutput.analysis)
sum(inputANDoutput.analysis$is.complete)

x.lim <-list(c(-0.02,0.04),c(0,0.3),c(0,0.07),c(0,0.9),
             c(0.1,0.9),c(0,1),c(0.1,0.9), c(0.1,0.45), c(0.1,0.45),c(3.5,6.5))
n.bars <-c(13,6,13,9,8,10,8,10,8,6)


inputANDoutput.analysis$inc.men.ok <- (inputANDoutput.analysis$inc.men.20.25 >0.011 &
                                         inputANDoutput.analysis$inc.men.20.25 <0.025)

inputANDoutput.analysis$inc.wom.ok <- (inputANDoutput.analysis$inc.wom.20.25 >0.033 &
                                         inputANDoutput.analysis$inc.wom.20.25 <0.056)

inputANDoutput.analysis$inc.complete <- inputANDoutput.analysis$inc.men.ok*inputANDoutput.analysis$inc.wom.ok

inputANDoutput.analysis$met.cat <- inputANDoutput.analysis$is.complete


inputANDoutput.analysis$met.cat[inputANDoutput.analysis$is.complete==TRUE &
                                  inputANDoutput.analysis$inc.complete==TRUE] <- 3

inputANDoutput.analysis.p <- inputANDoutput.analysis

inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat>2] <-
  inputANDoutput.analysis.p$sim.id[inputANDoutput.analysis.p$met.cat>2]

inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==6] <- 3
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==187] <- 4
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==54] <- 5

pair.plot <- pairs(inputANDoutput.analysis.p[, x.variables[1:4]],
                   col = 1+inputANDoutput.analysis.p$met.cat, pch = 16, cex = 2)


inputANDoutput.analysis <- subset(inputANDoutput.analysis, is.complete == TRUE)

simpact.z.analysis <- dplyr::select_(inputANDoutput.analysis,.dots=z.variables)
par(mfrow=c(1,1))
multi.hist(simpact.z.analysis)

sum(simpact.z.analysis$inc.complete)

par(mfrow=c(4,3)) # distribution of statistics
j <- 0
for (i in z.variables){
  j <- j + 1
  hist(simpact.z.analysis[,i], n.bars[j], xlim = x.lim[[j]],
       main = paste("Dist of ",i, sep = " "), xlab = "simpact.statistic")
  abline(v = as.numeric(targets[j]), col="red", lwd=3, lty=2)
}

















pair.plot <- pairs(inputANDoutput.analysis[, x.variables[1:8]],
                   col = 1+inputANDoutput.analysis$is.complete, pch = 16, cex = 2)

x.lim <-list(c(-0.02,0.04),c(0,0.3),c(0,0.07),c(0,0.9),
             c(0.1,0.9),c(0,1),c(0.1,0.9), c(0.1,0.45), c(0.1,0.45),c(3.5,6.5))
n.bars <-c(13,6,13,9,8,10,8,10,8,6)


inputANDoutput.analysis$inc.men.ok <- (inputANDoutput.analysis$inc.men.20.25 >0.011 &
                                         inputANDoutput.analysis$inc.men.20.25 <0.025)

inputANDoutput.analysis$inc.wom.ok <- (inputANDoutput.analysis$inc.wom.20.25 >0.033 &
                                         inputANDoutput.analysis$inc.wom.20.25 <0.056)

inputANDoutput.analysis$inc.complete <- inputANDoutput.analysis$inc.men.ok*inputANDoutput.analysis$inc.wom.ok

inputANDoutput.analysis$met.cat <- inputANDoutput.analysis$is.complete


inputANDoutput.analysis$met.cat[inputANDoutput.analysis$is.complete==TRUE &
                                  inputANDoutput.analysis$inc.complete==TRUE] <- 3

inputANDoutput.analysis.p <- inputANDoutput.analysis

inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat>2] <-
  inputANDoutput.analysis.p$sim.id[inputANDoutput.analysis.p$met.cat>2]

inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==6] <- 3
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==187] <- 4
inputANDoutput.analysis.p$met.cat[inputANDoutput.analysis.p$met.cat==54] <- 5

pair.plot <- pairs(inputANDoutput.analysis.p[, x.variables[1:8]],
                   col = 1+inputANDoutput.analysis.p$met.cat, pch = 16, cex = 2)


inputANDoutput.analysis <- subset(inputANDoutput.analysis, is.complete == TRUE)

simpact.z.analysis <- dplyr::select_(inputANDoutput.analysis,.dots=z.variables)
par(mfrow=c(1,1))
multi.hist(simpact.z.analysis)

sum(simpact.z.analysis$inc.complete)

par(mfrow=c(4,3)) # distribution of statistics
j <- 0
for (i in z.variables){
  j <- j + 1
  hist(simpact.z.analysis[,i], n.bars[j], xlim = x.lim[[j]],
       main = paste("Dist of ",i, sep = " "), xlab = "simpact.statistic")
  abline(v = as.numeric(targets[j]), col="red", lwd=3, lty=2)
}



