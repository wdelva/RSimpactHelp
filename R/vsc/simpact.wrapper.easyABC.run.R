#get the necessary libraries
pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr,lhs,data.table, dplyr, magrittr, exactci,
             nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
             igraph,lhs, GGally, emulator, multivator, tidyr)

#data file to read 
dirname <- getwd()
main.filename <- "INPUT.df-10Points10Par2016-10-31.csv" #Read the file produced by varying parameters *design.points
file.chunk.name.csv <-paste0(dirname, "/", main.filename) #### Input file name is produced from the .sh script
inANDout.df.complete <- read.csv(file = file.chunk.name.csv, header = TRUE, sep = ",")    


#Select a chunk to send to process
min.chunk <- 0
max.chunk <- 10
inANDout.df.chunk <- inANDout.df.complete[min.chunk:max.chunk,]

sim_repeat <- 10
ncluster.use <- 5

## In case you do not need all the target statistics
target.variables <-c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25", "prev.men.25.50",
                     "ART.cov.15.50", "incid.wom.15.30", "frac.degreeGT1.wom.15.30", "mean.degree",
                     "median.degree", "Q1.degree", "Q3.degree")

#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 1, ncol = length(target.variables)))
names(chunk.summary.stats.df) <- target.variables
chunk.summary.stats.df$sim.id <- NA

start.chunk.time <- proc.time()
for (chunk.sim.id in inANDout.df.chunk$sim.id){

  simpact.chunk.prior = list()
  
  for (i in preprior.names.chunk){
    
    #col.index <- which(colnames(preprior.names.chunk)==i)
    
    prior.chunk.val <- list(c("runif",1,as.numeric(inANDout.df.chunk[chunk.sim.id,i]),
                              as.numeric(inANDout.df.chunk[chunk.sim.id,i])), c("dunif",0,1))
    simpact.chunk.prior[[length(simpact.chunk.prior)+1]] <- prior.chunk.val
  }
  
  #invoke the ABC_rejection method repeating the number of simulation X* for each chunk row.
  ABC.chunk.result <- ABC_rejection(model = simpact.chunk.run,
                                        prior = simpact.chunk.prior,
                                        nb_simul= sim_repeat,
                                        use_seed = TRUE,
                                        seed_count = 0,
                                        n_cluster = ncluster.use)
  
  #Save the statistics results with the chunk row sim.id repeated X* from the ABC_rejection method 
  ABC.results.chunk.statistics <- data.frame(ABC.chunk.result$stats)
  names(ABC.results.chunk.statistics) <- target.variables
  ABC.results.chunk.statistics$sim.id <- chunk.sim.id
  
  chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)
  
}

write.csv(chunk.summary.stats.df, file =paste0("SummaryOutPut-inANDout.df.chunk-",min.chunk,"-",max.chunk,"-",Sys.Date(), 
                                               ".csv"), row.names = FALSE)

end.chunk.time <- proc.time() - start.chunk.time














