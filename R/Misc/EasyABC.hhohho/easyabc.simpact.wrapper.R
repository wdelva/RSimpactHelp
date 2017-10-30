#loading packages
pacman::p_load(RSimpactCyan, data.table, magrittr, dplyr, EasyABC,
               mice, tidyr, lhs, tidyr, exactci, readcsvcolumns, pbdMPI)

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
      dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
        dirname <- "~/Documents/RSimpactHelp"  #mac directory here
      }

###call simpact.wrapper.model
source("R/Misc/EasyABC.hhohho/easyabc.simpact.simulation.R")

training.df <- easyABC.simulation.wrapper(sim.seed = 1,   #what seed to use
                                          design.points = 4, #
                                          par.repeat = 1, #each row is repeated once.
                                          ncluster.use = 1, #how many cores to use
                                          min.sim = 1, max.sim = 4,
                                          datalist = NA, #Not set initially
                                          cal.simulation = FALSE)


#Write the final dataframe for analysis.
#rand.string.fin <- paste0(sample(c(LETTERS,letters), 10), collapse = "")
#filename.run.fun <- paste0(dirname,"/","ModelOutPutSimulated-df-",rand.string.fin,"-pdbMPI.csv")
#write.csv(training.df, file = filename.run.fun, row.names = FALSE)

#We can now do the analysis.





