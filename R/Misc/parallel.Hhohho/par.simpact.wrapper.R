#loading packages
#pacman::p_unload(all)
pacman::p_load(RSimpactCyan, data.table, magrittr, dplyr, lhs, tidyr,
               exactci, readcsvcolumns, parallel)

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/GIT_Projects/RSimpactHelp/"}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp/"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/output/"}else if(comp=="gent"){
      dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data/"}else{
        dirname <- "~/Documents/RSimpactHelp/"  #mac directory here
      }

## Set working directory
setwd("~/GIT_Projects/RSimpactHelp/")

###call simpact.wrapper.model
source("R/Misc/parallel.Hhohho/par.simpact.simulation.R")
#local source file#
#dirname <- "~/Documents/100000.Pre.Hhoho.OutFiles"

#indicate if initial simulation or from a file
sim.initiate <- TRUE #TRUE
ncluster.use <- 4

if(sim.initiate == TRUE){
  #change these to reflect the file to read
  ##simpact simulation traing dataset analysed and saved data
  cal.datalist <- "Final.Selected.Target.Simulation.SSE.df.csv"

  #Read the saved data from the Simpact simulation
  cal.datalist <- data.frame(read.csv(file = paste0(dirname, cal.datalist), header = TRUE))
  save.name <- paste0("par.simid.cal.",nrow(cal.datalist)) #file savename
  dp.val <- NA #No necessary replaced by cal.datalist
  rp.val <- 4 #default no repeat
  mins.val <- NA #Not necessary
  maxs.val <- NA #Not necessary

}else{

  dp.val <- 20 #create parameter space
  save.name <- paste0("par.sim.init.", dp.val) #file save name
  rp.val <- 1 #default no repeat
  cal.datalist <- NA #no necessary created by dp.val
  mins.val <- 1 #select from the parameter space min
  maxs.val <- 4 #select from teh parameter space max

}

#Run the simulation.
training.df <- par.simulation.wrapper(set.seed.id = 1,   #what seed to use
                                      design.points = dp.val, #
                                      par.repeat = rp.val, #each row is repeated once.
                                      ncluster.use = ncluster.use,
                                      min.sim = mins.val, max.sim = maxs.val,
                                      datalist = cal.datalist, #cal.datalist  #Not set initially
                                      cal.simulation = sim.initiate)

#We can now do the analysis.





