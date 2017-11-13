#loading packages
#pacman::p_unload(all)
pacman::p_load(RSimpactCyan, data.table, magrittr, dplyr, lhs, tidyr,
               exactci, readcsvcolumns, foreach, doParallel)

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/output"}else if(comp=="gent"){
      dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
        dirname <- "~/Documents/RSimpactHelp"  #mac directory here
      }

###call simpact.wrapper.model
source("R/Misc/foreach.Hhohho/foreach.simpact.simulation.R")
#local source file#
#dirname <- "~/Documents/100000.Pre.Hhoho.OutFiles"
##simpact simulation traing dataset analysed and saved data
#cal.datalist <- "PreHhohho.Repeat.1e5.SSE.df.csv"
save.name <- "foreach.test1" # ""
##Read the saved data from the Simpact simulation
#cal.datalist <- data.frame(read.csv(file = paste0(dirname,"/", cal.datalist), header = TRUE))

#cal.datalist <- head(cal.datalist, 4)

#Run the simulation.
training.df <- foreach.simulation.wrapper(set.seed.id = 1,   #what seed to use
                                          design.points = 4, #
                                          par.repeat = 2, #each row is repeated once.
                                          ncluster.use = 4,
                                          min.sim = 1, max.sim = 4,
                                          datalist = NA, #cal.datalist  #Not set initially
                                          cal.simulation = FALSE)

#We can now do the analysis.





