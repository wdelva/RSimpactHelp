#loading packages
#pacman::p_unload(all)
pacman::p_load(RSimpactCyan, data.table, magrittr, dplyr, EasyABC, foreach,
               mice, tidyr, lhs, tidyr, exactci, readcsvcolumns, pbdMPI, doParallel)

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
      dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
        dirname <- "~/Documents/RSimpactHelp"  #mac directory here
      }


#Initial simulation has not started
sim.cal.count <- 0
#Number of calibration runs
model.cal.count <- 2

################################################## IF DATA WAS SIMULATED BEFORE #########
##simpact simulation traing dataset (if saved data)
#sim <- 1
#training.df <- "SummaryOutPut-df-iSNweKFVoZ.csv"
##Read the saved data from the Simpact simulation
#training.df <- data.frame(read.csv(file = paste0(dirname,"/",training.df), header = TRUE))
##remove all the NA and Nan values from the resulting simulation
#complete.df.wna <- subset(training.df, complete.cases(training.df))

#call simpact.wrapper.model
source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.wrapper.easyABC.run.R")
source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/mice.calibration.R")


for (k in 1:model.cal.count){

  if(sim.cal.count == 0){

    training.df <- easyABC.simulation.wrapper(sim.seed = 1,   #what seed to use
                                              design.points = 10, #
                                              par.repeat = 1, #each row is repeated once.
                                              min.sim = 1, max.sim = 4,
                                              datalist = NA, #Not set initially
                                              cal.simulation = FALSE)
  }else{

    #check time
    calibrated.datalist  <- mice.calibration(training.df, method.id = 1, out.sets = 10)

      if(nrow(calibrated.datalist) > 5){
        #calibration returned enough rows to continue simulation

        training.df <- easyABC.simulation.wrapper(sim.seed = 1,   #what seed to use
                                                  design.points = 10, #
                                                  par.repeat = 1, #each row is repeated once.
                                                  min.sim = NA, max.sim = NA, #will use all rows
                                                  datalist = calibrated.datalist,
                                                  cal.simulation = TRUE)
      }else{
        sim.cal.count <- model.cal.count
      }

  }

  sim.cal.count <- sim.cal.count + 1

}

#Write the final dataframe for analysis.
rand.string.fin <- paste0(sample(c(LETTERS,letters), 10), collapse = "")
filename.run.fun <- paste0(dirname,"/","ModelOutPutCalibrated-df-",rand.string.fin,"-MICEpdbMPI.csv")
write.csv(inputANDoutput.chunk.df, file = filename.run.fun, row.names = FALSE)

#We can now do the analysis.





