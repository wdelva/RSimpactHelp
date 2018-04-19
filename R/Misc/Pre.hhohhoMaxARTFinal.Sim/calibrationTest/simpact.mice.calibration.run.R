#To run after obtaining the training dataset
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

#source simpact simulation function
#source("R/Misc/pbdMPI.hhohho/pbdmpi.simpact.simulation.R")

#Source all latest RSimpactHelper functions
#source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/hho.rsimpacthelper.R")

#call the calibration methods wrapper
source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/mice.calibration.R")

#simpact simulation TRAINING dataset will be read from file
training.df.file <- "Method.mpi.simid.fin.2001-4-SummaryOutPut-haRoMymzcf-pdbMPI.csv"
method.list.id <- c("pmm", "norm", "rf")

#Read the saved data from the Simpact simulation
training.df <- data.frame(read.csv(file = paste0(dirname,"/R/Misc/Pre.hhohhoMaxARTFinal.Sim/AdjAlteredMICE/",
                                                 training.df.file), header = TRUE))


#training.df$method <- "Init.df"
#write.csv(training.df, file = paste0(dirname,"/R/Misc/Pre.hhohhoMaxARTFinal.Sim/AdjAlteredMICE/Method.",
#                                     training.df.file), row.names = FALSE)

#Get only the rows that do not have NA from simpact simulations.
training.df <- subset(training.df, complete.cases(training.df))

training.df <- training.df %>%
  filter(year10.X1986.gr > 0, year10.X1997.gr > 0,
         year10.X2007.gr > 0, year10.X2012.gr > 0,
         meadian.FM.value.agediff > 0, mean.FM.value.agediff > 0,
         meadian.F.value.agediff > 0, mean.M.value.agediff > 0,
         meadian.M.value.agediff > 0 )


# training.df$method <- as.character(training.df$method)
#
# training.df$method[training.df$method == "pmm.fil"] <- "pmm"
# training.df$method[training.df$method == "pmm.all"] <- "pmm"
# training.df$method[training.df$method == "rf.all"] <- "rf"
# training.df$method[training.df$method == "rf.fil"] <- "rf"


save.name <- "filt.MICEpdbMPI."
#Calibrate on all the methods as a single df
method.val <- 1 #option, 1,2,3 c("pmm", "norm", "rf")
calibrated.datalist.1  <- mice.calibration(training.df, method.id = method.val, out.sets = 150)

#method.val <- 2
#calibrated.datalist.2  <- mice.calibration(training.df, method.id = method.val, out.sets = 150)

method.val <- 3
calibrated.datalist.3  <- mice.calibration(training.df, method.id = method.val, out.sets = 150)

# #Calibrate on the initial method of calibration
#
# save.name <- "MICEpdbMPI.mfilt."
#
# method.val <- 1
# training.df.f <- subset(training.df, method == method.list.id[method.val] )
# calibrated.datalist.1  <- mice.calibration(training.df.f, method.id = method.val, out.sets = 150)
#
# #method.val <- 2  #Not enough data to do the calibration further
# #training.df.f <- subset(training.df, method == method.list.id[method.val] )
# #calibrated.datalist.2  <- mice.calibration(training.df.f, method.id = method.val, out.sets = 150)
#
# method.val <- 3
# training.df.f <- subset(training.df, method == method.list.id[method.val] )
# calibrated.datalist.3  <- mice.calibration(training.df.f, method.id = method.val, out.sets = 150)
#

# #number of calibrations to do
# n.fin <- 3
# n.val <- 0
# save.name <- "pbdmpi.train0.cal." # ""
# method.val <- 1
#
# print("Simulation started")
#
# training.df <- pbdMPI.simulation.wrapper(set.seed.id = 1,   #what seed to use
#                                                  design.points = 100, #
#                                                  par.repeat = 1, #each row is repeated once.
#                                                  min.sim = 1, max.sim = 100,
#                                                  datalist = NA,  #Not set initially
#                                                  cal.simulation = FALSE)
#
# #close all the cores
# finalize()
#
# print("Initial simulation ended:")
# while(n.val < n.fin){
#
#   save.name <- paste0("MICEpdbMPI.pmm.",n.val,".")
#   print(paste0(n.val, " :caibration started"))
#   calibrated.datalist <- mice.calibration(training.df, method.id = method.val, out.sets = 100)
#
#
#   if(nrow(calibrated.datalist) > 2){
#     #calibration returned enough rows to continue simulation
#
#     training.df <- pbdMPI.simulation.wrapper(set.seed.id = 1,   #what seed to use
#                                              design.points = NA, # will use the nrow(datalist)
#                                              par.repeat = 1, #each row is repeated once.
#                                              min.sim = NA, max.sim = NA, #will use all rows
#                                              datalist = calibrated.datalist,
#                                              cal.simulation = TRUE)
#     #close all the cores
#     finalize()
#
#     #increment the cal number
#     n.val <- n.val + 1
#
#   }else{
#     n.val <- 4
#     save.name <- "NoOutputMICEpdbMPI"
#     #save the MICE inputfile
#     rand.string.fin <- paste0(sample(c(LETTERS,letters), 10), collapse = "")
#     filename.run.fun <- paste0(dirname,"/","ModelCalibrationOutput-df-",rand.string.fin,"-",
#                                method.list.id, "-", save.name,".csv")
#     write.csv(training.df, file = filename.run.fun, row.names = FALSE)
#
#   }
#
# }











