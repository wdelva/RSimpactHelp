#loading packages
pacman::p_unload(all)
pacman::p_load(magrittr, dplyr, EasyABC, foreach, mice,
               doParallel, lhs)

# comp <- "lin" #lin #mac #chpc #gent
#
# if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
#   dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
#     dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
#       dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
#         dirname <- "~/Documents/RSimpactHelp"  #mac directory here
#       }

#local source file
dirname <- "~/Documents/100000.Pre.Hhoho.OutFiles/"

##simpact simulation traing dataset analysed and saved data
#training.df <- "Cal-Filtered-1e5.csv" #original data
#training.df <- "pmm1.Cal-Filtered-1e5.csv"  #pmm1 first calibration
#training.df <- "pmm2.ppm1.Cal-Filtered-1e5.csv"  #pmm2 sec calibration
training.df <- "pmm3.pmm2.ppm1.Cal-Filtered-1e5.csv"  #pmm2 sec calibration


##Read the saved data from the Simpact simulation
training.input.df <- data.frame(read.csv(file = paste0(dirname, training.df), header = TRUE))
#call simpact.wrapper.model
source("R/Misc/SimpactCalibration/mice.calibration.Template.R")

#check time
calibrated.datalist  <- mice.calibration(training.input.df, method.id = 1, iter.num = 10,
                                         cores_2_use = 4, out.sets = 100, printflag = TRUE)

#Output is saved in dirname.





