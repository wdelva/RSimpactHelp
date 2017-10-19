#loading packages
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, magrittr, dplyr,
               ggplot2, mice, tidyr, lhs)

dirname <- "~/Documents/GIT_Projects/RSimpactHelp"

#dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"
#dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"

#simpact simulation traing dataset
training.df <- "SummaryOutPut-df-iSNweKFVoZ.csv"

#Read the saved data from the Simpact simulation
complete.df.wna <- data.frame(read.csv(file = paste0(dirname,"/",training.df), header = TRUE))

#remove all the NA and Nan values from the resulting simulation
complete.df.wna <- subset(complete.df.wna, complete.cases(complete.df.wna))

#Just to get the varied parameter names
source("R/Misc/SimpactWrapper/simpact.parameters.R")
inPUT.df.complete <- head(source.simpact.parameters(),1)
input.parameter.names <- names(dplyr::select(inPUT.df.complete, contains(".")))
input.parameter.names <- input.parameter.names[2:length(input.parameter.names)]

#get the target values and their names
source("R/Misc/SimpactWrapper/summary.statistics.creator.R")

#selection data for calibration
filter.names <- c(target.variables, input.parameter.names)

### target parameters with NA in simpact varied parameters
target.parameter.na <- c(target.values, rep(NA, length(input.parameter.names)))

#get only the df with the Target values and Simpact waried parameters
complete.df.wna <- complete.df.wna[ , (names(complete.df.wna) %in% filter.names)]

########      Now create the df for mice imputation #########################
complete.df.wna <- rbind(complete.df.wna, target.parameter.na)


mice.calibration <- function(df.wna, method.id = 1, out.sets = 5, iter.num = 5, seed.id = 2525){

  #selected methods used
  cal.list <- c("pmm", "norm", "norm.boot")

  #calibration
  mice.imputted.df <- mice(df.wna, m=out.sets, maxit = iter.num, method = cal.list[method.id],
                           seed = seed.id)

  #collect all the imputed sets
  mice.out.imputted.df <- tail(mice::complete(mice.imputted.df,1),1)

  for (i in 2:out.sets){
    complete.mice.wo.na <- mice::complete(mice.imputted.df,i)
    mice.out.imputted.df <- rbind(mice.out.imputted.df, tail(complete.mice.wo.na,1))
  }

  #remove all the rows with out of range imputed files
  mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.man.dist.normal.sigma > 0)
  mice.out.imputted.df <- subset(mice.out.imputted.df, hivtransmission.param.f1 > 0)

  mice.out.imputted.df$method <- cal.list[method.id]

  return(mice.out.imputted.df)

}

#check time
mice.startt <- as.numeric(proc.time()[3])
mice.calibration <- mice.calibration(complete.df.wna, method.id = 1, out.sets = 10)
save.string <- paste0(sample(c(LETTERS,letters), 5), collapse = "")
write.csv(mice.calibration, file = paste0(dirname,"/mice.",cal.list[method.id],".",save.string,
                                          ".",training.df), row.names = FALSE)

mice.end.pmm <- as.numeric(proc.time()[3]) - mice.startt

