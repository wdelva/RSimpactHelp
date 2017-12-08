#This function will calibrate from that simulated dataset

mice.calibration <- function(simpact.sim.df, method.id = 1, out.sets = 5, iter.num = 5, seed.id = 2525){

  #Just to get the varied parameter names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.parameters.R")
  inPUT.df.complete <- head(source.simpact.parameters(),1)
  input.parameter.names <- names(dplyr::select(inPUT.df.complete, contains(".")))
  input.parameter.names <- input.parameter.names[2:length(input.parameter.names)]

  #get the target values and their names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")

  #selection data for calibration
  filter.names <- c(target.variables, input.parameter.names)

  ### target parameters with NA in simpact varied parameters
  target.parameter.na <- c(target.values, rep(NA, length(input.parameter.names)))

  #Get only the rows that do not have NA from simpact simulations.
  simpact.sim.df <- subset(simpact.sim.df, complete.cases(simpact.sim.df))

  if(nrow(simpact.sim.df) > 10){

    #get only the df with the Target values and Simpact varied parameters
    simpact.sim.df <- simpact.sim.df[ , (names(simpact.sim.df) %in% filter.names)]

    ########      Now create the df for mice imputation #########################
    complete.df.wna <- rbind(simpact.sim.df, target.parameter.na)


    #selected methods used
    cal.list <- c("pmm", "norm", "norm.boot")

    #calibration
    mice.imputted.df <- mice(complete.df.wna, m=out.sets, maxit = iter.num,
                             method = cal.list[method.id], seed = seed.id)

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

    #save the data if needed.
    save.string <- paste0(sample(c(LETTERS,letters), 5), collapse = "")
    write.csv(mice.calibration, file = paste0(dirname,"/mice.",cal.list[method.id],".",save.string,
                                              ".",training.df), row.names = FALSE)
  }else{

    #If not enough simpact returned dataset
    mice.out.imputted.df <- as.data.frame(matrix(NA, 0, length(target.parameter.na)))
    names(mice.out.imputted.df) <- filter.names

  }

  return(mice.out.imputted.df)

}

