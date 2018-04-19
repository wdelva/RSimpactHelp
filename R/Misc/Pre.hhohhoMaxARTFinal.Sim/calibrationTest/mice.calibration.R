#This function will calibrate from that simulated dataset

mice.calibration <- function(simpact.sim.df, method.id = 1, out.sets = 5, iter.num = 10, seed.id = 1976){

  #Source all latest RSimpactHelper functions
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/hho.rsimpacthelper.R")
  #Just to get the varied parameter names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.parameters.R")

  inPUT.df.complete <- head(source.simpact.parameters(),1)
  input.parameter.names <- names(dplyr::select(inPUT.df.complete, contains(".")))
  input.parameter.names <- input.parameter.names[2:length(input.parameter.names)]

  #get the target values and their names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")

  #Not automated but this will reset the predictive matrix for the calibration

  # This is the default matrix.
  mice.predictorMatrix <- (1 - diag(1, length(c(input.parameter.names, target.variables))))
  # Modify the input rows , corresponding to the indicators of predictor variables for
  # the input variables. In brackets the values for the master model.

  # number of calibrated variables
  x.offset <- length(input.parameter.names)

  # First we "empty" the input parameter rows
  mice.predictorMatrix[1:x.offset, ] <- 0

  #Then we refill them: not allowing input variables to be predicted by other predictor variables.
  #Only via output variables.

  index.pM <- 0
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[1, x.offset + c(1:4, 12:13, 24:29)] <- 1 #conception.alpha_base predicted by ppn growth rate, mortality, agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 24:29)] <- 1 #dissolution.alpha_0: ppn growth rate and agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(24:29)] <- 1 #dissolution.alpha_4: agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(9:11)] <- 1 #monitoring.fraction.log_viralload: vlsup
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(5:11, 13:23)] <- 1 #person.art.accept.threshold.dist.fixed.value: ART init, ret, vlsup, prev
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 17:23)] <- 1 #person.eagerness.man.dist.gamma.a: ppn growth, men prev.
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 17:23)] <- 1 #person.eagerness.man.dist.gamma.b: ppn growth, men prev.
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[8, x.offset + c(1:4, 14:16, 20:23)] <- 1 #person.eagerness.woman.dist.gamma.a: ppn growth, wom prev.
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[9, x.offset + c(1:4, 14:23)] <- 1 #formation.hazard.agegapry.eagerness_diff: ppn growth, prev
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 14:16, 20:23)] <- 1 #person.eagerness.woman.dist.gamma.b: ppn growth, wom prev, FM prev
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 17:23, 26:29)] <- 1 #formation.hazard.agegapry.numrel_man: ppn growth, man prev, agediff
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[12, x.offset + c(1:4, 14:16, 20:25, 28:29)] <- 1 #formation.hazard.agegapry.numrel_woman: ppn growth, wom prev, agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 17:23, 26:29)] <- 1 #formation.hazard.agegapry.gap_factor_man_exp: ppn growth, prev, agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 14:16, 20:25, 28:29)] <- 1 #formation.hazard.agegapry.gap_factor_woman_exp:  ppn growth, prev, agediff
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[15, x.offset + c(1:4, 26:29)] <- 1 #person.agegap.man.dist.normal.mu: ppn growth, agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(24:25, 28:29)] <- 1 #person.agegap.woman.dist.normal.mu: age diff
  #index.pM <- index.pM + 1
  #mice.predictorMatrix[17, x.offset + c(26:29)] <- 1 #person.agegap.man.dist.normal.sigma: age diff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(1:4, 24:25, 28:29)] <- 1 #person.agegap.woman.dist.normal.sigma: ppn growth, agediff
  index.pM <- index.pM + 1
  mice.predictorMatrix[index.pM, x.offset + c(5:11, 13:29)] <- 1 #hivtransmission.param.f1: ART init, ret, vlsup, AIDS mort, prev, agediff

  #selecting data for calibration
  filter.names <- c(target.variables, input.parameter.names)

  ### target parameters with NA in simpact varied parameters
  target.parameter.na <- c(target.values, rep(NA, length(input.parameter.names)))

  if(nrow(simpact.sim.df) > 5){

    #get only the df with the Target values and Simpact varied parameters
    simpact.sim.df <- simpact.sim.df[ , (names(simpact.sim.df) %in% filter.names)]

    ########      Now create the df for mice imputation #########################
    complete.df.wna <- rbind(simpact.sim.df, target.parameter.na)

    #selected methods used
    cal.list <- c("pmm", "norm", "rf")

    #calibration
    mice.imputted.df <- mice(complete.df.wna, m=out.sets, maxit = iter.num,
                             predictorMatrix = mice.predictorMatrix,
                             method = cal.list[method.id], seed = seed.id)

    #collect all the imputed sets
    mice.out.imputted.df <- tail(mice::complete(mice.imputted.df,1),1)

    for (i in 2:out.sets){
      complete.mice.wo.na <- mice::complete(mice.imputted.df,i)
      mice.out.imputted.df <- rbind(mice.out.imputted.df, tail(complete.mice.wo.na,1))
    }

    #remove all the rows with out of range imputed files
    #mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.man.dist.normal.sigma > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, hivtransmission.param.f1 > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.woman.dist.normal.sigma > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, person.art.accept.threshold.dist.fixed.value > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, monitoring.fraction.log_viralload > 0)

    mice.out.imputted.df$method <- cal.list[method.id]
    mice.out.imputted.df <- cbind.data.frame(sim.id = 1:nrow(mice.out.imputted.df), mice.out.imputted.df)

    #save the data if needed.
    save.string <- paste0(sample(c(LETTERS,letters), 5), collapse = "")
    write.csv(mice.out.imputted.df, file = paste0(dirname,"/",save.name, cal.list[method.id],".",save.string,
                                              ".",training.df.file), row.names = FALSE)
  }else{

    #If not enough simpact returned dataset
    mice.out.imputted.df <- as.data.frame(matrix(NA, 0, length(target.parameter.na)))
    names(mice.out.imputted.df) <- filter.names

  }

  return(mice.out.imputted.df)

}

