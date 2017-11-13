#This function will calibrate from that simulated dataset
mice.calibration <- function(simpact.sim.df, method.id = 1, cores_2_use = 4,
                             out.sets = 20, iter.num = 15, seed.id = 2525,
                             printflag = FALSE){

  #Add the necessary functions
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/hho.rsimpacthelper.R")
  #Just to get the varied parameter names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.parameters.R")
  inPUT.df.complete <- head(source.simpact.parameters(),1)
  input.parameter.names <- names(dplyr::select(inPUT.df.complete, contains(".")))
  input.parameter.names <- input.parameter.names[2:length(input.parameter.names)]
  par.imputed <- length(input.parameter.names)

  #get the target values and their names
  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")

  #selection data for calibration
  filter.names <- c(target.variables, input.parameter.names)

  ### target parameters with NA in simpact varied parameters
  target.parameter.na <- c(target.values, rep(NA, length(input.parameter.names)))

  if(nrow(simpact.sim.df) > 10){

    #get only the df with the Target values and Simpact varied parameters
    simpact.sim.df <- simpact.sim.df[ , (names(simpact.sim.df) %in% filter.names)]

    ######## Now create the df for mice imputation #########################
    complete.df.wna <- rbind(simpact.sim.df, target.parameter.na)

    #selected methods used
    cal.list <- c("pmm", "norm", "norm.boot")

    out.sets <- out.sets%/%cores_2_use + out.sets%%cores_2_use

    cl <- makeCluster(cores_2_use)
    clusterSetRNGStream(cl, seed.id)
    registerDoParallel(cl)

    imp_merged <-
      foreach(no = 1:cores_2_use,
              .combine = ibind,
              .packages = "mice") %dopar%
              {
                mice(complete.df.wna, m = out.sets, printFlag = printflag,
                     maxit = iter.num, method = cal.list[method.id])
              }
    stopCluster(cl)

    mice.imputted.df <- unlist(imp_merged$imp) %>%
      matrix(., byrow = FALSE, ncol = par.imputed) %>%
      as.data.frame()

    #give the calibrated values names
    names(mice.imputted.df) <- input.parameter.names

    mice.out.imputted.df <- subset(mice.imputted.df, person.agegap.man.dist.normal.sigma > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.woman.dist.normal.sigma > 0)
    mice.out.imputted.df <- subset(mice.out.imputted.df, hivtransmission.param.f1 > 0)

    #save the data if needed.
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    save.string <- paste0(sample(c(LETTERS,letters), 5), collapse = "")
    write.csv(mice.out.imputted.df, file = paste0(dirname,"mice.",cal.list[method.id],".",save.string,
                                              ".",training.df), row.names = FALSE)
  }else{

    #If not enough simpact returned dataset
    mice.out.imputted.df <- as.data.frame(matrix(NA, 0, length(target.parameter.na)))
    names(mice.out.imputted.df) <- filter.names

  }

  return(mice.out.imputted.df)

}

