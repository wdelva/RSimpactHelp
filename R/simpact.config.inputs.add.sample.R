#' Resampling of Simpact config input parameters to be estimated.
#'
#' This is the first step in setting up for the multivator/emulator parameter estimation.
#' The user needs to list all the parameters that needs to be estimated and also specify
#' the boundaries for these parameters
#' The dataframe that is returned includes the design matrix ([0,1] random numbers for
#' each parameter) and the specific parameter value that will be used to run simpact
#' as used in the initial set up
#'
#' @param datalist the dataframe produced by simpact.wrapper.easyABC.run
#' @param resample.points Number of points to be added standard varied simpact parameters
#' @param set.seed.new In case cluster, you need to be using different seed for resampling
#' @param ... a set of parameter boundaries to be used for new re-sampling
#' @return dataframe with resampled parameter values to be used with simpact.run
#'
#@examples
#input.varied.params <- simpact.config.inputs.add.sample(complete.summary.output,
#formation.hazard.agegapry.baseline = c(1.5, 3), person.eagerness.man.dist.gamma.a = c(0.1, 2),
#'
#' @import dplyr
#' @importFrom lhs augmentLHS
#' @importFrom stats qunif
#' @export

simpact.config.inputs.add.sample <- function(datalist = datalist,
                                             resample.points = 10,
                                             set.seed.new = 1, ...){

  #Also use as sim.id.count to reflect the continuous sim.id count
  input.data <- nrow(datalist)

  #parameter boundaries
  input.boundaries <- list(...)

  if(length(input.data)==0){stop("No parameters have been specified.")}

  #continuous count on the sim.id range
  sim.id.start <- input.data + 1
  sim.id.end <- input.data + resample.points

  set.seed(set.seed.new)
  #select the previous set random numbers
  rlhs <- as.matrix(dplyr::select(datalist, contains("xdesign")))

  # Re-sample from the same grid as initially selected
  rlhs <- augmentLHS(rlhs, resample.points)
  rlhs <- tail(rlhs, resample.points)

  #find columncount of max xdesign
  max.xdesign.pos <- which(colnames(datalist)==paste0("xdesign",dim(rlhs)[2])) + 1
  input.parameters <- names(datalist[,max.xdesign.pos:length(datalist)])
  variables <- length(input.parameters)

  bound.variabled <- length(input.boundaries)

  if(variables != bound.variabled){stop("Input parameters their boundaries do not match.")}

  #Select the config parameters that will be varied from the input config
  lhs.df <- data.frame(matrix(NA, nrow = 1, ncol = variables))
  names(lhs.df) <- input.parameters

  #Repeat have to meet the design.points
  lhs.df <- as.data.frame(lapply(lhs.df, rep, resample.points))

  #Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
  x.index <- 0
  for (j in input.parameters){
    x.index <- x.index + 1
    min.var <- input.boundaries[j][[1]][1]
    max.var <- input.boundaries[j][[1]][2]
    col.index <- which(names(lhs.df)==j)
    lhs.df[col.index] <- qunif(rlhs[, x.index], min = as.numeric(min.var), max =as.numeric(max.var))
  }

  #Name the xdesign dataframe
  rlhs <- data.frame(rlhs)
  names(rlhs) <- paste0("xdesign",1:length(rlhs))

  ##This will create the input file for the simpact
  simpactInputParams <- data.frame(cbind.data.frame(sim.id = sim.id.start:sim.id.end, rlhs, lhs.df))


  return(simpactInputParams)

}


