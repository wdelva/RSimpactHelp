#' Simpact config input parameters to be estimated.
#'
#' This is the first step in setting up for the multivator/emulator parameter estimation.
#' The user needs to list all the parameters that needs to be estimated and also specify
#' the boundaries
#' for these parameters
#' The dataframe that is returned includes the design matrix ([0,1] random numbers for
#' each parameter) and the specific parameter value that will be used to run simpact
#'
#' @param datalist the selected parameters that produced minimal sum of squares
#'
#' @return dataframe with varied parameter values to be used with simpact.run
#'
#@examples
#input.varied.params <- simpact.config.inputs.par.select(datalist)
#' @export

simpact.config.inputs.par.select <- function(datalist = datalist){

  input.parameters <- length(datalist)

  if(length(input.parameters)==0){stop("No parameters have been specified.")}

  #remove columns that are not essential for the simulation
  simpactInputParams <- subset(datalist, select=-c(match, sum.square.df))
  #check when summary stats end
  sim.id.index <- which(names(simpactInputParams)=="sim.id")

  simpactInputParams <- simpactInputParams[,sim.id.index:length(simpactInputParams)]

  return(simpactInputParams)

}
