#' Simpact config input parameters to be estimated.
#'
#' This is the first step in setting up for the multivator/emulator parameter estimation.
#' The user needs to list all the parameters that needs to be estimated and also specify
#' the boundaries for these parameters
#' The dataframe that is returned includes the design matrix ([0,1] random numbers for
#' each parameter) and the specific parameter value that will be used to run simpact
#'
#' @param design.points the number simulation rows that needs to be perfomed by the simpact.run
#' @param resample.count replace the samples with these number of resample
#' @param ... standard varied simpact parameters
#' @return dataframe with varied parameter values to be used with simpact.run
#'
#' @examples
#' input.varied.params <- simpact.config.inputs(conception.alpha_base = c(-3.6, -1.2),
#' formation.hazard.agegapry.baseline = c(1.5, 3), person.eagerness.man.dist.gamma.a = c(0.1, 2),
#' person.eagerness.man.dist.gamma.b = c(5, 60),
#' formation.hazard.agegapry.numrel_man = c(-1.5, -0.1),
#' formation.hazard.agegapry.eagerness_diff = c(-0.1, 0),
#' formation.hazard.agegapry.gap_factor_man_exp = c(-1.5, -0.4),
#' person.agegap.man.dist.normal.mu = c(0, 4), person.agegap.woman.dist.normal.mu = c(0, 4),
#' person.agegap.man.dist.normal.sigma = c(0.5,2),person.agegap.woman.dist.normal.sigma =c(0.5,2))
#'
#' @importFrom lhs randomLHS
#' @importFrom lhs augmentLHS
#' @importFrom lhs improvedLHS
#' @importFrom stats qunif
#' @export

simpact.config.inputs <- function(design.points = 10, resample.count = 1, ...){

  input.parameters <- list(...)

  if(length(input.parameters)==0){stop("No parameters have been specified.
                                       Use e.g input.varied.params
                                       (conception.alpha_base = c(-3.6, -1.2))")}

  # Creating the LHS over the 0-1 uniform parameter space for the parameters to be estimated
  variables <- length(input.parameters)
  set.seed(1976)

  if(design.points > 4){
    add.design.points <- design.points %/% 4
    aug.points <- design.points %% 4

    rlhs <- improvedLHS(add.design.points, variables)

    for(i in 1:3){

      rlhs <- augmentLHS(rlhs, add.design.points)

    }
    if(aug.points > 0){

      rlhs <- augmentLHS(rlhs, aug.points)
    }

  }else{
    rlhs <- randomLHS(design.points, variables)
  }

  if (resample.count > 1){
    for(i in 1:resample.count){

      # Re-sample from the same grid as initially selected
      rlhs <- augmentLHS(rlhs, design.points)
    }
    rlhs <- tail(rlhs, design.points)
  }

  #Select the config parameters that will be varied from the input config
  lhs.df <- data.frame(matrix(NA, nrow = 1, ncol = variables))
  names(lhs.df) <- names(input.parameters)

  #Repeat have to meet the design.points
  lhs.df <- as.data.frame(lapply(lhs.df, rep, design.points))

  #Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
  x.index <- 0
  for (j in names(input.parameters)){
    x.index <- x.index + 1
    min.var <- input.parameters[j][[1]][1]
    max.var <- input.parameters[j][[1]][2]
    col.index <- which(names(lhs.df)==j)
    lhs.df[col.index] <- qunif(rlhs[, x.index], min = as.numeric(min.var), max =as.numeric(max.var))
  }

  #Name the xdesign dataframe
  rlhs <- data.frame(rlhs)
  names(rlhs) <- paste0("xdesign",1:length(rlhs))

  ##This will create the input file for the simpact
  simpactInputParams <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df)


  return(simpactInputParams)

}


