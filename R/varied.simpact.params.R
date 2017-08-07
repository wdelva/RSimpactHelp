#' Select the parameters that need to be calibrated to meet the target statistics.
#'
#' The default parameters are all the parameters are used during the set
#' up of \code{\link{simpact.config.inputs}}.
#' To exclude these and call for new ones, you will need to call the
#' function with params.replace=TRUE
#'
#' @param params.replace is to return all the parameters used during
#' @param ... standard varied simpact parameters
#'
#' the \code{\link{simpact.config.inputs}} setup.
#' @return a list of parameters to be calibrated
#' @examples input.varied.params.plus <- varied.simpact.params()[[1]]
#' @export

varied.simpact.params <- function(..., params.replace = FALSE){

  if(params.replace == FALSE){

    arg.list <- list(c("person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b",
                       "conception.alpha_base", "formation.hazard.agegapry.numrel_man",
                       "formation.hazard.agegapry.eagerness_diff",
                       "formation.hazard.agegapry.gap_factor_man_exp",
                       "person.agegap.man.dist.normal.mu", "person.agegap.woman.dist.normal.mu",
                       "person.agegap.man.dist.normal.sigma",
                       "person.agegap.woman.dist.normal.sigma",
                       "person.eagerness.woman.dist.gamma.a",
                       "person.eagerness.wooman.dist.gamma.b",
                       "formation.hazard.agegapry.numrel_woman",
                       "formation.hazard.agegapry.gap_factor_woman_exp", ...))

  }else{
    arg.list <- list(c(...))
  }
  return(arg.list)
}



