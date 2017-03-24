#' Select the summary statistics that you want from the simpact.run.
#'
#' The default parameters are all the parameters produced by output.summary.maker.
#' To have less of these you will need to call the function with summary.replace=TRUE
#'
#' @param ... standard summary statistics
#' @param summary.replace Set to True if you want to use different statistics
#' @return a list of parameters to be included in the summary statistics of the simulation
#' @examples
#' target.variables <- select.summary.params()[[1]]
#' @export

select.summary.params <- function(..., summary.replace = FALSE){

  if(summary.replace == FALSE){

    arg.list <- list(c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25",
                       "prev.men.25.50", "ART.cov.15.50", "incid.wom.15.30",
                       "frac.degreeGT1.wom.15.30", "mean.degree", "median.degree",
                       "Q1.degree", "Q3.degree",...))

  }else{
    arg.list <- list(c(...))
  }
  return(arg.list)
}

