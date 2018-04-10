#' Create a vector of ART coverage measurements.
#'
#' Calculate the ART coverage, for a specific age group, in yearly intervals,
#' from the start of the simulation until the end of the simulation.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper
#'   bound) that should be retained, e.g. c(15, 30)
#' @return A vector with ART coverage estimates, for the specific age group.
#' @examples
#' data(datalist)
#' art.cov.vector.df <- ART.coverage.vector.creator(datalist = datalist, agegroup = c(15,50))
#' art.cov.vector.df
#'
#' @export

ART.coverage.vector.creator <- function(datalist = datalist,
                                     agegroup = c(15, 50)){
  ART.cov.eval.timepoints <- seq(from = 0,
                                 to = datalist$itable$population.simtime[1])
  ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  for (art.cov.index in 1:length(ART.cov.vector)){
    ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist,
                                                                 agegroup = agegroup,
                                                                 timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
      sum(ART.coverage.calculator(datalist = datalist,
                                  agegroup = agegroup,
                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  }
  return(ART.cov.vector)
}
