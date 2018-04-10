#' Create a vector of HIV incidence measurements.
#'
#' Calculate the HIV incidence, for specific age groups, in yearly intervals,
#' from one until the end of the simulation.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper
#'   bound) that should be retained, e.g. c(15, 30)
#' @return A dataframe with HIV incidence estimates, for the specific age group.
#' @examples
#' data(datalist)
#' inc.vector.df <- incidence.vector.creator(datalist = datalist, agegroup = c(15,50))
#' inc.vector.df
#'
#' @export

incidence.vector.creator <- function(datalist = datalist,
                                     agegroup = c(15, 50)){
  inc.vector <- rep(NA, datalist$itable$population.simtime[1])
  for (inc.vector.index in 1:length(inc.vector)){
    inc.vector[inc.vector.index] <- incidence.calculator(datalist = datalist,
                                                         agegroup = agegroup,
                                                         timewindow = c((inc.vector.index - 1),
                                                                        inc.vector.index))$incidence[3]
  }
  return(inc.vector)
}

