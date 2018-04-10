#' Create a vector of mean ages among the HIV positive people
#'
#' Calculate the HIV prevalence, for specific age groups, in yearly intervals,
#' from zero until the end of the simulation.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timevector Vector of time point at which the mean age among HIV
#'   positive people must be calculated.
#' @return A vector of mean ages among the HIV positive people.
#' @examples
#' data(datalist)
#' mean.age.pos.vector <- mean.age.pos.vector.creator(datalist = datalist, timevector = 10:60)
#' prev.vector.df
#'
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @export

mean.age.pos.vector.creator <- function(datalist = datalist,
                                        timevector = datalist$itable$hivseed.time[1]:datalist$itable$population.simtime[1]){
  mean.age.pos.vector <- rep(NA, length(timevector))
  for (time.index in 1:length(timevector)) {
    alive.pos <- alive.infected(datalist = datalist,
                                timepoint = timevector[time.index],
                                site = "All") %>% dplyr::filter(Infected == TRUE)
    agevector.pos <- timevector[time.index] - as.numeric(alive.pos$TOB)
    mean.age.pos.vector[time.index] <- mean(agevector.pos, na.rm = TRUE)
  }
  return(mean.age.pos.vector)
}
