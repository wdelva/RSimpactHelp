#' Create a vector of HIV prevalence measurements.
#'
#' Calculate the HIV prevalence, for specific age groups, in yearly intervals,
#' from zero until the end of the simulation.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper
#'   bound) that should be retained, e.g. c(15, 30)
#' @return A dataframe with HIV prevalence estimates and surrounding confidence
#'   bounds, for the specific age group.
#' @examples
#' data(datalist)
#' prev.vector.df <- prevalence.vector.creator(datalist = datalist, agegroup = c(15,50))
#' prev.vector.df
#'
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export

prevalence.vector.creator <- function(datalist = datalist,
                                      agegroup = c(15, 50)){
  prev.agegroup <- agegroup
  timevect <- 0:datalist$itable$population.simtime[1]
  prev.vector.df <- data.frame(timevect = timevect, prev = NA,
                            ll = NA, ul = NA)
  rowindex <- 1
  for (timepoint in timevect) {
    prev.tibble <- prevalence.calculator(datalist = datalist,
                                         agegroup = prev.agegroup, timepoint = timepoint)
    prev.vector.df$prev[rowindex] <- as.numeric(prev.tibble$pointprevalence[3])
    prev.vector.df$ll[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ll[3])
    prev.vector.df$ul[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ul[3])
    rowindex <- rowindex + 1
  }
  prev.vector.df <- prev.vector.df %>% tidyr::gather(type, prevalence,
                                               prev:ul)
  prev.vector.df$pointestimate <- prev.vector.df$type == "prev"
  return(prev.vector.df)
}

