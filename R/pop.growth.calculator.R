#' Calculate population growth rate.
#'
#' Calculate the growth rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that
#' should be retained, e.g. c(20, 30)
#' @return population growth rate estimate, for the specified time window
#' @examples
#' data(datalist)
#' growth.rate <- pop.growth.calculator(datalist = datalist, timewindow = c(0, 20))
#' growth.rate
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export

pop.growth.calculator <- function(datalist = datalist,
                                  timewindow = c(0, 20)){

  end.popsize <- datalist$ltable %>%
    subset(Time==timewindow[2]) %>%
    dplyr::select(PopSize) %>%
    as.numeric()

  start.popsize <- ifelse(timewindow[1]==0,
                          yes = datalist$itable$population.nummen +
                            datalist$itable$population.numwomen,
                          no = datalist$ltable %>% subset(Time==timewindow[1]) %>%
                            select(PopSize) %>% as.numeric())

  growth.rate <- log(end.popsize/start.popsize) / diff(timewindow)

  return(growth.rate)
}
