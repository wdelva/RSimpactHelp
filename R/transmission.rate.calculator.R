#' Calculate the transmission rate
#'
#'Calculate the transmission rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#'
#' data(datalist)
#' transm.rate <- transmission.rate.calculator(datalist = datalist, timewindow = c(20, 40))
#' transm.rate
#'
#' @import dplyr

transmission.rate.calculator <- function(datalist = datalist,
                              timewindow = c(20, 40)){

  Infec.pop.table <- datalist$ptable[InfectType==1]

  numb.infec.pop <- nrow(Infec.pop.table %>%
                           subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))


  transm.rate <- numb.infec.pop / diff(timewindow)

  return(transm.rate)
}
