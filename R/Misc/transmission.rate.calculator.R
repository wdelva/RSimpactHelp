#' Calculate the transmission rate
#'
#'Calculate the transmission rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#' @param int Logic parameter which is FALSE by deafault meaning that only one value for the transmission rate
#' will be computed, it is TRUE when the time window is being split in small time interval and get different
#' transmission rates within these time intervals
#' @param by Time step unit to subdivide the time window
#' data(datalist)
#' transm.rate <- transmission.rate.calculator(datalist = datalist,
#' timewindow = c(20, 40), int = FALSE, by=1)
#' transm.rate
#'
#' @import dplyr
#' @export

transmission.rate.calculator <- function(datalist = datalist,
                              timewindow = c(20, 40), int = FALSE, by=1){

  if(int==FALSE){

    Infec.pop.table <- datalist$ptable[InfectType==1]

    numb.infec.pop <- nrow(Infec.pop.table %>%
                             subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))


    transm.rate <- numb.infec.pop / diff(timewindow)

    return(transm.rate)
  }

  if(int==TRUE){

    Infec.pop.table <- datalist$ptable[InfectType==1]

    upper.limit <- ceiling(diff(timewindow)/by)

    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

    infec.pop.int <- vector()
    trans.rate.int <- vector()


    for(i in 0:(upper.limit-2)){

      timewindow.int <- c(interv.time[1+i], interv.time[2+i])

      infec.pop.numb <- nrow(Infec.pop.table %>%
                               subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))


      infec.pop.int <- c(infec.pop.int,  infec.pop.numb)


      trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))

    }

    return(trans.rate.int)
  }

}
