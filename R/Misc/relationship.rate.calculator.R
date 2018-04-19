#' Relationship rate
#'
#'Calculate the relationship rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#' @param int Logic parameter which is FALSE by deafault and means tha only one value for the relationship rate
#' needed, it is TRUE when the time window is being split in small time interval and get different
#' relationship rates within these time intervals
#' @param by Time step unit to subdivide the time window
#'
#' data(datalist)
#' rels.rate <- relationship.rate.calculator(datalist = datalist,
#' timewindow = c(20, 40), int = FALSE, by=1)
#' rels.rate
#'
#' @import dplyr


relationship.rate.calculator <- function(datalist = datalist,
                      timewindow = c(20, 40), int = FALSE, by=1){

  if(int==FALSE){
    Rels.table <- datalist$rtable

    Rels.table.window <- Rels.table %>%
      subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

    numb.rels <- nrow(Rels.table.window)

    numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

    numb.rels.women <- length(unique(Rels.table.window$ID2))

    rels.rate <- (numb.rels)/ diff(timewindow)

    return(rels.rate)
  }

  if(int==TRUE){

    Rels.table <- datalist$rtable

    upper.limit <- ceiling(diff(timewindow)/by)

    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

    rels.int <- vector()
    rels.rate.int <- vector()


    for(i in 0:(upper.limit-2)){

      timewindow.int <- c(interv.time[1+i], interv.time[2+i])

      rels.numb <- nrow(Rels.table %>%
                          subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))


      rels.int <- c(rels.int,  rels.numb)


      rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))

    }

    return(rels.rate.int)
  }


}
