#' Relationship rate
#'
#'Calculate the relationship rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#'
#' data(datalist)
#' rels.rate <- relationship.rate(datalist = datalist, timewindow = c(20, 40))
#' rels.rate
#'
#' @import dplyr


relationship.rate <- function(datalist = datalist,
                      timewindow = c(20, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.rate <- vector("list", length(c(0,1,2)))

  for(i in 1:3){
    rels.rate$overall <- (numb.rels)/ diff(timewindow)
    rels.rate$men <- (numb.rels.men)/diff(timewindow)
    rels.rate$women <- (numb.rels.women)/diff(timewindow)
  }
  return(rels.rate)
}
