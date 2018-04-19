#' Transmission gender ratio
#'
#' Calculate the transmission gender ratio (ratio of men & women in relationship) over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#' data(datalist)
#' transm.ratio <- transmission.gender.ratio(datalist = datalist, timewindow = c(20, 40))
#' transm.ratio
#'
#' @import dplyr

# Transmission ratio for men and women

transmission.gender.ratio <- function(datalist = datalist,
                                timewindow = c(20, 40)){

  pers.table.hiv <- datalist$ptable[InfectType==1]

  pers.table.hiv.window <- pers.table.hiv %>%
    subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1])

  numb.transm <- nrow(pers.table.hiv.window)

  numb.transm.men <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID))) # Gender 0 men

  numb.transm.women <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID)))# Gender 1 women

  transm.ratio <- vector("list", length(c(0,1)))

  for(i in 1:2){
    transm.ratio$men <- (numb.transm.men/numb.transm)  #/diff(timewindow)
    transm.ratio$women <- (numb.transm.women/numb.transm)  #/diff(timewindow)
  }
  return(transm.ratio)
}
