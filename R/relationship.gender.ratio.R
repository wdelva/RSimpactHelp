#' Relationship gender ratio
#'
#' Calculate the relationship gender ratio (ratio of men & women in relationship) over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 40)
#' data(datalist)
#' rels.ratio <- relationship.gender.ratio(datalist = datalist, timewindow = c(20, 40))
#' rels.ratio
#'
#' @import dplyr
#'
# Relationship ratio for men and women

relationship.gender.ratio <- function(datalist = datalist,
                              timewindow = c(20, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.ratio <- vector("list", length(c(0,1)))

  for(i in 1:2){
    rels.ratio$men <- (numb.rels.men/numb.rels)  #/diff(timewindow)
    rels.ratio$women <- (numb.rels.women/numb.rels)  #/diff(timewindow)
  }
  return(rels.ratio)
}

