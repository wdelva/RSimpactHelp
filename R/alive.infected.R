#' Subset the datalist ptable to people alive at a point in time and add HIV status.
#'
#' Subset the datalist ptable to include only people who were alive at a point in time and add their HIV status at that point in time.
#'
#' @param DT The person table (ptable) that is produced by \code{\link{readthedata()}}
#' @param timepoint Point in time at which the subset should be created and HIV status should be evaluated.
#' @return a data.table that only includes people who were alive at the timepoint and that records their HIV status.
#' @examples
#' alive.twenty.dt <- alive.infected(DT = datalist$ptable, timepoint = 20)

alive.infected <- function(DT = datalist$ptable,
                           timepoint = 20){ # arguments are the personlog data.table and a point in time
  DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  DTalive$Infected <- timepoint >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people
  return(DTalive)
}
