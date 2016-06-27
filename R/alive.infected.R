#' Subset the datalist ptable to people alive at a point in time and add HIV status.
#'
#' Subset the datalist ptable to include only people who were alive at a point in time and add their HIV status at that point in time.
#'
#' @param DT The person table (ptable) that is produced by \code{\link{readthedata()}}
#' @param timepoint Point in time at which the subset should be created and HIV status should be evaluated.
#' @param dec.places Number of decimals for the Age variable that is also added to the dataset.
#' @return a data.table that only includes people who were alive at the timepoint and that records their HIV status and age at that timepoint.
#' @examples
#' alive.twenty.dt <- alive.infected(DT = datalist$ptable, timepoint = 20, dec.places = 0)

alive.infected <- function(DT = datalist$ptable,
                           timepoint = 20,
                           dec.places = 0){ # arguments are the personlog data.table and a point in time
  DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  DTalive$Age <- time -round(ceiling(DTalive$TOB), dec.places) # Next we allocate them in discrete age bins with bin size as wide as timestep
  DTalive$Infected <- time >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people
  return(DTalive)
}
