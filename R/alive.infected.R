#' Subset the datalist ptable to include only people who were alive at a point
#'  in time and add their HIV status at that point in time.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timepoint Point in time at which the subset should be created and HIV status
#' should be evaluated.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a data.table that only includes people who were alive at the timepoint and that
#' records their HIV status.
#' @examples
#' data(datalist)
#' alive.at.twenty.dt <- alive.infected(datalist = datalist, timepoint = 20,
#' site = "All")
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

alive.infected <- function(datalist = datalist,
                           timepoint = 40,
                           site = "All") {
  # arguments are the personlog data.table and a point in time
  DT <- datalist$ptable
  DT$pfacility <- "NA"
  pf.index <- which(colnames(DT)=="pfacility") #person.facility.index

  if (site == "All") {
    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  } else{
    facilities.df <- datalist.test$ftable
    colnames(facilities.df) <- c("facility.xy", "XCoord", "YCoord")

    for (i in 1:nrow(DT)) {
      fc.id <- which.min(sqrt((DT[i, XCoord] - facilities.df$XCoord)^2 +
                         (DT[i, YCoord] - facilities.df$YCoord)^2 ))

      DT[i, pf.index] <- facilities.df[fc.id, facility.xy]
    }

    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint & pfacility == site)
  }

  # Now we allocate infection status to all alive people in our table
  DTalive <-  DTalive %>%  dplyr::mutate(Infected = (timepoint >= InfectTime))

  return(DTalive)
}
