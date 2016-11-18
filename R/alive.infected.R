#' Subset the datalist ptable to include only people who were alive at a point in time and add their HIV status at that point in time.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timepoint Point in time at which the subset should be created and HIV status should be evaluated.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a data.table that only includes people who were alive at the timepoint and that records their HIV status.
#' @examples
#' alive.twenty.dt <- alive.infected(DT = datalist, timepoint = 20)

alive.infected <- function(datalist = datalist,
                           timepoint = 40, site="All"){ # arguments are the personlog data.table and a point in time
  DT <- datalist$ptable
  if(site=="All"){
    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  }else{
    facilities.df <- read.cv(datalist$itable$facilities.geo.coords)
    facilities.df <- filter(facilities.df, Facility = site[1])
    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint
                      & XCoord==facilities.df$Longitude
                      & YCoord==facilities.df$Latitude)
  }

  #DTalive$Infected <- timepoint >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people

  DTalive <- DTalive %>%  mutate(Infected = (timepoint >= InfectTime))
  return(DTalive)
}
