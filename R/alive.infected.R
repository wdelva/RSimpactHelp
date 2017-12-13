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
  df <- datalist$ptable 
  
  if (site == "All") {
    
    df.alive <- df %>%
      dplyr::filter(TOB <= timepoint,
                    TOD > timepoint)
    
  } else{
    
    df.alive <- df %>%
      dplyr::filter(TOB <= timepoint,
                    TOD > timepoint,
                    pfacility == site)
  }
  
  # Now we allocate infection status to all alive people in our table
  df.alive <- df.alive %>%
    dplyr::mutate(Infected = (timepoint >= InfectTime))
  
  return(df.alive)
}
