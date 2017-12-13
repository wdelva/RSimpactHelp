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
  
  # Use personlog of the datalist
  df <- datalist$ptable 
  
  # Check to see if a specific site has been specified
  if (site == "All") {
    
    df.alive <- df %>%
      dplyr::filter(TOB <= timepoint, # Keep those born before timepoint
                    TOD > timepoint) # Keep those died after timepoint
    
  } else{
    
    df.alive <- df %>%
      dplyr::filter(TOB <= timepoint,
                    TOD > timepoint,
                    pfacility == site) # Keep those at the facility specified
  }
  
  # Now we allocate infection status to all alive people in our table
  df.alive <- df.alive %>%
    dplyr::mutate(Infected = (timepoint >= InfectTime)) # TRUE if time point is after infection time
  
  return(df.alive)
}
