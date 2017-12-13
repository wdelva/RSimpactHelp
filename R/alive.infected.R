#' Subset people alive at a point in time and adds their HIV status.
#' 
#' Produces a data frame that contains persons from the simulation that were
#' alive at a specified point in time. The new data frame also contains a
#' variable that indicates whether or not the person was infected with HIV at
#' that point in time.
#'
#' @param datalist The list object that is produced by \code{\link{readthedata}}
#' @param timepoint Point in time at which the subset should be created and HIV 
#'   status should be evaluated.
#' @param site Users can specify a particular site they are interested in. The 
#'   default is "All", which means that persons from all sites should be 
#'   included.
#' @return A data frame that contains a newly generated variable called 
#'   "Infected" that records the HIV status of persons at the time point
#'   indicated.
#' @examples
#' cfg <- list()
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp")
#' dl <- readthedata(modeloutput)
#' alive.at.twenty.df <- alive.infected(datalist = dl, timepoint = 20)
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
