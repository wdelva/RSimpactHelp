#' A function that returns the total number of people between two points in simulation
#' time for a particular age group whose VL is within a threshold
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param agegroup alive people within this agegroup.
#' @param viralload log10 set point viral load threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL at ART initiation are within a given threshold
#' @examples
#' vl.atARTinit <- vl.atARTinit(datalist = datalist, agegroup = c(15, 30), timewindow = c(15, 30), viralload = c(3, 4), site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr

vl.atARTinit <- function(datalist = datalist, agegroup = c(15, 30),
                             timewindow = c(15, 30), viralload = c(3, 4), site="All"){

  vl.atARTinit.Count <- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  vl.atARTinit.Count <- subset(vl.atARTinit.Count, TreatTime !=Inf) #HIV positive individuals

  vl.atARTinit.Count <- vl.atARTinit.Count %>% mutate(VLThresholdatARTinit = (viralload[1] <= log10SPVL &  log10SPVL <= viralload[2]))

  # log10SPVL is not affected by treatment
  #provide a summary of those that are on treatmen and those that started below a threshold
  vl.atARTinit.Count <- data.frame(dplyr::summarise(dplyr::group_by(vl.atARTinit.Count, Gender),
                                                    TotalCases = n(),
                                                    VLatARTinitThreshold =sum(VLThresholdatARTinit)))


  return(vl.atARTinit.Count)
}

