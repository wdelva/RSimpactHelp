#' A function that returns the total number of people between two points in simulation
#' time for a particular age group whose VL is within a threshold
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup alive people within this agegroup.
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param viralload log10 set point viral load threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL at ART initiation are within a given threshold
#' @examples
#' data(datalist)
#' vl.atARTinit <- vl.atARTinit(datalist = datalist, agegroup = c(15, 30),
#' timewindow = c(15, 30), viralload = c(3, 4))
#'
#' @importFrom  magrittr %>%
#' @import dplyr
#' @export

vl.atARTinit <- function(datalist = datalist, agegroup = c(15, 30),
                         timewindow = c(15, 30),
                         viralload = c(3, 4), site="All"){

  vl.atARTinit.Count <- age.group.time.window(datalist = datalist,
                                              agegroup = agegroup,
                                              timewindow = timewindow, site="All")

  #HIV positive individuals
  vl.atARTinit.Count <- subset(vl.atARTinit.Count, TreatTime !=Inf)

  #VL within a selected threshhold
  vl.at.treatment.df <- datalist$vltable %>%
    subset(ID %in% vl.atARTinit.Count$ID) %>%
    filter(Desc == "Started ART", Time >= timewindow[1],
           Time <= timewindow[2]) %>%
    group_by(ID) %>%
    filter(row_number() == n()) %>%
    data.frame()

  vl.at.treatment.df <- vl.at.treatment.df %>%
    dplyr::select(ID, Log10SPVL) %>%
    dplyr::mutate(vl.atARTinit = (viralload[1] <= Log10SPVL &
                                    Log10SPVL <= viralload[2]))

  #provide a summary of those that are on treatment and those that started within a threshold
  vl.atARTinit.Count <- vl.at.treatment.df %>%
    group_by(Gender) %>%
    summarise(
      TotalCases = n(),
      VLatARTinitThreshold =sum(VLThresholdatARTinit, na.rm = TRUE)) %>%
    as.data.frame()

  return(vl.atARTinit.Count)
}

