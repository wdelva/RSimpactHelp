#' A function that returns the total number of people between two points in simulation
#' time for a particular age group and gender whose vl was between the given threshold
#' and CD4 count was between the given threshold
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param agegroup alive people within this agegroup.
#' @param viralload at the time of ART initiation.
#' @param cd4count threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL and CD4 count at ART initiation are given
#' @examples
#' data(datalist)
#' cd4.vl.atARTinit <- cd4.vl.atARTinit(datalist = datalist,
#' agegroup=c(15,40), timewindow=c(15,40), viralload=c(3,4), cd4count=c(350,500), site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

cd4.vl.atARTinit <- function(datalist = datalist, agegroup = c(15,40),
                             timewindow = c(15,40), viralload = c(3,4),
                             cd4count = c(350,500), site="All"){

  #HIV positive individuals
  cd4.vl.atARTinit <- subset(cd4.vl.atARTinit, TreatTime !=Inf)

  raw.df <- data.frame(cd4.vl.atARTinit)
  art.df <- subset(datalist$ttable, ID %in% cd4.vl.atARTinit$ID &
                     TStart >= timewindow[1] & TStart <= timewindow[2])

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #select those who started their treatment when their CD4 count was between the given threshold
  raw.df <- raw.df %>%
    dplyr::mutate(cd4count.atARTInit = (cd4count[1] <= CD4atARTstart &
                                          CD4atARTstart <= cd4count[2]))

  #select those who started their treatment when their vl was between the given threshold
  vl.at.treatment.df <- datalist$vltable %>%
    subset(ID %in% raw.df$ID) %>%
    filter(Desc == "Started ART", Time >= timewindow[1],
                  Time <= timewindow[2]) %>%
    group_by(ID) %>%
    filter(row_number() == n()) %>%
    data.frame()

  vl.at.treatment.df <- vl.at.treatment.df %>%
    dplyr::select(ID, Log10SPVL) %>%
    dplyr::mutate(vl.atARTinit = (viralload[1] <= Log10SPVL &
                                    Log10SPVL <= viralload[2]))

  # Now we apply the left_join to add the Log10SPVl at ART init status to raw.df.
  raw.df <- left_join(raw.df, vl.at.treatment.df, by = "ID")

  #provide a summary of those that are on treatment and those that started below a threshold
  cd4.vl.atARTinit <- raw.df %>%
    group_by(Gender) %>%
    summarise(
              TotalCases = n(),
              cd4count.atARTInit = sum(cd4count.atARTInit, na.rm = TRUE),
              vl.atARTinit = sum(vl.atARTinit, na.rm = TRUE)) %>%
    as.data.frame()

  cd4.vl.atARTinit.all <- raw.df %>%
    summarise(
      Gender = NA,
      TotalCases = n(),
      cd4count.atARTInit = sum(cd4count.atARTInit, na.rm = TRUE),
      vl.atARTinit = sum(vl.atARTinit, na.rm = TRUE)) %>%
    as.data.frame()

  cd4.vl.atARTinit <- rbind(cd4.vl.atARTinit, cd4.vl.atARTinit.all)

  return(cd4.vl.atARTinit)
}

