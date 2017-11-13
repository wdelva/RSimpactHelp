#' A function that returns the percentage of clients who are alive and are
#' virally suppressed 6 or more months
#' after ART initiation, at a paticular point in the simulation time
#'
#' TO REVISE (Need to get the first time that ART was started then reset when
#' drop out and look at the VL six month after, Ignore everthing else
#' after this)
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timepoint alive people within this simulation time e.g timewindow = 30.
#' @param vlcutoff viral load below this threshold e.g <1000 copies/mL is defined as suppressed
#' @param lessmonths months from which time point viral load suppression is sort e.g 6 or
#' more months after ART initiation
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number and percentage of people who are virally suppressed 6 or more
#' months after ART initiation at one timepoint
#' @examples
#' data(datalist)
#' vl.suppressed <- vl.suppressed(datalist = datalist, timepoint=40,
#' vlcutoff=1000, lessmonths = 6, site="All")
#' vl.suppressed
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @export

vl.suppressed <- function(datalist = datalist,
                          timepoint = 35, vlcutoff = 1000,
                          lessmonths = 6, site="All"){

  #Get the individual who started treatment before the time timepoint
  if (site == "All") {
    DT.infected <- subset(datalist$ptable,
                          TreatTime < timepoint)
  } else{

    DT.infected <- subset(datalist$ptable,
                          TreatTime < timepoint,
                          pfacility == site)
  }

  #log10 vlcutoff
  vl.cutoff <- log10(vlcutoff)

  #cut off in months
  .months <- timepoint - lessmonths/12

  #look at the treatment table
  #get the treatment started before the timepoint
  #was on treatment within the time window to check VL

  last.treatment <- datalist$ttable %>%
    subset(ID %in% DT.infected$ID) %>%
    group_by(ID) %>%
    dplyr::filter(row_number()==n()) %>%
    as.data.frame()

  #Get the VL on value due to starting treatment
  vl.table <-  datalist$vltable %>%
    group_by(ID) %>%
    dplyr::filter(Desc =="Started ART",
                  row_number()==n()) %>%
    as.data.frame()

  last.treatment.vl <- dplyr::left_join(x = last.treatment,
                                        y = vl.table,
                                        by = c("ID") )

  #Assuming that we know when Tend.
  #As long as Treatment was started before timepoint and AFTER "lessmonths"
  #clients still on still treatment, we keep the record

  last.treatment.vl <- last.treatment.vl %>%
    mutate(on.treatment.lessmonths = (.months < TEnd),
           vl.suppressed = (Log10VL < vl.cutoff),
           vl.treat.suppressed = on.treatment.lessmonths * vl.suppressed
    )


  if(nrow(last.treatment.vl)==0){

    vlSuppressed.TP <- as.data.frame(matrix(NA, 3, 4))

  }else{
    vlSuppressed.TP <- last.treatment.vl %>%
      group_by(Gender) %>%
      summarise(TotalCases = n(),
      VLSuppressed = sum(vl.treat.suppressed, na.rm = TRUE),
      Percentage = sum(vl.treat.suppressed, na.rm = TRUE)/n() *100) %>%
      as.data.frame()


     vlSuppressed.TP <- rbind(vlSuppressed.TP,
                  c(NA, nrow(last.treatment.vl),
                  sum(last.treatment.vl$vl.treat.suppressed, na.rm = TRUE),
                  sum(last.treatment.vl$vl.treat.suppressed, na.rm = TRUE)/nrow(last.treatment.vl)*100
                  ))

  }

  names(vlSuppressed.TP) <- c("Gender", "TotalCases", "VLSuppressed", "Percentage")


  return(vlSuppressed.TP)
}

