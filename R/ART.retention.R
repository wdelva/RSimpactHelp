#' A function that returns the percentage of clients who are on treatment
#' at 6 or more months after ART initiation,
#' at a paticular point in the simulation time
#'
#' Definition. Percentage of adults and children with HIV alive and on ART 12, 24, 36 (etc.)
#' months after initiating treatment among patients initiating ART during a specified time period.
#'
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param ARTtimewindow the window which treatment was started e.g timepoint = 30.
#' @param retentiontimeMonths months from which retention is to be computed e.g 6 or
#' more months after ART initiation
#' @param site Select only the particular site from the study, if all use all sites.
#' @return fraction of people who started treatment and are still on treatment after stated time
#' @examples
#' data(datalist)
#' ART.retention <- ART.retention(datalist = datalist, ARTtimewindow= c(30,31),
#' retentiontimeMonths = 6, site="All")
#' ART.retention
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @export

ART.retention <- function(datalist = datalist,
                          agegroup = c(15,150),
                          ARTtimewindow = c(20,34),
                          retentiontimeMonths = 12, #6 months default
                          site="All"){

  #get the treatment table
  DT.treatment <- datalist$ttable

  #Get people who have started treatment within timepoint
  if (site == "All") {
    DT.treatment <- subset(DT.treatment,
                           (TStart >= ARTtimewindow[1] &
                           TStart < ARTtimewindow[2])
                           )
  } else{
    DT.treatment <- subset(DT.treatment,
                           (TStart >= ARTtimewindow[1] &
                              TStart < ARTtimewindow[2]) &
                               pfacility == site )
  }

  #Get the right age within the time window
  age.group.ART <- age.group.time.window(datalist = datalist,
                            agegroup = agegroup,
                            timewindow = ARTtimewindow, site = site)


  if(nrow(DT.treatment)>0){
    #get everyone last known treatment status
    DT.treatment <- DT.treatment %>%
      group_by(ID) %>%
      dplyr::filter(row_number()==n()) %>%
      as.data.frame

    DT.treatment$retention.time <- DT.treatment$TStart + retentiontimeMonths/12

    #restrict the list to the right agegroup
    DT.treatment <- subset(DT.treatment, ID %in% age.group.ART$ID)



    DT.treatment$retention <- DT.treatment$TEnd > DT.treatment$retention.time

    DT.treatment.retention <- DT.treatment %>%
      group_by(Gender) %>%
      summarise(TotalCases = n(),
                ART.retention = sum(retention),
                percentage = ART.retention / TotalCases * 100 ) %>%
      as.data.frame

  }else{ DT.treatment.retention <- as.data.frame(matrix(NA, 0, 4))}


  if(nrow(DT.treatment.retention)==0){
    DT.treatment.retention <- as.data.frame(matrix(NA, 3, 4))
  }else{
  DT.treatment.retention <- rbind(DT.treatment.retention,
                                  c(NA, nrow(DT.treatment), sum(DT.treatment$retention),
                                    sum(DT.treatment$retention)/nrow(DT.treatment)*100) )
  }

  names(DT.treatment.retention) <- c("Gender","TotalCases","ART.retention", "percentage")

  return(DT.treatment.retention)
}

