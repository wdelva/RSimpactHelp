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

  #get everyone last known treatment status
  DT.treatment <- DT.treatment %>%
    group_by(ID) %>%
    filter(row_number()==n()) %>%
    as.data.frame

  DT.treatment$retention.time <- DT.treatment$TStart + retentiontimeMonths/12

  #restrict the list to the right agegroup
  DT.treatment <- subset(DT.treatment, ID %in% age.group.ART$ID)

  if(nrow(DT.treatment)>0){

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

