#' A function that returns the percentage of clients who are alive and are virally suppressed 6 or more months
#' after ART initiation, at a paticular point in the simulation time
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timepoint alive people within this simulation time e.g timewindow = 30.
#' @param vlcutoff viral load below this threshold e.g <1000 copies/mL is defined as suppressed
#' @param lessmonths months from which time point viral load suppression is sort e.g 6 or more months after ART initiation
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number and percentage of people who are virally suppressed 6 or more months after ART initiation
#' at one timepoint
#' @examples
#' vl.suppressed <- vl.suppressed(datalist = datalist, timepoint=40, vlcutoff=1000, lessmonths = 6, site="All")

vl.suppressed <- function(datalist = datalist, timepoint = 30, vlcutoff = 1000, lessmonths = 6, site="All"){

  DTalive.infected <- datalist$ptable
  DTalive.infected$pfacility <- "NA"
  pf.index <- which(colnames(DTalive.infected)=="pfacility") #person.facility.index

  if (site == "All") {
    DTalive.infected <- subset(DTalive.infected, TOB <= timepoint & TOD > timepoint & InfectTime!=Inf)
  } else{
    facilities.df <- datalist.test$ftable
    colnames(facilities.df) <- c("facility.xy", "XCoord", "YCoord")

    for (i in 1:nrow(DTalive.infected)) {
      fc.id <-
        which.min(sqrt((DTalive.infected[i, XCoord] - facilities.df$XCoord)^2 + (DTalive.infected[i, YCoord] - facilities.df$YCoord)^2
        ))
      DTalive.infected[i, pf.index] <- facilities.df[fc.id, facility.xy]
    }

    DTalive.infected <- subset(DTalive.infected, TOB <= timepoint & TOD > timepoint & pfacility == site & InfectTime!=Inf)
  }

  vl.cutoff <- log10(vlcutoff)
  six.months <- timepoint - lessmonths/12

  DTalive.infected <- DTalive.infected %>%  mutate(ARTLess6mnths = (TreatTime <= six.months))
  raw.df <- data.frame(DTalive.infected)

  InfectedOnTreatment <- subset(raw.df, ARTLess6mnths==TRUE)

  #Check the vl information of these individuals
  VLevent.df <- subset(datalist$vltable, ID %in% InfectedOnTreatment$ID )

  yaxis <- vl.cutoff + 0.2

  #Visualise the VL points
  q <- ggplot()
  q <- q + geom_point(data=VLevent.df, aes(x=VLTimeLog, y=Log10VL, group=ID, colour = Description))
  q <- q + geom_line(data=VLevent.df, aes(x=VLTimeLog, y=Log10VL, group=ID, colour = Description))
  q <- q + geom_hline(yintercept=vl.cutoff) + annotate("text", datalist$itable$hivseed.time[1], vl.cutoff + 0.2, label = "  VLCutoff")
  q <- q + geom_vline(xintercept=timepoint-0.5, colour = "blue") + annotate("text", timepoint, max(VLevent.df$Log10VL), label = "  TimePoint")
  q <- q +theme_bw() + theme(legend.position = "right") +
    theme(panel.background = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.ontop = TRUE) +
    theme(legend.background = element_rect(colour = "black"))

  #Get the last recorded VL and the desc
  #if StartedART and below vl.cutoff then suppressed otherwise NOT suppressed

  ##VLevent.df <- ddply(VLevent.df,.(ID), tail,1)
  VLevent.df <- VLevent.df[, .SD[c(.N)], by=ID]

  VLevent.df <- VLevent.df %>% mutate(VL.Suppressed.Timepoint = (Description=="Started ART" & Log10VL < vl.cutoff))

  # Now we apply the left_join dplyr function to add the VL status to raw.df.
  InfectedOnTreatment <- left_join(x = InfectedOnTreatment, y = VLevent.df, by = c("ID"))


  vlSuppressedSixmonthLessTP <- data.frame(dplyr::summarise(dplyr::group_by(InfectedOnTreatment, Gender),
                                                            TotalCases = n(),
                                                            VLSuppressed6TP = sum(VL.Suppressed.Timepoint),
                                                            Percentage = sum(VL.Suppressed.Timepoint)/n()))

  vlSuppressedSixmonthLessTP$Gender[vlSuppressedSixmonthLessTP$Gender==0] <- "Woman"
  vlSuppressedSixmonthLessTP$Gender[vlSuppressedSixmonthLessTP$Gender==1] <- "Man"


  return(vlSuppressedSixmonthLessTP)
}

