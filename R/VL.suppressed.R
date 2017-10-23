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
                          timepoint = 30, vlcutoff = 1000,
                          lessmonths = 6, site="All"){

  DTalive.infected <- datalist$ptable

  if (site == "All") {
    DTalive.infected <- subset(DTalive.infected, TOB <= timepoint &
                                 TOD > timepoint & InfectTime!=Inf)
  } else{

    DTalive.infected <- subset(DTalive.infected,
                               TOB <= timepoint &
                                 TOD > timepoint &
                                 pfacility == site &
                                 InfectTime!=Inf)
  }

  vl.cutoff <- log10(vlcutoff)
  six.months <- timepoint - lessmonths/12

  DTalive.infected <- DTalive.infected %>%
    dplyr::mutate(ARTLess6mnths = (TreatTime <= six.months))
  raw.df <- as.data.frame(DTalive.infected)

  InfectedOnTreatment <- subset(raw.df, ARTLess6mnths==TRUE)

  #Check the vl information of these individuals
  VLevent.df <- data.frame(subset(datalist$vltable, ID %in% InfectedOnTreatment$ID ))

  # yaxis <- vl.cutoff + 0.2 # add the axis to make visual
  #
  # #Visualise the VL points
  # q <- ggplot()
  # q <- q + geom_point(data=VLevent.df, aes(x=Time, y=Log10VL, group=ID, colour = Desc))
  # #q <- q + geom_line(data=VLevent.df, aes(x=Time, y=Log10VL, group=ID, colour = Desc))
  # q <- q + geom_hline(yintercept=vl.cutoff) +
  #   annotate("text", datalist$itable$hivseed.time[1], vl.cutoff + 0.2, label = "  VLCutoff")
  # q <- q + geom_vline(xintercept=timepoint-0.5, colour = "blue") +
  #   annotate("text", timepoint, max(VLevent.df$Log10VL), label = "  TimePoint")
  # q <- q +theme_bw() + theme(legend.position = "right") +
  #   theme(panel.background = element_blank(),panel.grid.major.x = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         panel.grid.major.y = element_blank(), panel.ontop = TRUE) +
  #   theme(legend.background = element_rect(colour = "black"))

  #Get the last recorded VL and the desc
  #if StartedART and below vl.cutoff then suppressed otherwise NOT suppressed

  ##VLevent.df <- dplyr(VLevent.df,.(ID), tail,1)
  ##VLevent.df <- VLevent.df[, .SD[c(.N)], by=ID]

  VLevent.df <- VLevent.df %>% mutate(VL.Suppressed.Timepoint =
                                        (Desc=="Started ART" & Log10VL < vl.cutoff))

  # Now we apply the left_join dplyr function to add the VL status to raw.df.
  InfectedOnTreatment <- dplyr::left_join(x = InfectedOnTreatment, y = VLevent.df, by = c("ID"))


  if(nrow(InfectedOnTreatment)==0){

    vlSuppressed.TP <- as.data.frame(matrix(NA, 3, 4))

  }else{
    vlSuppressed.TP <- data.frame(dplyr::summarise(
      dplyr::group_by(InfectedOnTreatment, Gender),
      TotalCases = n(),
      VLSuppressed = sum(VL.Suppressed.Timepoint, na.rm = TRUE),
      Percentage = sum(VL.Suppressed.Timepoint, na.rm = TRUE)/n() *100))


     vlSuppressed.TP <- rbind(vlSuppressed.TP,
                  c("NA", nrow(InfectedOnTreatment),
                  sum(InfectedOnTreatment$VL.Suppressed.Timepoint, na.rm = TRUE),
                  sum(InfectedOnTreatment$VL.Suppressed.Timepoint, na.rm = TRUE)/nrow(InfectedOnTreatment)*100
                  ))

  }

  names(vlSuppressed.TP) <- c("Gender", "TotalCases",
                              "VLSuppressed", "Percentage")


  return(vlSuppressed.TP)
}

