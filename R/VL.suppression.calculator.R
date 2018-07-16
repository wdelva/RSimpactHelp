#' Calculate viral suppression fraction.
#'
#' Calculate the fraction of HIV infected people of a particular age group who
#' are virally suppressed at a point in time.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower.bound <= age <
#'   upper.bound) that should be retained, e.g. agegroup = c(15, 30)
#' @param timepoint Point in time at which the viral suppression fraction should
#'   be calculated.
#' @param vl.cutoff viral load below this threshold e.g <1000 copies/mL is
#'   defined as suppressed
#' @param site Select only the particular site from the study, if all ignore
#'   site/use all sites.
#' @return a dataframe with VL suppression estimates and surrounding confidence
#'   bounds, for the specified time point and age group, overall, and stratified
#'   by gender.
#' @examples
#' data(datalist)
#' VL.suppression.df <- VL.suppression.calculator(datalist = datalist,
#' agegroup = c(15, 30), timepoint = 30, vl.cutoff = 1000, site="All")
#'
#' @importFrom magrittr %>%
#' @importFrom stats binom.test
#' @import dplyr
#' @export

VL.suppression.calculator <- function(datalist = datalist,
                                    agegroup = c(15, 30),
                                    timepoint = 30,
                                    vl.cutoff = 1000,
                                    site = "All"){

  # First we only take the data of people who were alive at the timepoint

  DTalive.infected <- alive.infected(datalist = datalist,
                                     timepoint = timepoint,
                                     site = site)

  DTalive.infected <- DTalive.infected %>%
    dplyr::filter(Infected==TRUE)


  #Limit the list to those who match the age group.
  DTalive.infected.agegroup <- DTalive.infected %>%
    dplyr::filter(TOB <= timepoint - agegroup[1] &
                    TOB > timepoint - agegroup[2])

  raw.df <- data.frame(DTalive.infected.agegroup)

  vl.df <- datalist$vltable %>%    # a dataframe with the most recent VL for each person.
    dplyr::select(c(Time, ID, Log10VL)) %>%
    dplyr::filter(Time <= timepoint) %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(Time) %>%
    slice(n())

  # Now we apply the left_join dplyr function to add the VL status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = vl.df, by = "ID")

  if (nrow(raw.df) > 0 & sum(raw.df$Infected)>0) {
    #Now we apply some dplyr function to get the sum of cases and population size per gender.
    raw.df$Gender <- as.character(raw.df$Gender)
    raw.df <- mutate(raw.df,
                     vl.suppr = Log10VL < log10(vl.cutoff))


    VL.suppression.df <- dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                        popsize = length(Gender),
                                        sum.cases = sum(Infected),
                                        sum.vl.suppr = sum(vl.suppr),
                                        vl.suppr.frac = sum(vl.suppr) / sum(Infected),
                                        vl.suppr.frac.95.ll = as.numeric(
                                          binom.test(x=sum(vl.suppr),n=sum(Infected))$conf.int)[1],
                                        vl.suppr.frac.95.ul = as.numeric(
                                          binom.test(x=sum(vl.suppr),n=sum(Infected))$conf.int)[2]
    )

    #ART.coverage.df <- ART.coverage.all.df #rbind(ART.coverage.df, ART.coverage.all.df)
  } else {
    VL.suppression.df <- data.frame(Gender = NA,
                                  popsize = NA,
                                  sum.cases = NA,
                                  sum.vl.suppr = NA,
                                  vl.suppr.frac = NA,
                                  vl.suppr.frac.95.ll = NA,
                                  vl.suppr.frac.95.ul = NA
    )
  }
  return(VL.suppression.df)
}
