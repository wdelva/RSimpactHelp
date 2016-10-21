#' Calculate HIV prevalence and ART coverage, overall and stratified.
#'
#' Calculate the HIV prevalence and ART coverage at a point in time, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the ART coverage should be calculated.
#' @return a dataframe with HIV prevalence estimates and ART coverage estimate and surrounding confidence bounds,
#' for the specified time point and age group, overall, and stratified by gender
#' @examples
#' ART.coverage.df <- ART.coverage.calculator(datalist = datalist, agegroup = c(15, 30), timepoint = 30)
#' @import dplyr

ART.coverage.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){
  DTP <- datalist$ptable
  DTalive.infected <- alive.infected(DT = DTP, timepoint = timepoint) # First we only take the data of people who were alive at the timepoint
  DTalive.infected.agegroup <- subset(DTalive.infected, TOB <= timepoint - agegroup[1] & TOB > timepoint - agegroup[2])
  raw.df <- data.frame(DTalive.infected.agegroup)
  art.df <- subset(datalist$ttable, TStart <= timepoint & TEnd > timepoint)
  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  if (nrow(raw.df > 0) & sum(raw.df$Infected)>0) {
    raw.df$onART <- !is.na(raw.df$TStart)

    # Now we apply some dplyr function to get the sum of cases and population size per gender.
    # ART.coverage.df <- dplyr::summarise(group_by(raw.df, Gender),
    #                                  popsize = length(Gender),
    #                                  sum.cases = sum(Infected),
    #                                  sum.onART = sum(onART),
    #                                  pointprevalence = sum(Infected) / length(Gender),
    #                                  pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
    #                                  pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2],
    #                                  ART.coverage = sum(onART) / sum(Infected),
    #                                  ART.coverage.95.ll = as.numeric(binom.test(x = sum(onART), n = sum(Infected))$conf.int)[1],
    #                                  ART.coverage.95.ul = as.numeric(binom.test(x = sum(onART), n = sum(Infected))$conf.int)[2]
    #                                  )
    ART.coverage.all.df <- cbind(Gender = NA,
                                 dplyr::summarise(raw.df,
                                                  popsize = length(Gender),
                                                  sum.cases = sum(Infected),
                                                  sum.onART = sum(onART),
                                                  pointprevalence = sum(Infected) / length(Gender),
                                                  pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                                  pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2],
                                                  ART.coverage = sum(onART) / sum(Infected),
                                                  ART.coverage.95.ll = as.numeric(binom.test(x = sum(onART), n = sum(Infected))$conf.int)[1],
                                                  ART.coverage.95.ul = as.numeric(binom.test(x = sum(onART), n = sum(Infected))$conf.int)[2]

                                 )
    )
    ART.coverage.df <- ART.coverage.all.df #rbind(ART.coverage.df, ART.coverage.all.df)
  } else {
    ART.coverage.df <- data.frame(Gender = NA,
                                  popsize = NA,
                                  sum.cases = NA,
                                  sum.onART = NA,
                                  pointprevalence = NA,
                                  pointprevalence.95.ll = NA,
                                  pointprevalence.95.ul = NA,
                                  ART.coverage = NA,
                                  ART.coverage.95.ll = NA,
                                  ART.coverage.95.ul = NA)
  }
  return(ART.coverage.df)
}
