#' Calculate HIV prevalence, overall and stratified.
#'
#' Calculate the HIV prevalence at a point in time, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the HIV prevalence should be calculated.
#' @return a dataframe with prevalence estimate and surrounding confidence bounds,
#' for the specified time point and age group, overall, and stratified by gender
#' @examples
#' prevalence.df <- prevalence.calculator(datalist = datalist, agegroup = c(15, 30), timepoint = 30)
#' @import dplyr

prevalence.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){

  # First we only take the data of people who were alive at the timepoint
  DTalive.infected <- alive.infected(datalist = datalist, timepoint = timepoint)

  DTalive.infected.agegroup <- subset(DTalive.infected, TOB <= timepoint - agegroup[1] & TOB > timepoint - agegroup[2])

  raw.df <- data.frame(DTalive.infected.agegroup)
  if(nrow(raw.df)>0 & sum(raw.df$Infected)>0){
    # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
    prevalence.df <- dplyr::summarise(group_by(raw.df, Gender),
                                     popsize = length(Gender),
                                     sum.cases = sum(Infected),
                                     pointprevalence = sum(Infected) / length(Gender),
                                     pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                     pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
                                     )
    prevalence.all.df <- cbind(Gender = NA,
                              dplyr::summarise(raw.df,
                                               popsize = length(Gender),
                                               sum.cases = sum(Infected),
                                               pointprevalence = sum(Infected) / length(Gender),
                                               pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                               pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
                                               )
                              )
    prevalence.df <- rbind(prevalence.df, prevalence.all.df)
  } else {

    prevalence.df <- data.frame(Gender = c(NA,NA,NA),
                                popsize = c(NA,NA,NA),
                                sum.cases = c(NA,NA,NA),
                                pointprevalence = c(NA,NA,NA),
                                pointprevalence.95.ll = c(NA,NA,NA),
                                pointprevalence.95.ul = c(NA,NA,NA)
                                )
  }

  return(prevalence.df)
}

prevalence.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){
  #DTP <- datalist$ptable
  DTalive.infected <- alive.infected(datalist = datalist, timepoint = timepoint) # First we only take the data of people who were alive at the timepoint
  DTalive.infected.agegroup <- subset(DTalive.infected, TOB <= timepoint - agegroup[1] & TOB > timepoint - agegroup[2])

  raw.df <- data.frame(DTalive.infected.agegroup)
  if(nrow(raw.df)>0 & sum(raw.df$Infected)>0){

  # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
  prevalence.df <- dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                   popsize = length(Gender),
                                   sum.cases = sum(Infected),
                                   pointprevalence = sum(Infected) / length(Gender),
                                   pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                   pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
                                   )
  prevalence.all.df <- cbind(Gender = NA,
                            dplyr::summarise(raw.df,
                                             popsize = length(Gender),
                                             sum.cases = sum(Infected),
                                             pointprevalence = sum(Infected) / length(Gender),
                                             pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                             pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
                                             )
                            )
  prevalence.df <- rbind(prevalence.df, prevalence.all.df)
  } else {

    prevalence.df <- data.frame(Gender = c(NA,NA,NA),
                                popsize = c(NA,NA,NA),
                                sum.cases = c(NA,NA,NA),
                                pointprevalence = c(NA,NA,NA),
                                pointprevalence.95.ll = c(NA,NA,NA),
                                pointprevalence.95.ul = c(NA,NA,NA)
                                )
  }

  return(prevalence.df)
}
