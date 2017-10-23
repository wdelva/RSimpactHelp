#' Calculate concurrent point prevalence, overall and stratified.
#'
#' Calculate the concurrent point prevalence, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should
#' be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the HIV prevalence should be calculated.
#' @param hivstatus HIV status at the time of the survey. Options are 2, means all; 0
#' means only HIV-negative, 1 means only HIV-positive.
#' @return A dataframe with concurrent point prevalence estimate and surrounding confidence bounds,
#' for the specified time point, age group and hiv status, overall, and stratified by gender
#' @examples
#' data(datalist)
#' concurrent.pointprev.df <- concurr.pointprev.calculator(datalist = datalist, agegroup = c(15,50),
#' timepoint = 10, hivstatus = 2)
#' concurrent.pointprev.df
#'
#' @importFrom stats binom.test
#' @import dplyr
#' @importFrom tibble as_tibble
#' @export

concurr.pointprev.calculator <- function(datalist = datalist,
                                         agegroup = c(15, 50),
                                         timepoint,
                                         hivstatus = 2){
  newtimepoint = timepoint - 0.5
  # We only take the data of people who were alive at timepoint
  DTalive.infected <- alive.infected(datalist = datalist, timepoint = newtimepoint, site = "All")
  agemix <- agemix.df.maker(datalist = datalist)

  # for males
  DTalive.infected.agegroup.men <- subset(DTalive.infected,
                                          TOB <= newtimepoint - agegroup[1] &
                                            TOB > newtimepoint - agegroup[2] &
                                            Gender == 0)

  degree.df <- degree.df.maker(dataframe.df = agemix,
                               agegroup = agegroup,
                               hivstatus = hivstatus,
                               survey.time = newtimepoint,
                               window.width = 0,
                               gender.degree = "male",
                               only.new = FALSE)

  number.males.with.cps <- sum(degree.df$Degree > 1)
  popsize.males <- nrow(DTalive.infected.agegroup.men)

  if(!is.na(number.males.with.cps)){
   concurr.pointprev.males <- number.males.with.cps / popsize.males
    pointprevalence.95.ll <- binom.test(x = number.males.with.cps, n = popsize.males)$conf.int[1]
    pointprevalence.95.ul <- binom.test(x = number.males.with.cps, n = popsize.males)$conf.int[2]

    concurr.pointprev.males.df <- cbind.data.frame(Gender = "male",
                                                   concurr.pointprev = concurr.pointprev.males,
                                                   pointprevalence.95.ll = pointprevalence.95.ll,
                                                   pointprevalence.95.ul = pointprevalence.95.ul)
  } else{

    concurr.pointprev.males.df <- cbind.data.frame(Gender = "male",
                                                   concurr.pointprev = NA,
                                                   pointprevalence.95.ll = NA,
                                                   pointprevalence.95.ul = NA)
  }

  # for females
  DTalive.infected.agegroup.female <- subset(DTalive.infected,
                                             TOB <= newtimepoint - agegroup[1] &
                                               TOB > newtimepoint - agegroup[2] &
                                               Gender == 1)

  degree.df <- degree.df.maker(dataframe.df = agemix,
                               agegroup = agegroup,
                               hivstatus = hivstatus,
                               survey.time = newtimepoint,
                               window.width = 0,
                               gender.degree = "female",
                               only.new = FALSE)

  number.females.with.cps <- sum(degree.df$Degree > 1)
  popsize.females <- nrow(DTalive.infected.agegroup.female)

  if(!is.na(number.females.with.cps)){
    concurr.pointprev.females <- number.females.with.cps / popsize.females
    pointprevalence.95.ll <- binom.test(x = number.females.with.cps, n = popsize.females)$conf.int[1]
    pointprevalence.95.ul <- binom.test(x = number.females.with.cps, n = popsize.females)$conf.int[2]

    concurr.pointprev.females.df <- cbind.data.frame(Gender = "female",
                                                     concurr.pointprev = concurr.pointprev.females,
                                                     pointprevalence.95.ll = pointprevalence.95.ll,
                                                     pointprevalence.95.ul = pointprevalence.95.ul)
  } else{

    concurr.pointprev.females.df <- cbind.data.frame(Gender = "female",
                                                     concurr.pointprev = NA,
                                                     pointprevalence.95.ll = NA,
                                                     pointprevalence.95.ul = NA)
  }

  # for all

  number.people.with.cps <- sum(number.males.with.cps, number.females.with.cps)
  popsize <- sum(popsize.males, popsize.females)

  if(!is.na(number.people.with.cps)){
    concurr.pointprev <- number.people.with.cps / popsize
    pointprevalence.95.ll <- binom.test(x = number.people.with.cps, n = popsize)$conf.int[1]
    pointprevalence.95.ul <- binom.test(x = number.people.with.cps, n = popsize)$conf.int[2]

    concurr.pointprev.all.df <- cbind.data.frame(Gender = "all",
                                                 concurr.pointprev = concurr.pointprev,
                                                 pointprevalence.95.ll = pointprevalence.95.ll,
                                                 pointprevalence.95.ul = pointprevalence.95.ul)
  } else{

    concurr.pointprev.all.df <- cbind.data.frame(Gender = "all",
                                                 concurr.pointprev = NA,
                                                 pointprevalence.95.ll = NA,
                                                 pointprevalence.95.ul = NA)
  }
  concurr.pointprev.df <- rbind.data.frame(concurr.pointprev.males.df, concurr.pointprev.females.df,concurr.pointprev.all.df)
  return(as_tibble(concurr.pointprev.df))
}


