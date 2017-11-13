#' Calculate overall HIV prevalence and ART coverage aggregated by gender.
#'
#' Calculate the HIV prevalence and ART coverage at a point in time, for specific
#' age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower.bound <= age < upper.bound) that
#' should be retained, e.g. agegroup = c(15, 30)
#' @param timepoint Point in time at which the ART coverage should be calculated.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a dataframe with HIV prevalence estimates and ART coverage estimate and surrounding
#' confidence bounds,
#' for the specified time point and age group, overall, and stratified by gender
#' @examples
#' data(datalist)
#' ART.coverage.df <- ART.coverage.calculator(datalist = datalist,
#' agegroup = c(15, 30), timepoint = 30, site="All")
#'
#' @importFrom magrittr %>%
#' @importFrom stats binom.test
#' @import dplyr
#' @export

ART.coverage.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30, site="All"){

  # First we only take the data of people who were alive at the timepoint

  DTalive.infected <- alive.infected(datalist = datalist,
                                     timepoint = timepoint, site = site)

  DTalive.infected <- DTalive.infected %>%
                      dplyr::filter(Infected==TRUE)


  #Limit the list to the one that match the age group.
  DTalive.infected.agegroup <- DTalive.infected %>%
                              dplyr::filter(TOB <= timepoint - agegroup[1] &
                                            TOB > timepoint - agegroup[2])

  raw.df <- data.frame(DTalive.infected.agegroup)

  art.df <- datalist$ttable %>%
            dplyr::filter(TStart <= timepoint & TEnd > timepoint)

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  if (nrow(raw.df) > 0 & sum(raw.df$Infected, na.rm = TRUE)>0) {
    raw.df$onART <- !is.na(raw.df$TStart)

    #Now we apply some dplyr function to get the sum of cases and population size per gender.
    raw.df$Gender <- as.factor(raw.df$Gender)

    ART.coverage.df <- raw.df %>%
      dplyr::group_by(Gender) %>%
      dplyr::summarise(popsize = length(Gender),
                       sum.cases = sum(Infected, na.rm = TRUE),
                       sum.onART = sum(onART, na.rm = TRUE),
                       pointprevalence = sum(Infected, na.rm = TRUE) / length(Gender),
                       pointprevalence.95.ll = as.numeric(
                         binom.test(x=sum(Infected, na.rm = TRUE),
                                    n=length(Gender))$conf.int)[1],
                       pointprevalence.95.ul = as.numeric(
                         binom.test(x=sum(Infected, na.rm = TRUE),
                                    n=length(Gender))$conf.int)[2],
                       ART.coverage = sum(onART, na.rm = TRUE) / sum(Infected, na.rm = TRUE),
                       ART.coverage.95.ll = as.numeric(
                         binom.test(x=sum(onART, na.rm = TRUE),
                                    n=sum(Infected, na.rm = TRUE))$conf.int)[1],
                       ART.coverage.95.ul = as.numeric(
                         binom.test(x=sum(onART, na.rm = TRUE),
                                    n=sum(Infected, na.rm = TRUE))$conf.int)[2]
                                     )

    #ART.coverage.df <- ART.coverage.all.df #rbind(ART.coverage.df, ART.coverage.all.df)
  } else {
    ART.coverage.df <- data.frame(Gender = c(NA, NA, NA),
                                  popsize = c(NA, NA, NA),
                                  sum.cases = c(NA, NA, NA),
                                  sum.onART = c(NA, NA, NA),
                                  pointprevalence = c(NA, NA, NA),
                                  pointprevalence.95.ll = c(NA, NA, NA),
                                  pointprevalence.95.ul = c(NA, NA, NA),
                                  ART.coverage = c(NA, NA, NA),
                                  ART.coverage.95.ll = c(NA, NA, NA),
                                  ART.coverage.95.ul = c(NA, NA, NA))
  }
  return(ART.coverage.df)
}
