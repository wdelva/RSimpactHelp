#' Calculate HIV prevalence, overall and stratified.
#'
#' Calculate the HIV prevalence at a point in time, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should
#' be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the HIV prevalence should be calculated.
#' @return a dataframe with prevalence estimate and surrounding confidence bounds,
#' for the specified time point and age group, overall, and stratified by gender
#' @examples
#' data(datalist)
#' prevalence.df <- prevalence.calculator(datalist = datalist, agegroup = c(15, 30),
#' timepoint = 30)
#' prevalence.df
#'
#' @importFrom stats binom.test
#' @import dplyr
#' @export

prevalence.calculator <- function(datalist = datalist, 
                                  agegroup = c(15, 30), 
                                  timepoint = 30){
  
  # Using alive.infect() take a subset who were alive at the timepoint
  df.alive.infect <- alive.infected(datalist = datalist, 
                                    timepoint = timepoint)
  
  # Retain only those who are in the specified age groups
  df.alive.infect <- df.alive.infect %>%
    dplyr::mutate(age = timepoint - TOB) %>%
    dplyr::filter(age >= agegroup[1],
                  age < agegroup[2])
  
  
  if(nrow(df.alive.infect) > 0){
    
    # Create summary table of prevalence by gender
    prevalence.df <- df.alive.infect %>%
      dplyr::group_by(Gender) %>% # Stratify data by gender
      dplyr::summarise(popsize = n(), # Total observations for each gender
                       sum.cases = sum(Infected)) %>% # Total cases for each gender
      dplyr::mutate(pointprevalence = sum.cases / popsize) %>% # New var for prevalence
      dplyr::group_by(Gender, popsize, sum.cases) %>% # Need to group, so each row is a group. 
      dplyr::mutate(pointprevalence.95.ll = binom.test(sum.cases, popsize)$conf.int[1], # binom.test is not vectorized
                    pointprevalence.95.ul = binom.test(sum.cases, popsize)$conf.int[2]) # Upper limit
    
    prevalence.all.df <- df.alive.infect %>%
      dplyr::summarise(Gender = NA, # Create column for gender
                       popsize = n(),
                       sum.cases = sum(Infected)) %>%
      dplyr::mutate(pointprevalence = sum.cases / popsize) %>%
      dplyr::group_by(popsize, sum.cases) %>%
      dplyr::mutate(pointprevalence.95.ll = binom.test(sum.cases, popsize)$conf.int[1],
                    pointprevalence.95.ul = binom.test(sum.cases, popsize)$conf.int[2]) 
    
    
    prevalence.df <- bind_rows(prevalence.df, prevalence.all.df) # Combine stratified, and overall prev
    
  } else { # In the rare event that a simulated data set has zero observations
    
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

