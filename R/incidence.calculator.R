#' Calculate HIV incidence, overall and stratified.
#'
#' Produces a data frame that contains the overall, and gender-stratified,
#' HIV incidence at a specified time window and for a specific
#' age group.
#'
#' @param datalist The list object that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper 
#'   bound) that should be retained. Should be expressed as a vector of two 
#'   integers. e.g. c(15, 30)
#' @param timewindow of the time window (lower bound < time <= upper bound) that
#'   should be retained. Should be expressed as a vector of two integers. e.g.
#'   c(20, 30)
#' @param only.active Should only women who are in sexual relationships 
#'   contribute exposure time? If "Strict", all time spent being single will be
#'   excluded from exposure time. If "Harling", time will be excluded from
#'   exposure time in blocks of one year, if the person spent that entire block
#'   not in any relationship. If "No", exposure time is being contributed, even
#'   while not in any relationships. WARNING: The "Harling" and "Strict" 
#'   estimates currently only work for women. So the resulting data frame will
#'   still contain incidence estimates for men as if you selected "No".
#' @return A dataframe with cases, exposure times, incidence estimates and 
#'   surrounding confidence bounds, for the specified time window and age group,
#'   overall, and stratified by gender.
#' @examples
#' cfg <- list() 
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp") 
#' dl <- readthedata(modeloutput) 
#' incidence.df <- incidence.calculator(datalist = datalist, agegroup = c(15, 30), timewindow = c(20, 30))
#' incidence.df
#'
#' @importFrom exactci poisson.exact
#' @import dplyr
#' @export

incidence.calculator <- function(datalist = datalist,
                                 agegroup = c(15, 30),
                                 timewindow = c(20, 30),
                                 only.active = "No"){
  
  # Incidence is calculated in 3 steps.
  # 1. Calculate PY of exposure per person
  # 2. Calculate whether the person had the event or not
  # 3. Divide events by sum of PY.
  
  # 1. Calculation of PY of exposure per person. Exposure time starts at the max
  # of the lower bound of the timewindow and the time at which the person reaches
  # the lower bound of the age group. Exposure time ends at the min of the upper
  # bound of the timewindow, the time at which the person reaches the upper bound
  # of the age group, and the time at which the person gets infected with HIV.
  # Any negative exposure time will be reset to zero.
  
  # Scalars
  lwr.agegroup <- agegroup[1]
  upr.agegroup <- agegroup[2]
  lwr.timewindow <- timewindow[1]
  upr.timewindow <- timewindow[2]
  
  # Data frames from the list
  df <- datalist$ptable
  
  # Calculate naive exposure time
  # 1. Determine the time person is X yo
  # 2. Determine the time person is Y yo
  # 3. Their exposure starts whichever comes last: lwr bound age, or lwr window
  # 4. Their exposure ends whichever comes 1st: upr bound age, upr window, or
  # infected with HIV
  # 5. Subtract exposure start from exposure end
  df1 <- df %>%
    dplyr::mutate(time.lwr.agegroup = TOB + lwr.agegroup,
                  time.upr.agegroup = TOB + upr.agegroup,
                  exposure.start = pmax(time.lwr.agegroup, lwr.timewindow),        
                  exposure.end = pmin(InfectTime, 
                                      time.upr.agegroup, 
                                      upr.timewindow),
                  exposure.time = exposure.end - exposure.start)
  
  # If the exposure time is negative, then they had no exposure during the
  # timewindow or agegroup of interest
  df2 <- df1 %>%
    dplyr::mutate(had.exposure.time = exposure.time > 0,
                  exposure.time = ifelse(exposure.time < 0, 
                                         0, 
                                         exposure.time))
  
  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, 
  # and before or at exposure.end
  df3 <- df2 %>%
    dplyr::mutate(incident.case = (InfectTime > exposure.start & 
                                     InfectTime <= exposure.end))
  
  # Create a copy for the special case below
  df3.plus <- df3 
  
  # If the user specifies that they want "Strict" or "Harling" definition
  # of exposure time 
  if (only.active != "No"){
    
    # timespentsingle.calculator() calculates time not in a relationship
    # for each woman
    norels.timespent.df <- timespentsingle.calculator(datalist = datalist,
                                                      agegroup = agegroup,
                                                      timewindow = timewindow,
                                                      type = only.active)
    
    # Combine the person-exposure dataset to the time spent single dataset
    df3.plus <- df3.plus %>%
      dplyr::left_join(norels.timespent.df,
                       by = c("ID" = "woman.ID"))
    
    # Modify variables for:
    # 1. If there wasn't any time spent single, replace with 0
    # 2. Update exposure time by subtracting the time spent single
    df3.plus <- df3.plus %>%
      dplyr::mutate(sum.norels.timespent = ifelse(is.na(sum.norels.timespent), # 1
                                                  0, 
                                                  sum.norels.timespent),
                    exposure.time = exposure.time - sum.norels.timespent) # 2
    
  }
  
  # Keep people who had exposure time
  df4.plus <- df3.plus %>%
    dplyr::filter(exposure.time >= 0)
  
  if(nrow(df4.plus) > 0){
    
    # Create summary table of incidence by gender
    incidence.df <- df4.plus %>%
      dplyr::group_by(Gender) %>% # Stratify by gender
      dplyr::summarise(sum.exposure.time = sum(exposure.time), # Add all exposure time by gender
                       sum.incident.cases = sum(incident.case)) %>% # Add all incident cases by gender
      dplyr::mutate(incidence = sum.incident.cases / sum.exposure.time) %>% # Calc incidence
      dplyr::group_by(Gender, sum.exposure.time, sum.incident.cases) %>% # Need to group, so each row is a group.
      dplyr::mutate(incidence.95.ll = exactci::poisson.exact(sum.incident.cases, sum.exposure.time)$conf.int[1], # Because poisson.exact is 
                    incidence.95.ul = exactci::poisson.exact(sum.incident.cases, sum.exposure.time)$conf.int[2]) # not a vectorized function
    
    # Now overall incidence 
    incidence.all.df <- df4.plus %>%
      dplyr::summarise(Gender = NA, # Create column for gender
                       sum.exposure.time = sum(exposure.time),
                       sum.incident.cases = sum(incident.case)) %>%
      dplyr::mutate(incidence = sum.incident.cases / sum.exposure.time) %>%
      dplyr::group_by(sum.exposure.time, sum.incident.cases) %>%
      dplyr::mutate(incidence.95.ll = exactci::poisson.exact(sum.incident.cases, 
                                                             sum.exposure.time)$conf.int[1],
                    incidence.95.ul = exactci::poisson.exact(sum.incident.cases, 
                                                             sum.exposure.time)$conf.int[2])
    
    
    # Combine stratified, and overall incidence
    incidence.df <- bind_rows(incidence.df, incidence.all.df) %>%
      ungroup()
    
  }else{ # In the rare event that a simulated data set has zero observations
    
    incidence.df <- data.frame(Gender = c(NA,NA,NA),
                               sum.exposure.time = c(NA,NA,NA),
                               sum.incident.cases = c(NA,NA,NA),
                               incidence = c(NA,NA,NA),
                               incidence.95.ll = c(NA,NA,NA),
                               incidence.95.ul = c(NA,NA,NA))
  }
  
  return(incidence.df)
  
}