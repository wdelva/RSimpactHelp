#' We expect this to grow and the returned summary statistics to include a lot more.
#' Will include an Ignore variable to each of the summary statistics. Only computed if IGNORE is set to FALSE
#'
#' After the simpact simulation a number of statistics can be computed to summarise the output
#'
#' art.intro list should include the time when the diagnosis to take effect.
#'
#'
#' The function take in a list of paramenters set at intervals of cd4 count threshhold used to
#' decide if a person can be offerd treatment
#'
#' @param datalist.test The dataframe that is produced by \code{\link{readthedata()}}
#' @param growth.rate parameters used  by \code{\link{pop.growth.calculator()}}
#' @param agemix.maker parameters used  by \code{\link{pattern.modeller()}} with input from \code{\link{agemix.df.make()}}
#' @param prev.15.25 parameters used  by \code{\link{prevalence.calculator()}} of the 15-25 age group
#' @param prev.25.50 parameters used  by \code{\link{prevalence.calculator()}} of the 25-50 age group
#' @param art.coverage parameters used  by \code{\link{ART.coverage.calculator()}}
#' @param inc.15.30 parameters used  by \code{\link{incidence.calculator()}}
#' @param partner.degree  parameters used  by \code{\link{degree.df.maker()}}
#' @return a matrix with summary statistics segregated by age and gender and different time points within the simulation period
#' @examples
#' out.test <- output.summary.maker(datalist.test, growth.rate=list(0,20), agemix.maker=list(15,30,30,1,FALSE,female),
#' prev.15.25 = list(15,25,35,men), prev.25.50 = list(25,50,35,men), art.coverage = list(15,50,35,men), inc.15.30=list())

output.summary.maker <- function(datalist = datalist.test, growth.rate=list(timewindow.min = 0, timewindow.max = 20),
                                   agemix.maker=list(agegroup.min = 15, agegroup.max=30, timepoint =30,
                                                     timewindow = 1, start=FALSE, gender = "female"),
                                   prev.15.25 = list(age.group.min=15, age.group.max=25, timepoint = 35, gender = "men"),
                                   prev.25.50 = list(age.group.min=25, age.group.max=50, timepoint = 35, gender = "men"),
                                   art.coverage = list(age.group.min=15, age.group.max=50, timepoint = 35, gender = "men"),
                                   inc.15.30 = list(age.group.min=15, age.group.max=30, timewindow.min = 30,
                                                      timewindow.max = 40, gender = "men", only.active = "No"),
                                   partner.degree = list(age.group.min=15, age.group.max=30, hivstatus = 0, survey.time = 30,
                                                         window.width = 1, gender="female", only.new = FALSE)){

  prev.15.25.gender.index = 1
  prev.25.50.gender.index = 1
  inc.15.30.gender.index = 1
  art.coverage.gender.index = 1

  if(prev.15.25$gender!="men"){prev.15.25.gender.index = 2}
  if(prev.25.50$gender!="men"){prev.25.50.gender.index = 2}
  if(inc.15.30$gender!="men"){inc.15.30.gender.index = 2}
  #if(agemix.maker$gender!="female"){agemix.gender = "male"}
  if(art.coverage$gender!="men"){art.coverage.gender.index = 2}

  #summary statistics #################
  # 1. Population growth
  growth.rate <- pop.growth.calculator(datalist = datalist,
                                       timewindow = c(growth.rate$timewindow.min, growth.rate$timewindow.max))
  # 2. Median age difference in most recent relationships (target 2-5 years)
  agemix.df <- agemix.df.maker(datalist)
  pattern <- pattern.modeller(dataframe = agemix.df,
                              agegroup = c(agemix.maker$agegroup.min,agemix.maker$agegroup.max),
                              timepoint = agemix.maker$timepoint,
                              timewindow = agemix.maker$timewindow,
                              start = agemix.maker$start)
  # agedifmedtable$median[2] # Median age difference as reported by women
  median.AD <- as.numeric(median(pattern[[1]]$AgeGap[pattern[[1]]$Gender == agemix.maker$gender]))

  # 3. IQR age difference in most recent relationships (target IQR 2 - 6)
  Q1.AD <- as.numeric(summary(pattern[[1]]$AgeGap[pattern[[1]]$Gender == agemix.maker$gender])[2])
  Q3.AD <- as.numeric(summary(pattern[[1]]$AgeGap[pattern[[1]]$Gender == agemix.maker$gender])[5])
  #IQR.AD <- agedifmedtable$Q3[2] - agedifmedtable$Q1[2] # IQR as reported by women

  # 4. HIV prevalence in the window among men 15-25 (target 5% - 10%)
  prev <- prevalence.calculator(datalist = datalist, agegroup = c(prev.15.25$age.group.min,
                                                                       prev.15.25$age.group.max), timepoint = prev.15.25$timepoint)
  prev.men.15.25 <- prev$pointprevalence[prev.15.25.gender.index] # among men

  # 5. HIV prevalence in the window among men 25-50 (target 10% - 30%)
  prev <- prevalence.calculator(datalist = datalist, agegroup = c(prev.25.50$age.group.min,
                                                                       prev.25.50$age.group.max), timepoint = prev.25.50$timepoint)
  prev.men.25.50 <- prev$pointprevalence[prev.25.50.gender.index] # among men

  # 6. Point prevalence of concurrency. postsim function to be converted to RSimpactHelper function

  # 7. ART coverage among all HIV+ people aged 15-50 (target 15% - 40% in 2011)
  ARTcov <- ART.coverage.calculator(datalist = datalist, agegroup = c(art.coverage$age.group.min,
                                                                           art.coverage$age.group.max), timepoint = art.coverage$timepoint) # 2011 is 34 years after 1977
  ART.cov.15.50 <- ARTcov$ART.coverage[art.coverage.gender.index] # among men and women

  # 8. HIV incidence among women 15 <= x < 30 in the window 20-30 years after HIV introduction
  inc <- incidence.calculator(datalist = datalist, agegroup = c(inc.15.30$age.group.min, inc.15.30$age.group.max),
                              timewindow = c(inc.15.30$timewindow.min, inc.15.30$timewindow.max), only.active = inc.15.30$only.active)
  incid.wom.15.30 <- inc$incidence[inc.15.30.gender.index] # among women

  # NOTE: We may want to also calculate a different type of HIV incidence, more in line with the Harling paper:
  # Where we only accumulate exposure time for the periods that a woman was in (a) relationship(s).
  # For now, this (correct) HIV incidence measure suffices

  # 9. % of women who had > 1 partner in the past 12 months
  degree.df <- degree.df.maker(agemix.df, agegroup = c(partner.degree$age.group.min, partner.degree$age.group.max),
                              hivstatus = partner.degree$hivstatus, survey.time = partner.degree$survey.time,
                              window.width = partner.degree$window.width, gender =  partner.degree$gender,
                              only.new = partner.degree$only.new)
  frac.degreeGT1.wom.15.30 <- mean(degree.df$Degree > 1)


  out.test <- matrix(c(growth.rate, median.AD, Q1.AD, Q3.AD, prev.men.15.25, prev.men.25.50, ART.cov.15.50,
                       incid.wom.15.30, frac.degreeGT1.wom.15.30), nrow = 1,
                     dimnames = list(NULL, c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25",
                                             "prev.men.25.50", "ART.cov.15.50", "incid.men.15.30",
                                             "frac.degreeGT1.wom.15.30"))) # 9 summary statistics

  if(prev.15.25$gender!="men"){
    colnames(out.test)[colnames(out.test)=="prev.men.15.25"] <- "prev.wom.15.25"
  }
  if(prev.25.50$gender!="men"){
    colnames(out.test)[colnames(out.test)=="prev.men.25.50"] <- "prev.wom.25.50"
  }
  if(inc.15.30$gender!="men"){
    colnames(out.test)[colnames(out.test)=="incid.men.15.30"] <- "incid.wom.15.30"
  }

  return(out.test)
}


