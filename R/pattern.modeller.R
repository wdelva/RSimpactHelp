#' Model age-mixing pattern.
#'
#' Models the age-mixing pattern for the population simulated in Simpact. The
#' user specifies a time point, time window, and age group for which they would
#' like to obtain a summary of the age-mixing pattern.
#'
#' The \code{pattern.modeller()} function uses the \code{link[nlme]{lme}}
#' function to build a linear mixed effects model regressing the partner's age
#' at the time the relationship started on his/her own age at the beginning of
#' the relationship. The models are stratified by gender. Each person can have
#' more than one relationship so a random intercept at the level of the person
#' is added to the model.
#'
#' The models also explicitly account for heteroskedastic variance, since the
#' variance in partner ages tends to grow for older ages.
#'
#' The function produces a list containing two elements. The first is a
#' subsetted dataframe of relationships that meet the inclusion criteria
#' specified in the function's arguments. This dataframe includes
#' population-average partner age predictions from the model. The second element
#' is a dataframe that contains key age-mixing pattern outputs from the model.
#' These outputs can be found in the second dataframe:
#'
#' \describe{ \item{slope} {Also known as the beta-coefficient, this is the
#' change in partner age for each year increase in age. It represents how fast
#' the average age differences grow in the population.} \item{intercept} {This
#' is the average partner age for the first age in the specified interval.}
#' \item{power} {This represents how fast the variance in partner ages grows
#' with increasing age.} \item{lowerpower} {This is the lower limit of the 95%
#' confidence interval for the power coefficient.}
#' \item{upperpower} {This is
#' the upper limit of the 95% confidence interval for the power coefficient.
#' }
#' \item{bvar} {This is the between-subject variance from the model.}
#' \item{wvar} {Also known as the residual error, this represents the amount of
#' within-subject variance from the model.}}
#'
#' @param df The dataframe that is produced by \code{\link{agemix.df.maker()}}
#' @param agegroup Boundaries of the age group that should be retained, e.g.
#'   c(15, 30). The interval is closed on the left and open on the right.
#' @param timepoint Point in time during the simulation to be used in the
#'   calculation.
#' @param timewindow The length of time before the timepoint for which
#'   relationships should be included, e.g. 1, representing one year before the
#'   timepoint. This should be a whole number.
#' @param start This is a logical indicating that only relationships starting
#'   after the beginning of the window should be used. If start = FALSE
#'   relationships could start before the time window. This is the default.
#'
#' @return returns a list with two elements: 1. Dataframe of relationships
#'   inputted into model, and 2. dataframe of model outputs.
#'
#' @examples
#' load(persreldf)
#' agemixpatdat <- pattern.modeller(df = persreldf, agegroup = c(15, 30), timewindow = 1, timepoint = 30)
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom nlme lme
#' @importFrom nlme predict
#' @importFrom nlme intervals
#' @importFrom nlme VarCorr

pattern.modeller <- function(df,
                            agegroup,
                            timepoint,
                            timewindow,
                            start = FALSE) {
  #Warnings
  if (!is.data.frame(df)) {
    stop("df wrong type")
  }

  if (length(agegroup) != 2) {
    stop("Need an upper and lower bound")
  }

  if (timepoint < 1) {
    stop("Time must be at least 1 year")
  }

  if (timewindow < 0) {
    stop("Window must be a whole number")
  }

  #Constants
  time <- timepoint
  window <- timepoint - timewindow
  lwrage <- agegroup[1]
  uprage <- agegroup[2]


  if (start == TRUE) {
    # This only includes relationships that started
    # in the time window

    men <- df %>%
      mutate(age = time - TOB) %>%
      filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "male" &
               TOD > time) %>%
      mutate(agerelform0 = agerelform - lwrage)

    women <- df %>%
      mutate(age = time - TOB) %>%
      filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "female" &
               TOD > time) %>%
      mutate(agerelform0 = agerelform - lwrage)

  } else {
    # This includes all relationships that were ongoing
    # at somepoint during the time window, but may have
    # started long before the time window.

    men <- df %>%
      mutate(age = time - TOB) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "male" &
               TOD > time) %>%
      distinct(ID, relid) %>%
      mutate(agerelform0 = agerelform - lwrage)

    women <- df %>%
      mutate(age = time - TOB) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "female" &
               TOD > time) %>%
      distinct(ID, relid) %>%
      mutate(agerelform0 = agerelform - lwrage)

  }

  malemodel <- lme(pagerelform ~ agerelform0,
                   data = men,
                   random = ~1 | ID,
                   weights = varPower(value = 0.5, form = ~agerelform0 + 1),
                   method = "REML")

  femalemodel <- lme(pagerelform ~ agerelform0,
                     data = women,
                     random = ~1 | ID,
                     weights = varPower(value = 0.5, form = ~agerelform0 + 1),
                     method = "REML")

  men$pred <- predict(malemodel, men, level = 0)
  women$pred <- predict(femalemodel, women, level = 0)

  comb <- bind_rows(men, women)

  modoutput <- matrix(nrow = 1, ncol = 14) %>%
    as.data.frame()
  colnames(modoutput) <- c("slopem", "slopew", "interceptm", "interceptw",
                           "powerm",  "lowerpowerm", "upperpowerm", "powerw",
                           "lowerpowerw", "upperpowerw", "bvarm", "bvarw",
                           "wvarm", "wvarw")

  # Extract slope/Beta-coefficent of the fixed effects from model
  modoutput$slopem <- malemodel$coefficients$fixed[2]
  modoutput$slopew <- femalemodel$coefficients$fixed[2]

  # Extract population intercepts
  modoutput$interceptm <- malemodel$coefficients$fixed[1]
  modoutput$interceptw <- femalemodel$coefficients$fixed[1]

  # Extract power
  modoutput$powerm <- (attributes(malemodel$apVar)$Pars["varStruct.power"])
  modoutput$lowerpowerm <- intervals(malemodel)$varStruct[,1]
  modoutput$upperpowerm <- intervals(malemodel)$varStruct[,3]

  modoutput$powerw <- (attributes(femalemodel$apVar)$Pars["varStruct.power"])
  modoutput$lowerpowerw <- intervals(femalemodel)$varStruct[, 1]
  modoutput$upperpowerw <- intervals(femalemodel)$varStruct[, 3]

  # Extract between-participant variance
  modoutput$bvarm <- VarCorr(malemodel)[1] %>%
    as.numeric()

  modoutput$bvarw <- VarCorr(femalemodel)[1] %>%
    as.numeric()

  # Extract residual variance
  modoutput$wvarm <- VarCorr(malemodel)[2] %>%
    as.numeric()

  modoutput$wvarw <- VarCorr(femalemodel)[2] %>%
    as.numeric()

  agemix.pieces <- list(comb, modoutput)

  return(agemix.pieces)

}
