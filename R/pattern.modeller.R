#' Age-mixing pattern
#'
#'Models the age-mixing pattern in the population for relationships that
#'occurred within a specified time window at a specified point in the simulation.
#'
#' @param df The dataframe that is produced by \code{\link{agedif.df.maker()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the age-mixing metrics should be calculated.
#' @param timewindow The length of time before the timepoint for which relationships should be included,
#' e.g. 1, representing one year before the timepoint. This should be a whole number.
#' @param start This is a logical indicating that only relationships starting after the beginning of the window
#' should be used. If start = FALSE relationships could start before the time window.
#'
#' @return a list that contains two elements. One, a subsetted dataframe of relationships that
#' meet criteria of the function's arguments, including partner age predictions from
#' the model, and two, a dataframe that contains key age-mixing pattern
#' outputs from the model.
#'
#' @examples
#' amp <- pattern.modeller(df = dataframe, agegroup = c(15, 30), timewindow = 1, timepoint = 30)

# dplyr, magrittr, nlme

pattern.modeller <- function(df,
                            agegroup,
                            timepoint,
                            timewindow,
                            start = FALSE) {
  #Warnings
  if (!is.data.frame(df)) {
    stop("Dataframe wrong type")
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

  # Create subset df of relationships that meet input criteria
  # Also remove duplicated relationships (from multiple episodes)
  # by gender only
  # The same relationship will still be represented twice —once by the female
  # partner and once by the male partner.
  # However, there will only be one episode for each.

  if (start == TRUE) {

    men <- df %>%
      arrange(ID, relid, FormTime) %>%
      group_by(ID, relid) %>%
      mutate(episodeorder = row_number(),
             agerelform = TOB + FormTime,
             agerelform = first(agerelform),
             pagerelform = agerelform - AgeGap) %>%
      ungroup() %>%
      mutate(age = TOB + time) %>%
      filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "male") %>%
      mutate(agerelform0 = agerelform - 18)

    women <- df %>%
      arrange(ID, relid, FormTime) %>%
      group_by(ID, relid) %>%
      mutate(episodeorder = row_number(),
             agerelform = TOB + FormTime,
             agerelform = first(agerelform),
             pagerelform = agerelform + AgeGap) %>%
      ungroup() %>%
      mutate(age = TOB + time) %>%
      filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "female") %>%
      mutate(agerelform0 = agerelform - 18)

  } else {

    men <- df %>%
      arrange(ID, relid, FormTime) %>%
      group_by(ID, relid) %>%
      mutate(episodeorder = row_number(),
             agerelform = TOB + FormTime,
             agerelform = first(agerelform),
             pagerelform = agerelform - AgeGap) %>%
      ungroup() %>%
      mutate(age = TOB + time) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "male") %>%
      distinct(ID, relid) %>%
      mutate(agerelform0 = agerelform - 18)

    women <- df %>%
      arrange(ID, relid, FormTime) %>%
      group_by(ID, relid) %>%
      mutate(episodeorder = row_number(),
             agerelform = TOB + FormTime,
             agerelform = first(agerelform),
             pagerelform = agerelform + AgeGap) %>%
      ungroup() %>%
      mutate(age = TOB + time) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "female") %>%
      distinct(ID, relid) %>%
      mutate(agerelform0 = agerelform - 18)

  }

  # Now fit separate models for men and women
  # varPower variance structure for each model
  # The formula to calculate the weights for variance is |v|^(2*t)
  # Age can't be at 0 in varPower formula because then it will
  # evaluate to 0 variance for the first level (18 year olds)
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

  # Add the predicted values to the dataset
  # Level 0 mean population level predictions
  men$pred <- predict(malemodel, men, level = 0)
  women$pred <- predict(femalemodel, women, level = 0)

  # Combine both datasets again
  comb <- bind_rows(men, women)

  # Create list to store model output
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

  # Return datalist with 1. Dataframe based upon input arguments
  # and 2. Dataframe of model output
  agemix.pieces <- list(comb, modoutput)

  return(agemix.pieces)

}
