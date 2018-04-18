#' Model age-mixing pattern.
#'
#' Models the age-mixing pattern for the population simulated in Simpact. The
#' user specifies a time point, time window, and age group for which they would
#' like to obtain a summary of the age-mixing pattern. The user also specifies
#' whether only relationships that started within the time window should be
#' used.
#'
#' The \code{amp.modeller} function can use either the \code{link[nlme]{lme}} or
#' the \code{link[lme4]{lmer}} function to build a linear mixed effects model
#' regressing the partner's age at the time the relationship started on his/her
#' own age at the beginning of the relationship. The models are stratified by
#' gender. Each person can have more than one relationship so a random intercept
#' at the level of the person is added to the model.
#'
#' The models built with \code{link[nlme]{lme}} also explicitly account for
#' heteroskedastic variance, since the variance in partner ages may grow for
#' older ages.
#'
#' The function produces a fitted model.
#'
#' @param dataframe The dataframe that is produced by
#'   \code{\link{agemix.episodes.df.maker}}
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
#' @param SHIMS Should only the three most recent relationships be kept, like in
#'   the SHIMS survey?
#' @param method Should the \code{link[nlme]{lme}} or \code{link[lme4]{lmer}}
#'   function be used? The lme function can be used to model heteroskedastic
#'   residual error.
#'
#' @return returns a model fitted with \code{link[nlme]{lme}} or
#'   \code{link[lme4]{lmer}}
#'
#' @examples
#' data(persreldf)
#' agemixpatdat <- amp.modeller(dataframe = persreldf, agegroup = c(18, 50),
#' timepoint = 30, timewindow = 0.5, start = FALSE, SHIMS = TRUE, method = "lmer")
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom nlme lme
#' @importFrom nlme varPower
#' @importFrom lme4 lmer
#' @export
#'
#@importFrom nlme lme
#@importFrom lme4 lmer
#@importFrom nlme predict
#@importFrom nlme intervals
#@importFrom nlme VarCorr

amp.modeller <- function(dataframe,
                            agegroup,
                            timepoint,
                            timewindow,
                            start = FALSE,
                         SHIMS = TRUE,
                            method = "lmer") {

  agemix.rels.df <- agemix.rels.df.maker(dataframe = dataframe,
                                         agegroup = agegroup,
                                         timepoint = timepoint,
                                         timewindow = timewindow,
                                         start = start)
  if (SHIMS == TRUE){
    agemix.rels.df <- dplyr::group_by(agemix.rels.df, ID) %>%
      dplyr::top_n(-3, FormTime) # male age must be at least 15 at the start of the relationship, but that is always the case in the default Simpact simulations
    agemix.rels.df$agerelform0 <- agemix.rels.df$agerelform - 15
  }
  male.rels.df <- dplyr::filter(agemix.rels.df, Gender =="male")
  female.rels.df <- dplyr::filter(agemix.rels.df, Gender =="female")

  if (method == "lmer") {
    male.model <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                       data = male.rels.df,
                       REML = TRUE)
    female.model <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                       data = female.rels.df,
                       REML = TRUE)
  }
  if (method == "lme") {
    male.model <- lme(pagerelform ~ agerelform0,
                     data = male.rels.df,
                     random = ~1 | ID,
                     weights = nlme::varPower(value = 0.5, form = ~agerelform0 + 1),
                     method = "REML")
    female.model <- lme(pagerelform ~ agerelform0,
                      data = female.rels.df,
                      random = ~1 | ID,
                      weights = nlme::varPower(value = 0.5, form = ~agerelform0 + 1),
                      method = "REML")
  }
  models.list <- list(male.model = male.model, female.model = female.model)
  return(models.list)
}
