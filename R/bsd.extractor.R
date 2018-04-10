#' Compute between-subject standard deviation from agemixing model.
#'
#' Returns the between-subject standard deviation from the agemixing model that
#' was fitted with \code{amp.modeller}.
#'
#'@param models The models object from which the between-subject standard deviation must
#'  be extracted.
#'@return The between-subject standard deviation
#'@export
#' @examples
#'agemix.episodes.df <- agemix.episodes.df.maker(datalist.agemix)
#'amp.models <- amp.modeller(dataframe = agemix.episodes.df, agegroup = c(15, 30),
#'                           timepoint = 15, timewindow = 3, start = FALSE, method = "lme4")
#'

#' @importFrom nlme VarCorr

bsd.extractor <- function(models, gender = "male"){
  model.index <- ifelse(gender == "male", 1, 2)
  if (attr(models[[model.index]], "class") == "lme") {
    bsd <- as.numeric(nlme::VarCorr(models[[model.index]])[1,2])
  } else {
    bsd <- as.numeric(as.data.frame(lme4::VarCorr(models[[model.index]]))[1,5])
  }
  return(bsd)
}
