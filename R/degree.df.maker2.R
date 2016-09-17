#' Degree distribution of recent sexual relationships.
#'
#' \code{degree.df.maker()} creates a table of two columns, "ID" and "degree" from the df dataframe,
#' with specific degrees of new relationships or ongoing relationships among women within
#' specified age group. The degree distribution shows the total number of relationships that
#' were ongoing (if only.new == FALSE) at some time in the survey window
#' or that were newly formed (if only.new == TRUE) during the survey window.
#'
#' @param dataframe.df The dataframe that is produced by \code{\link{agemix.df.maker()}}
#' @param survey.time Time point of the cross-sectional survey.
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param hivstatus HIV status at the time of the survey. Options are >= 0, means all; 0 means only HIV-negative, 1 means only HIV-positive.
#' @param window.width Time period before the survey e.g 1 year before the survey.
#' @param only.new Logical indicator. If TRUE, only relationships that were newly started during window.width are counted
#' (i.e. the individual was NEVER in a relationship with these partners before the start of the window).
#' If FALSE, all relationships that were ongoing at some point during the window.width are counted.
#'
#' @return A dataframe, degree.df that has the columns "ID", "degree", and a fraction of women
#' with >1 sexual partner in the last year.
#' @examples
#' load(dataframe.df)
#' degree.df <- degree.df.maker(dataframe.df, agegroup = c(15, 30), hivstatus = 0, survey.time = 10, window.width = 1, only.new = TRUE)

#' @importFrom magrittr %>%
#' @import dplyr

degree.df.maker <- function(dataframe.df,
                            agegroup = c(15, 30),
                            hivstatus = 0,
                            survey.time = 10,
                            window.width = 1,
                            only.new = TRUE){


  # Before filtering dataframe.df, we need to reshape the dataframe.df,
  # so that there is one row (duplicated because of gender) per relid (relationship)
  # instead of one row (duplicated) per relationship episode.

dataframe.rels.df <- dataframe.df # The dataframe of relationship episodes must be turned into a dataframe of relationships.
dataframe.rels.df <- subset(dataframe.rels.df, select = -c(episodeorder, FormTime, DisTime))
dataframe.rels.df <- unique.data.frame(dataframe.rels.df)

rels.form.dis.df <- dplyr::summarise(group_by(dataframe.df, relid),
                 FormTime = min(FormTime),
                 DisTime = max(DisTime))

dfnew <- left_join(x = dataframe.rels.df,
                   y = rels.form.dis.df,
                   by = "relid")

  # Also only subset the relationships for the people that were still alive at the time of the survey.
dfnew <- subset(dfnew, TOD > survey.time)


  # Now we filter dfnew, based on hivstatus at the time of the survey.

{if (hivstatus==0){
  dfnew <- dplyr::filter(dfnew,
                         InfectTime == "Inf")
}
  else if
  (hivstatus==1){
    dfnew <- dplyr::filter(dfnew,
                           InfectTime!= "Inf")
  }
}

  # newly formed relationships "else" ongoing relationships.
  {if(only.new)
  {dfnew <- dplyr::filter(dfnew,
                                 FormTime >= survey.time-window.width,
                                 FormTime < survey.time,
                                 DisTime > survey.time-window.width,
                                 #InfectTime > survey.time; InfectTime <= survey.time
                                 Gender=='female',survey.time-TOB>=agegroup[1], survey.time-TOB<agegroup[2])}
    else
    {dfnew <- dplyr::filter(dfnew,
                                   FormTime < survey.time,
                                   DisTime > survey.time-window.width,
                                   Gender=='female',survey.time-TOB>=agegroup[1], survey.time-TOB<agegroup[2])}
  }
  # unique relid's(those in the relationship) of each ID.
  uniqueOut <- dfnew %>%
    dplyr::select(ID, relid) %>%
    distinct %>%
    rename(Degree= relid)
  # counting total number of relationships each ID has.
  degreedata.df <- aggregate(Degree ~ ID, data = uniqueOut, length)

  return(degreedata.df)
}
