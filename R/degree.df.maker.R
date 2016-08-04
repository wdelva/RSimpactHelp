#' Degree distribution of recent sexual relationships.
#'
#'  \code{degree.df.maker()} creates a table of two columns, "ID" and "degree" from the df dataframe,
#'  with specific degrees of new relationships and ongoing relationships.
#'
#'  The degree distribution shows the total number of relationships that
#'  were ongoing (if only.new == FALSE) or that were newly formed (if only.new == TRUE)
#'  at some time in the survey window.
#'
#'  @param dataframe.df The dataframe that is produced by \code{\link{agemix.df.maker()}}
#'  @param survey.time period of the survey or simulation.
#'  @param window.width period before the end of the survey e.g 1 year before the end of the survey.
#'  @param only.new relationships that were newly started during window.width (i.e. the individual
#'  was NEVER in a relationship with these partners before the start of the window).
#'
#'  @return A dataframe, degree.df that has the columns "ID", "degree".
#'  @examples
#'  load(dataframe.df)
#'  degree.df <- degree.df.maker(dataframe = dataframe.df, survey.time = 10, window.width = 1, only.new = TRUE)
#'  @importFrom magrittr %>%
#'  @import dplyr

degree.df.maker <- function(dataframe.df,
                            survey.time = 20,
                            window.width = 1,
                            only.new = TRUE){

# newly formed relationships "else" ongoing relationships.
  {if(only.new)
  {dataframe.df <- dplyr::filter(dataframe.df, FormTime >=survey.time-window.width, DisTime >= survey.time)}
    else
    {dataframe.df <- dplyr::filter(dataframe.df, FormTime < survey.time-window.width, DisTime >= survey.time)}
  }
# unique relid's(those in the relationship) of each ID.
  uniqueOut <- dataframe.df %>%
    dplyr::select(ID, relid) %>%
    distinct %>%
    rename(Degree= relid)
# counting total number of relationships each ID has.
  degree.df <- aggregate(Degree ~ ID, data = uniqueOut, length)
  return(degree.df)
}
