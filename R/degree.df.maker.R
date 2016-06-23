degree.df.maker <- function(datalist = datalist, survey.time = 10, window.width = 1, only.new = TRUE){
  # For the people in the datalist,
  # make a dataframe that has the variables of the ptable,
  # but in addition, also the number of relationships that were ongoing (if only.new == FALSE)
  # or that were newly formed (if only.new == TRUE) at some time in the survey window.
  # (E.g. The one-year period from 9 to 10 years into the simulation).

  return(degree.df)
}
