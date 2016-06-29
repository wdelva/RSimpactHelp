#' Create dataframe to use for age-mixing calculations.
#'
#' Create a useable dataframe from a simulated population that can then be used for age-mixing metrics.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @return A dataframe of class..
#' @examples
#' agemix.df <- agedif.df.maker(datalist = datalist)


agemix.df.maker <- function(datalist) {

  # Create dataframe that merges all the necessary variables from
  # different tables
  # Each row represents a relationship episode
  # Each relationship could have multiple rows if a relationship
  # The same relationship episode will also be represented twice
  # —once by the female partner and once by the male partner

  if(!is.list(datalist) | length(datalist) != 5) {
    stop("Datalist wrong type")
  }

  dfrmale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDm) %>%
    mutate(relid = paste0(ID, IDw)) # Unique id for the relationship

  dfrfemale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDw) %>%
    mutate(relid = paste0(IDm, ID)) # Male id first, as above

  dfmale <- datalist$ptable %>%
    data.frame() %>%
    select(ID, Gender, TOB) %>%
    filter(Gender == 0) %>%
    left_join(dfrmale, by = "ID") %>%
    select(-IDw)

  dffemale <- datalist$ptable %>%
    data.frame() %>%
    select(ID, Gender, TOB) %>%
    filter(Gender == 1) %>%
    left_join(dfrfemale, by = "ID") %>%
    select(-IDm)

  df <- bind_rows(dfmale, dffemale) %>%
    mutate(Gender = factor(Gender, labels = c("male", "female")))

  return(df)
}
