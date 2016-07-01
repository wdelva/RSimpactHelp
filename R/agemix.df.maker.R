#' Create dataframe for age-mixing calculations.
#'
#' \code{agemix.df.maker()} creates a useable dataframe from Simpact- simulated
#' tables that can then be used for age-mixing calculations.
#'
#' The dataframe that is produced keeps columns from the relationship table
#' and merges them with person-level attributes. A unique id —"relid"— for each
#' relationship is generated that combines the unique ID for the male partner
#' (first) and the female partner (second) in the relationship.
#'
#' The resulting dataframe is in a long format where each row/observation
#' represents an episode from each relationship. Each relationship episode is
#' represented twice: once for each man in the dataset who was the male partner
#' in that relationship, and once for each woman who was the female partner in
#' the relationship. This allows the user to conduct stratified analyses based
#' on gender. For relationships where there is more than one episode, there will
#' be multiple rows. The order in which the episodes occur are numbered in
#' ascending order with the column "episodeorder".
#'
#' @section Additional columns generated:
#'
#'   \describe{
#'     \item{agerelform} {The age of the person when when the relationship began. In
#'       other words, the age at the first episode}
#'     \item{pagerelform} {The age of their partner when the relationship began.}
#'    }
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}.
#' @return A dataframe of class "tbl_df" "tbl" "data.frame"
#' @examples
#' load(datalist)
#' agemix.df <- agemix.df.maker(datalist = datalist)
#'
#' @importFrom magrittr %>%
#' @import dplyr


agemix.df.maker <- function(datalist) {

  if(!is.list(datalist) |
     names(datalist)[[1]] != "ptable" |
     names(datalist)[[2]] != "rtable") {
    stop("Datalist wrong type")
  }

  dfrmale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDm) %>%
    mutate(relid = paste0(ID, IDw))

  dfrfemale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDw) %>%
    mutate(relid = paste0(IDm, ID))

  dfmale <- datalist$ptable %>%
    data.frame() %>%
    filter(Gender == 0) %>%
    left_join(dfrmale, by = "ID") %>%
    select(-IDw)

  dffemale <- datalist$ptable %>%
    data.frame() %>%
    filter(Gender == 1) %>%
    left_join(dfrfemale, by = "ID") %>%
    select(-IDm)

  df <- bind_rows(dfmale, dffemale) %>%
    arrange(Gender, ID, relid, FormTime) %>%
    group_by(Gender, ID, relid) %>%
    mutate(episodeorder = row_number(),
           agerelform = FormTime - TOB,
           agerelform = first(agerelform)) %>%
    ungroup() %>%
    mutate(Gender = factor(Gender, labels = c("male", "female")),
           pagerelform = ifelse(Gender == "male",
                                agerelform - AgeGap,
                                agerelform + AgeGap))

  return(df)
}
