#' Create dataframe for age-mixing calculations.
#'
#' The dataframe that is produced keeps columns from the relationship table
#' and merges them with person-level attributes. A unique id relid for each
#' relationship is generated that combines the unique ID for the male partner
#' first and the female partner second in the relationship.
#'
#' The resulting dataframe is in a long format where each row/observation
#' represents an episode from each relationship. Each relationship episode is
#' represented twice: once for each man in the dataset who was the male partner
#' in that relationship, and once for each woman who was the female partner in
#' the relationship. This allows the user to conduct stratified analyses based
#' on gender. For relationships where there is more than one episode, there will
#' be multiple rows. The order in which the episodes occur are numbered in
#' ascending order with the column episodeorder.
#'
#' Additional columns generated:
#' agerelform - The age of the person when when the relationship began. In
#'       other words, the age at the first episode
#' pagerelform - The age of their partner when the relationship began.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @return A dataframe
#' @examples
#' data(datalist)
#' agemix.df <- agemix.df.maker(datalist = datalist)
#'
#' @importFrom  magrittr %>%
#' @import dplyr
#' @export

agemix.df.maker <- function(datalist) {

  if(!is.list(datalist) |
     names(datalist)[[1]] != "ptable" |
     names(datalist)[[2]] != "rtable") {
    stop("Datalist wrong type")
  }

  dfrmale <- datalist$rtable %>%
    data.frame() %>%
    dplyr::rename(ID = ID1) %>%
    dplyr::mutate(relid = paste0(ID, ID2))

  dfrfemale <- datalist$rtable %>%
    data.frame() %>%
    dplyr::rename(ID = ID2) %>%
    dplyr::mutate(relid = paste0(ID1, ID))

  dfmale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 0) %>%
    dplyr::left_join(dfrmale, by = "ID") %>%
    dplyr::select(-ID2)

  dffemale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 1) %>%
    dplyr::left_join(dfrfemale, by = "ID") %>%
    dplyr::select(-ID1)

  df <- dplyr::bind_rows(dfmale, dffemale) %>%
    dplyr::arrange(Gender, ID, relid, FormTime) %>%
    dplyr::group_by(Gender, ID, relid) %>%
    dplyr::mutate(episodeorder = row_number(),
           agerelform = FormTime - TOB,
           agerelform = dplyr::first(agerelform)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Gender = factor(Gender, labels = c("male", "female")),
           pagerelform = ifelse(Gender == "male",
                                agerelform - AgeGap,
                                agerelform + AgeGap))

  return(df)
}
