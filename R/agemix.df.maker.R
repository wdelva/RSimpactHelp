#' Create dataframe for relationship analyses.
#'
#' Produces a dataframe of relationships from list created by \code{\link{readthedata}}. 
#' 
#' The resulting dataframe is in a long format where each row/observation 
#' represents an episode from each relationship. Each relationship episode is 
#' represented twice: once for each man in the dataset who was the male partner 
#' in that relationship, and once for each woman who was the female partner in 
#' the relationship. This allows the user to conduct stratified analyses based 
#' on gender. For relationships where there is more than one episode, there will
#' be multiple rows. The order in which the episodes occur are numbered in 
#' ascending order with the column episodeorder.

#' agerelform - The age of their partner when the relationship began.
#'
#' @param datalist The list object that is produced by \code{\link{readthedata}}
#' @return A dataframe, with the following additional columns generated:
#' \describe{
#'     \item{relid}{Unique relationship identier. It combines the unique ID for
#'     the male parter first, and the female partner second in the relationship.}
#'     \item{agerelform}{The age of the person when when the relationship began. In
#'       other words, the age at the first episode} 
#'  }
#' @examples
#' cfg <- list()
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp")
#' dl <- readthedata(modeloutput)
#' agemix.df <- agemix.df.maker(datalist = dl)
#'
#' @importFrom  magrittr %>%
#' @import dplyr
#' @import tibble
#' @export

agemix.df.maker <- function(datalist) {
  
  # Check to make sure the list has ptable and rtable in it
  if(!is.list(datalist) |
     names(datalist)[[1]] != "ptable" |
     names(datalist)[[2]] != "rtable") {
    stop("Datalist wrong type")
  }

  # Create relationship table for men
  dfrmale <- datalist$rtable %>% # Extract relationship table
    data.frame() %>% # Format as tibble
    dplyr::rename(ID = ID1) %>% # Rename the male ID variable
    dplyr::mutate(relid = paste0(ID, ID2)) # Create unique rel id
  
  # Create relationship table for women
  dfrfemale <- datalist$rtable %>% 
    data.frame() %>%
    dplyr::rename(ID = ID2) %>% # Rename the female ID variable
    dplyr::mutate(relid = paste0(ID1, ID))

  # Create person-level table for men
  dfmale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 0) %>% # Keep only men
    dplyr::left_join(dfrmale, by = "ID") %>% # Merge person to relationship data
    dplyr::select(-ID2) # Remove the female ID variable
  
  # Create person-level table for women
  dffemale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 1) %>% # Keep only women
    dplyr::left_join(dfrfemale, by = "ID") %>% # Merge person to relationship data
    dplyr::select(-ID1) # Remove the male ID variable

  df <- dplyr::bind_rows(dfmale, dffemale) %>% # Combine male and female relationships
    dplyr::arrange(Gender, ID, relid, FormTime) %>% # Sort by timing of episode
    dplyr::group_by(Gender, ID, relid) %>% # By gender, person, and relationship
    dplyr::mutate(episodeorder = row_number(), # Order of episodes
           agerepisform = FormTime - TOB, # Age at start of episode
           agerelform = dplyr::first(agerepisform)) %>% # Age at first episode
    dplyr::ungroup() %>%
    dplyr::mutate(Gender = factor(Gender, labels = c("male", "female")),
           pagerelform = ifelse(Gender == "male",
                                agerelform - AgeGap,
                                agerelform + AgeGap))

  return(df)
}
