#' A function that determines the facility the client attended
#' Assumption is that the client facility will never change once asssigned
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @return the datalist with facility attached to each client
#'
#' @examples
#' data(datalist)
#' datalist <- client.facility(datalist = datalist, site="All")
#' datalist
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

client.facility <- function(datalist = datalist, site="All"){

  datalist.clients <- datalist$ptable
  datalist.clients$pfacility <- NA
  datalist.clients$pfacility.value <- NA

  if (site == "All") {
    datalist.clients <- datalist.clients
  } else{
    facilities.df <- datalist$ftable
    colnames(facilities.df) <- c("facility.xy", "XCoord", "YCoord")

    for (i in 1:nrow(datalist.clients)) {

      diff.distance <- sqrt((datalist.clients[i, XCoord] - facilities.df$XCoord)^2 +
                              (datalist.clients[i, YCoord] - facilities.df$YCoord)^2 )

      fc.id <- which.min(diff.distance)

      fc.val <- min(diff.distance)

      datalist.clients$pfacility[i] <- facilities.df[fc.id, facility.xy]
      datalist.clients$pfacility.value[i] <- fc.val
    }

  }


  # facilities.df$proportion <- c()
  # datalist.clients <- as.data.frame( datalist.clients %>%
  #   group_by(pfacility) %>%
  #   sample_n(size = 3) %>%
  #   mutate(study.client = "Yes"))
  #
  #
  #
  # datalist.clients[sample(nrow(datalist.clients), 2), ]




  return(datalist.clients)
}

