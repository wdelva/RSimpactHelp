#' Plot the histogram of a degree distribution.
#'
#' Plot the histogram of the degree distribution that is stored in degreedata.df$Degree, after running \code{\link{degree.df.maker()}}
#'
#' @param dataframe The dataframe that is produced by \code{\link{degree.df.maker()}}
#' @return a ggplot2 object
#' @examples
#' degree.hist.maker <- degree.hist.maker(dataframe = degreedata.df)
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2

degree.hist.maker <- function(dataframe = degreedata.df){
  fig <- dataframe %>%
    ggplot(aes(x=Degree)) +
    geom_histogram(binwidth = 1, aes(y=..density..)) +
    xlab("Degree") +
    ylab("Density") +
    #facet_wrap(~Gender, nrow=1)
    theme_bw()
  return(fig)
}
