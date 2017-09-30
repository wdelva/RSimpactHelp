#' Create an age distribution for the population at the start of the simulation.
#'
#' Create an age distribution, consistent with a Weibull survival distribution in
#' the absence of HIV-related mortality, for the population at the start of the simulation.
#'
#' @param shape The shape parameter for the Weibull distribution function (default is 5)
#' @param scale The scale parameter for the Weibull distribution function (default is 65)
#'
#'
#' @return a data.frame with Age, Percent.Male and Percent.Female as variables, and 101
#' rows (0.5 to 100.5 years old). Nobody gets older than 100 years old.
#' @examples
#' agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)
#' agedist.data.frame
#' @importFrom stats pweibull
#' @export
agedistr.creator <- function(shape = 5, scale = 65){
  agebins <- seq(0.5, 99.5, 1)
  probofstillalive <- 1 - pweibull(agebins, shape = shape, scale = scale)
  fractionsinagebins <- 100 * probofstillalive/sum(probofstillalive)
  agedist.data.frame <- data.frame(Age = c(agebins, 100.5),
                                   Percent.Male = c(fractionsinagebins, 0),
                                   Percent.Female = c(fractionsinagebins, 0))
  return(agedist.data.frame)
}
