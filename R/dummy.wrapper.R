#' Fast and simple wrapper function to test MICE-assisted calibration functions.
#'
#' A made-up model with no meaning whatsoever.
#'
#' @param inputvector Vector of 4 strictly positive input parameter
#'   values, e.g. c(1000, 2, 3, 4)
#' @return A vector of model features
#' @examples
#' input.vector <- c(1000, 2, 3, 4)
#' dummy.features <- dummy.wrapper(input.vector)
#' @export


dummy.wrapper <- function(inputvector = input.vector){ # Let's start with 4 input parameters, e.g. 1000, 2, 3, 4
  set.seed(inputvector[1])
  random.gamma.data <- rgamma(n = ceiling(inputvector[2]), shape = inputvector[3], scale = inputvector[4])
  random.norm.data <- rnorm(n = ceiling(inputvector[2]), mean = inputvector[3]*inputvector[4], sd = inputvector[5])
  feature1 <- median(random.gamma.data)
  feature2 <- sd(random.gamma.data * random.norm.data)
  feature3 <- IQR(random.gamma.data + random.norm.data)
  feature4 <- diff(range(random.gamma.data * random.norm.data))
  outputvector <- c(feature1, feature2, feature3, feature4)
  return(outputvector)
}
