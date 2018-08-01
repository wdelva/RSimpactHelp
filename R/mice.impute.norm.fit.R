#' adjusted mice norm fitter
#'
#' Augmented to also return details of the regression model and residual error sampled.
#'
#' @param y y
#' @param ry ry
#' @param x x
#' @param wy wy
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice
#' @export

mice.impute.norm.fit <- function(y, ry, x, wy = NULL, ...)
{
  if (is.null(wy))
    wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw(y, ry, x, ...)
  returned.list <- list()
  rnorm.values <- rnorm(sum(wy))
  returned.list$imputation <- x[wy, ] %*% parm$beta + rnorm.values * parm$sigma
  returned.list$parm <- parm
  returned.list$rnorm.values <- rnorm.values
  return(returned.list)
}
