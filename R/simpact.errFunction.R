#' Handling Simpact errors.
#'
#' \code{simpact.errFunction} returns an empty list if the error message of a
#' problematic Simpact run contains \code{NA}.
#'
#' @param e The error object returned by the Simpact run.
#' @return An empty list.
#' @export
#' @examples
#'cfg <- list()
#'results <- tryCatch(simpact.run(configParams = cfg, destDir = "temp"),
#'                    error = simpact.errFunction)

simpact.errFunction <- function(e){
  if (length(grep("NaN",e$message)) != 0){
    return(list())
  }
}
