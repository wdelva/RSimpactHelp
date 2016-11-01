#' This function will control the errors that can be produced by simpact run.
#' Error function catches as well the case when population.maxevents is reached before population.simtime is reached
#' @param e is the function that runs simpact.run
#' @return a string indicating that the program terminated premature
#' @examples
#' out.test <- tryCatch(simpact.run(), error = err.function)

err.function <- function(e){
  if (length(grep("MAXEVENTS",e$message)) != 0)
    return(chunk.summary.stats = rep(NA,length(select.summary.params()[[1]])))
  if (length(grep("internal event time",e$message)) != 0)
    return(chunk.summary.stats = rep(NA,length(select.summary.params()[[1]])))
  stop(e)
}


