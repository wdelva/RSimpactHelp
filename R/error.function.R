#' We expect this to grow and the returned summary statistics to include a lot more.
#' Will include an Ignore variable to each of the summary statistics. Only computed if IGNORE is set to FALSE
#'
#' After the simpact simulation a number of statistics can be computed to summarise the output
#'
#' art.intro list should include the time when the diagnosis to take effect.
#'
#'
#' The function take in a list of paramenters set at intervals of cd4 count threshhold used to
#' decide if a person can be offerd treatment
#'
#' @param the output message from simpact.run \code{\link{readthedata()}}
#' @return a string indicating that the program terminated premature
#' @examples
#' out.test <- tryCatch(simpact.run(), error = errFunction)

# Creating an error function to catch the case when population.maxevents is reached before population.simtime is reached
errFunction <- function(e){
  if (length(grep("MAXEVENTS",e$message)) != 0)
    return("Not Complete")
  if (length(grep("internal event time",e$message)) != 0)
    return("Not Complete")
  stop(e)
}


