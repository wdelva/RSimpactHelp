#' Run mice in parallel
#'
#' A short description here...
#'
#' @param mice.model Wrapper function for mice (mice.wrapper)
#' @param df.give.to.mice Training dataset to give to mice
#' @param m Number of imputations
#' @param predictorMatrix Predictor matrix
#' @param method Method for imputations
#' @param defaultMethod Default method for imputations
#' @param maxit Number of times we cycle through chained equations
#' @param printFlag print flag TRUE of FALSE
#' @param seed_count Origin of random number seed
#' @param n_cluster Number of cores available for parallel running of mice
#' @param nb_simul Number of guesses to be produced by mice
#' @return a matrix of mice guesses
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapplyLB
#' @importFrom magrittr %>%
#' @export

mice.parallel <- function(mice.model = mice.wrapper,
                          df.give.to.mice = df.give.to.mice,
                          m = 1,
                          predictorMatrix = predictorMatrix,
                          method = "norm",
                          defaultMethod = "norm",
                          maxit = 5,
                          printFlag = FALSE,
                          seed_count = 0,
                          n_cluster = 8,
                          nb_simul = n.experiments){
  cl <- makeCluster(getOption("cl.cores", n_cluster))


  mice.input.lists <- vector("list", length = nb_simul) #list(NULL)

  mice.input.list <- list() # An element of the mice.input.lists object. Used as input for mice
  mice.input.list$data <- df.give.to.mice # Not dependent on i
  mice.input.list$m <- m # number of imputations
  mice.input.list$predictorMatrix <- predictorMatrix
  mice.input.list$method <- method
  mice.input.list$defaultMethod <- defaultMethod
  mice.input.list$maxit <- maxit
  mice.input.list$printFlag <- printFlag

  tab_simul_guess = NULL

  # nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    mice.input.list$seed <- seed_count + i
    mice.input.lists[[i]] <- mice.input.list
  }
  list_simul_guesses = parLapplyLB(cl, mice.input.lists,
                                       mice.model)
  tab_simul_guess <- do.call(rbind, list_simul_guesses)
  stopCluster(cl)
  return(tab_simul_guess)
}

