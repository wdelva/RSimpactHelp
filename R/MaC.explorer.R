#' Distribute simpact runs from MaC output over multiple cores
#'
#' A short description here...
#'
#' @param model Wrapper function for running simpact (simpact.wrapper)
#' @param MaC.output Output from \link{MaC} function.
#' @param n_cluster Number of cores available for parallel running of Simpact
#' @return a matrix of model features and the seed of the random number
#'   generator
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapplyLB
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export

MaC.explorer <- function(model = model,
                         MaC.output = test.MaC,
                         n_cluster = n_cluster){
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  tab_simul_summarystat = NULL
  list_param <- list(NULL)
  tab_param <- NULL
  paramtemp <- NULL
  simultemp <- NULL

  last.wave <- length(test.MaC$selected.experiments)
  experiments <- dplyr::select(test.MaC$selected.experiments[[last.wave]],
                   contains("x."))
  seeds <- dplyr::select(test.MaC$selected.experiments[[last.wave]],
                                          contains("seed")) %>%
    unlist() %>% as.numeric()
  nb_simul <- nrow(experiments)

  for (i in 1:nb_simul) {
    l <- ncol(experiments)
    param <- c(seeds[i], experiments[i, ])
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
    paramtemp <- rbind(paramtemp, param[2:(l + 1)])
  }
  list_simul_summarystat = parLapplyLB(cl, list_param,
                                       model)
  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
  stopCluster(cl)
  return(cbind(tab_simul_summarystat, seeds))
}

