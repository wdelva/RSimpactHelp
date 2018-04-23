#' Parallelize simpact with phylodynamic component
#'
#' @param model Wrapper function computing epidemiological features
#' @param actual.input.matrix Inputs parameters for simpact setup
#' @param seed_count Initiation of seed for reproducibility
#' @param n_cluster Number of cores
#' @import ape
#' @import treedater
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapplyLB
#' @importFrom parallel stopCluster
#' @export


phylo.simpact.parallel <- function(model = phylo.simpact.wrapper,
                             actual.input.matrix = inputmatrix,
                             #nb_simul = 16,
                             seed_count = 0,
                             n_cluster = 4){
  #cl <- makeCluster(getOption("cl.cores", n_cluster)) # makeCluster: library(snow) - library(parallel)

  cl <- makeCluster(n_cluster)

  #initiate the parameter set
  list_param <- list(NULL)
  tab_simul_summarystat <- tab_param <- paramtemp <- simultemp <-  NULL

  nb_simul <- nrow(actual.input.matrix)

  print("Just to test if this printing ---- YES")

  #the commented variables in  here do not seem to be used.
  for (i in 1:nb_simul) {
    #l <- ncol(actual.input.matrix)
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
    #tab_param <- rbind(tab_param, param[2:(l + 1)])
    #paramtemp <- rbind(paramtemp, param[2:(l + 1)])
  }

  #no need to make a wrapper ftn input this will do
  #source("R/Misc/David/Parallel_Simpact/min.phylo.simpact.wrapper.R")
  #clusterExport(cl, "phylo.simpact.wrapper")

  list_simul_summarystat <- parLapplyLB(cl, list_param,
                                        model)

  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)

  stopCluster(cl)

  return(tab_simul_summarystat)
}

