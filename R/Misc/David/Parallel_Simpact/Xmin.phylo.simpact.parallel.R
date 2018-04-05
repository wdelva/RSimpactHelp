
pacman::p_load(snow, parallel)
# Phylo-Simpact in parallel

simpact.parallel <- function(model = phylo.simpact.wrapper,
                             actual.input.matrix = inputmatrix,
                             #nb_simul = 16,
                             seed_count = 0,
                             n_cluster = 4){
  cl <- makeCluster(getOption("cl.cores", n_cluster)) # makeCluster: library(snow) - library(parallel)
  tab_simul_summarystat = NULL
  list_param <- list(NULL)
  tab_param <- NULL
  paramtemp <- NULL
  simultemp <- NULL

  nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    l <- ncol(actual.input.matrix)
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
    paramtemp <- rbind(paramtemp, param[2:(l + 1)])
  }
  list_simul_summarystat = parLapplyLB(cl, inputvector = list_param,
                                       model)
  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
  stopCluster(cl)
  return(tab_simul_summarystat)
}


source("R/Misc/David/Parallel_Simpact/min.phylo.simpact.wrapper.R")


## Try this


# Phylo-Simpact in parallel


simpact.parallel <- function(model = phylo.simpact.wrapper,
                             actual.input.matrix = inputmatrix,
                             #nb_simul = 16,
                             seed_count = 0,
                             n_cluster = 8){
  cl <- makeCluster(getOption("cl.cores", n_cluster)) # makeCluster: library(snow) - library(parallel)
  tab_simul_summarystat = NULL
  list_param <- list(NULL)
  tab_param <- NULL
  paramtemp <- NULL
  simultemp <- NULL

  nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    l <- ncol(actual.input.matrix)
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
    paramtemp <- rbind(paramtemp, param[2:(l + 1)])
  }
  list_simul_summarystat = parLapplyLB(cl, list_param,
                                       model)
  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
  stopCluster(cl)
  return(tab_simul_summarystat)
}


# Input vector: paratemet combination for Simpact


# inputvector <- c(123, 1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
#                  -0.3, -0.3,
#                  -2.7, # conception
#                  -0.52, -0.05)

inputvector <- c(1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)

reps <- 4


# Input parameters in matrix form reps times (rows).
inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)


features.matrix <- simpact.parallel(model = phylo.simpact.wrapper,
                                    actual.input.matrix = inputmatrix,
                                    seed_count = 123,
                                    n_cluster = 8)
save(features.matrix, file = "features.matrix.csv")
