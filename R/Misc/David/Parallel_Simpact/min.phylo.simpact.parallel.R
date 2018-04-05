
pacman::p_load(snow, parallel)

source("R/Misc/David/Parallel_Simpact/min.phylo.simpact.wrapper.R")
# Phylo-Simpact in parallel
simpact.parallel <- function(model = phylo.simpact.wrapper,
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
                                        phylo.simpact.wrapper)

  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)

  stopCluster(cl)

  return(tab_simul_summarystat)
}


# Input vector: paratemet combination for Simpact
inputvector <- c(1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)
#
# inputvector <- c(151,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
#                  -0.3, -0.3,
#                  -2.7, # conception
#                  -0.52, -0.05)
#
reps <- 240

# Input parameters in matrix form reps times (rows).
inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)

#the seed can be pre-pended like this without a loop
#inputmatrix <- cbind(c(1:reps),inputmatrix)

sim.start.time <- proc.time()[3]

features.matrix <- simpact.parallel(model = phylo.simpact.wrapper,
                                    actual.input.matrix = inputmatrix,
                                    seed_count = 123,
                                    n_cluster = 4)

sim.end.time <- proc.time()[3] - sim.start.time
print(paste0("Simulation time: ", round(sim.end.time/60,2), " minutes"))

save(features.matrix, file = "features.matrix.csv")

