#' MICE-assisted calibration
#'
#' Produces a list with multiple waves of proposed input parameter values to
#' match a vector of target features.
#'
#' @param targets.empirical The vector of target features
#' @param RMSD.tol.max Tolerance for the root mean squared distance between target features and model output features
#' @param min.givetomice Minimal number of observations in the training dataset to which MICE is applied
#' @param n.experiments Number of parameter combinations in each wave of model runs
#' @param lls Vector of lower limits of the prior distribution of input parameter values
#' @param uls Vector of upper limits of the prior distribution of input parameter values
#' @param model Wrapper function for the simulation model. See details for a description of the required format.
#' @param strict.positive.params Vector of indices that indicate which of the input parameters are strictly positive
#' @param probability.params Vector of indices that indicate which of the input parameters are strictly between 0 and 1
#' @param predictorMatrix Matrix of indices that indicate which variables are included in the chained equations in MICE
#' @param maxit The maxit argument used in MICE (number of times that the chained equations are cycled through)
#' @param maxwaves The maximum number of waves of model runs
#' @param n_cluster The number of cores available for parallel model runs
#' @return A list with multiple waves of proposed input parameter values to
#' match a vector of target features.
#'
#' @import mice
#' @importFrom gsubfn strapplyc
#' @importFrom data.table setDT
#' @importFrom data.table setnames
#' @importFrom readcsvcolumns read.csv.columns
#' @importFrom randtoolbox sobol
#' @importFrom pcaPP l1median
#' @export

MaC <- function(targets.empirical = dummy.targets.empirical,
                            RMSD.tol.max = 2,
                            min.givetomice = 64,
                            n.experiments = 256,
                            lls,
                            uls,
                            model = VEME.wrapper, # simpact.wrapper,
                            strict.positive.params,
                            probability.params,
                            predictorMatrix,
                            maxit = 50,
                            maxwaves = 4,
                            n_cluster = n_cluster){
  # 0. Start the clock
  ptm <- proc.time()
  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored
  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  rel.dist.cutoff <- Inf # initially it is infinitely large, but in later iterations it shrinks
  #sim.results.with.design.df <- NULL # Will be growing with each wave (appending)
  #sim.results.with.design.df.selected <- NULL
  #final.intermediate.features <- NULL

  modelstring <- unlist(paste0(deparse(model), collapse = " "))
  input.vector.length <- max(unique(na.omit(as.numeric(gsubfn::strapplyc(modelstring, "[[](\\d+)[]]", simplify = TRUE))))) - 1 # minus one because the first input parameter is the random seed


  # input.vector.length <- max(unique(na.omit(as.numeric(unlist(strsplit(unlist(paste0(deparse(model), collapse = " ")), "[^0-9]+"))))))
  output.vector.length <- length(targets.empirical)


  x.names <- paste0("x.", seq(1:input.vector.length))
  y.names <- paste0("y.", seq(1:output.vector.length))
  x.offset <- length(x.names)
  sim.results.with.design.df <- data.frame(matrix(vector(), 0, (input.vector.length + output.vector.length),
                                                  dimnames = list(c(), c(x.names, y.names))),
                                           stringsAsFactors = FALSE) # Will be growing with each wave (appending)
  # sim.results.with.design.df.selected I think this does not need to be initiated
  final.intermediate.features <- rep(NA, times = length(targets.empirical))

  # 1. Start loop of waves, based on comparing intermediate features with targets.empirical
  while (wave <= maxwaves & !identical(final.intermediate.features, targets.empirical)){
    print(c("wave", wave), quote = FALSE)

    if (wave == 1){
      # 2. Initial, naive results, based on Sobol sequences
      range.width <- uls - lls
      ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
      range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
      sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
      experiments <- ll.mat + sobol.seq.0.1 * range.width.mat
    }

    sim.results.simple <- simpact.parallel(model = model,
                                           actual.input.matrix = experiments,
                                           seed_count = 0,
                                           n_cluster = n_cluster)
    # save(sim.results.simple, file = "/Users/delvaw/Documents/MiceABC/sim.results.simple.RData")
    # load(file = "sim.results.simple.RData")

    new.sim.results.with.design.df <- as.data.frame(cbind(experiments,
                                                          sim.results.simple))

    names(new.sim.results.with.design.df) <- c(x.names, y.names)

    new.sim.results.with.design.df <- new.sim.results.with.design.df %>% dplyr::filter(complete.cases(.))

    if (wave==1){ # TRUE for the first wave only
      sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                          new.sim.results.with.design.df)
    } else {
      sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                        -contains("RMSD")),
                                          new.sim.results.with.design.df)
    }

    ######## CONTEXT
    # sim.results.with.design.df contains all simulations from previous waves (sim.results.with.design.df)
    # and the simulations from the most recent wave (new.sim.results.with.design.df)
    ########

    # save(sim.results.with.design.df, file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")
    # load(file = "/Users/delvaw/Documents/MiceABC/sim.results.with.design.df.RData")

    # sim.results.with.design.df.median.features <- l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    ######## CONTEXT
    # sim.results.with.design.df.median.features are the median features of ALL experiments for all waves done so far
    ########

    # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
    # targets.diff <- targets.empirical - sim.results.with.design.df.median.features # experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
    candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

    # Initiate n.close.to.targets
    n.close.to.targets <- 0 # This will be overwritten.
    candidate.intermediate.features <- targets.empirical # We start with the empirical target features
    RMSD.tol <- 0 # This will be increased if n.close.to.targets < min.givetomice for this tolerance level

    while (n.close.to.targets < min.givetomice & RMSD.tol <= RMSD.tol.max){
      sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
      for (i in 1:length(candidate.intermediate.features)) { # This for loop can be taken out of the while loop, to increase speed.
        name.dist <- paste0("y.", i, ".sq.rel.dist")
        value.dist <- ((sim.results.with.design.df[,i + x.offset] - candidate.intermediate.features[i]) / candidate.intermediate.features[i])^2
        assign(name.dist, value.dist)
        sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
      }
      RMSD <- sqrt(sum.sq.rel.dist / length(candidate.intermediate.features))
      n.close.to.targets <- sum(RMSD <= RMSD.tol, na.rm = TRUE)
      #n.close.to.targets.mat[(1+steps.intermediate.targets), (1+steps.RMSD.tol)] <- n.close.to.targets
      #large.enough.training.df <- n.close.to.targets >= min.givetomice
      RMSD.tol <- RMSD.tol + 0.01  # Increasing RMSD.tol
    }
    sim.results.with.design.df$RMSD <- RMSD
    final.intermediate.features <- candidate.intermediate.features

    # 5. Select n.close.to.targets shortest distances
    dist.order <- order(RMSD) # Ordering the squared distances from small to big.
    selected.distances <- dist.order[1:n.close.to.targets]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]
    # 5.a. Record intermediate features
    #calibration.list$intermediate.features[[wave]] <- experim.median.features
    # 5.aa. Record final.intermediate.features to be given to mice in that wave
    #calibration.list$final.intermediate.features[[wave]] <- final.intermediate.features
    # 5.aaa. Record experimental parameter values and associated output for that wave (non-complete cases filtered out)
    calibration.list$new.sim.results.with.design.df[[wave]] <- new.sim.results.with.design.df

    # 5.aaaa Keeping track of medians
    # The median across all simulations so far: sim.results.with.design.df.median.features <- l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    calibration.list$sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    # The median of the simulations in the lastest wave
    calibration.list$new.sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(new.sim.results.with.design.df, contains("y.")))
    # The median of the simulations to give to mice
    calibration.list$sim.results.with.design.df.selected.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df.selected, contains("y.")))

    # 5.b. Record highest RMSD value for that the selected experiments
    calibration.list$max.RMSD[[wave]] <- max(sim.results.with.design.df.selected$RMSD)
    # 5.c. Record n.close.target
    calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets

    # 6. Record selected experiments to give to mice for this wave
    calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

    # 7. Put intermediate features in dataframe format
    final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
    names(final.intermediate.features.df) <- y.names

    # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
    df.give.to.mice <- dplyr::full_join(dplyr::select(sim.results.with.design.df.selected,
                                                      -contains("RMSD")), # adding target to training dataset
                                        final.intermediate.features.df[rep(1:nrow(final.intermediate.features.df),
                                                                           each = n.experiments), ],
                                        by = names(final.intermediate.features.df)) # "by" statement added to avoid printing message of the variables were used for joining


    # We are transforming parameters that are necessarily strictly positive: sigma, gamma.a, gamma.b.
    # We could also consider a similar transformation for input parameters that we think should be negative (e.g. formation.hazard.agegapry.gap_factor_man_exp) but for now not yet
    # strict.positive.params <- c(1:4) # Move to the list of arguments of the function
    #print(strict.positive.params)

    #print(df.give.to.mice)
    df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
    # probability.params <- 14
    df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation

    # 9. Override default predictorMatrix with a sparser matrix
    # Let's think a bit more carefully about which variables should be allowed as input for which input parameters.
    # IN THE FUTURE THIS COULD BE AUTOMATED WITH VARIABLE SELECTION ALGORITHMS.
    # predictorMatrix <- (1 - diag(1, ncol(df.give.to.mice))) # This is the default matrix.

    print(c(nrow(df.give.to.mice) - n.experiments, "nrows to give to mice"), quote = FALSE)
    # do imputation
    mice.test <- tryCatch(mice::mice(df.give.to.mice,
                                     m = 1,
                                     method = "norm",
                                     defaultMethod = "norm",
                                     predictorMatrix = predictorMatrix,
                                     maxit = maxit,
                                     printFlag = FALSE),
                          error = function(mice.err) {
                            return(list())
                          })

    print(c(length(mice.test), "this is length of mice.test", quote = FALSE))
    if (length(mice.test) > 0){

      # 11. Turn mice proposals into a new matrix of experiments

      experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names)) # %>% data.frame()
      #colnames(mice.guesses3.df) <- imputed.params.names

      # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
      experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])
      # And we must also backtransform the logit-transformed values
      experiments[, probability.params] <- exp(experiments[, probability.params]) / (1 + exp(experiments[, probability.params]))
      wave <- wave + 1
    } else {
      wave <- maxwaves + 1
    }
  }

  # 13. Check if intermediate features have converged to empirical targests (if not, then maxwaves is reached)
  if (identical(final.intermediate.features, targets.empirical)){
    # 14. Do more "traditional mice optimisation" for the remaining waves
    while (wave <= maxwaves){
      print(c("wave", wave), quote = FALSE)

      sim.results.simple <- simpact.parallel(model = model,
                                             actual.input.matrix = experiments,
                                             seed_count = 0,
                                             n_cluster = n_cluster)

      new.sim.results.with.design.df <- as.data.frame(cbind(experiments,
                                                            sim.results.simple))
      x.names <- paste0("x.", seq(1:ncol(experiments)))
      y.names <- paste0("y.", seq(1:ncol(sim.results.simple)))
      x.offset <- length(x.names)
      names(new.sim.results.with.design.df) <- c(x.names, y.names)

      new.sim.results.with.design.df <- new.sim.results.with.design.df %>% dplyr::filter(complete.cases(.))

      if (nrow(sim.results.with.design.df)==0){ # TRUE for the first wave only
        sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                            new.sim.results.with.design.df)
      } else {
        sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                          -contains("RMSD")),
                                            new.sim.results.with.design.df)
      }


      # Initiate n.close.to.targets
      n.close.to.targets <- 0 # This will be overwritten.
      candidate.intermediate.features <- targets.empirical # We start with the empirical target features
      RMSD.tol <- 0 # This will be increased if n.close.to.targets < min.givetomice for this tolerance level

      while (n.close.to.targets < min.givetomice & RMSD.tol <= RMSD.tol.max){
        sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
        for (i in 1:length(candidate.intermediate.features)) {
          name.dist <- paste0("y.", i, ".sq.rel.dist")
          value.dist <- ((sim.results.with.design.df[,i + x.offset] - candidate.intermediate.features[i]) / candidate.intermediate.features[i])^2
          assign(name.dist, value.dist)
          sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
        }
        RMSD <- sqrt(sum.sq.rel.dist / length(candidate.intermediate.features))
        n.close.to.targets <- sum(RMSD <= RMSD.tol, na.rm = TRUE)
        #n.close.to.targets.mat[(1+steps.intermediate.targets), (1+steps.RMSD.tol)] <- n.close.to.targets
        #large.enough.training.df <- n.close.to.targets >= min.givetomice
        RMSD.tol <- RMSD.tol + 0.01  # Increasing RMSD.tol
      }
      sim.results.with.design.df$RMSD <- RMSD
      final.intermediate.features <- candidate.intermediate.features

      # 5. Select n.close.to.targets shortest distances
      dist.order <- order(RMSD) # Ordering the squared distances from small to big.
      selected.distances <- dist.order[1:n.close.to.targets]
      sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]

      # 5.a. Record intermediate features
      # calibration.list$intermediate.features[[wave]] <- experim.median.features
      # 5.aa. Record final.intermediate.features to be given to mice in that wave
      # calibration.list$final.intermediate.features[[wave]] <- final.intermediate.features
      # 5.aaa. Record experimental parameter values and associated output for that wave (non-complete cases filtered out)
      calibration.list$new.sim.results.with.design.df[[wave]] <- new.sim.results.with.design.df

      # 5.aaaa Keeping track of medians
      # The median across all simulations so far: sim.results.with.design.df.median.features <- l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
      calibration.list$sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
      # The median of the simulations in the lastest wave
      calibration.list$new.sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(new.sim.results.with.design.df, contains("y.")))
      # The median of the simulations to give to mice
      calibration.list$sim.results.with.design.df.selected.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df.selected, contains("y.")))



      # 5.b. Record highest RMSD value for that the selected experiments
      calibration.list$max.RMSD[[wave]] <- max(sim.results.with.design.df.selected$RMSD)
      # 5.c. Record n.close.target
      calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets

      # 6. Record selected experiments to give to mice for this wave
      calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

      # 7. Put intermediate features in dataframe format
      final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
      names(final.intermediate.features.df) <- y.names

      # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
      df.give.to.mice <- dplyr::full_join(dplyr::select(sim.results.with.design.df.selected,
                                                        -contains("RMSD")), # adding target to training dataset
                                          final.intermediate.features.df[rep(1:nrow(final.intermediate.features.df),
                                                                             each = n.experiments), ],
                                          by = names(final.intermediate.features.df)) # "by" statement added to avoid printing message of the variables were used for joining

      # We are transforming parameters that are necessarily strictly positive: sigma, gamma.a, gamma.b.
      # We could also consider a similar transformation for input parameters that we think should be negative (e.g. formation.hazard.agegapry.gap_factor_man_exp) but for now not yet
      # strict.positive.params <- c(1:4)
      df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
      df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation
      # 9. Override default predictorMatrix with a sparser matrix
      # Let's think a bit more carefully about which variables should be allowed as input for which input parameters.
      # IN THE FUTURE THIS COULD BE AUTOMATED WITH VARIABLE SELECTION ALGORITHMS.
      # predictorMatrix <- (1 - diag(1, ncol(df.give.to.mice))) # This is the default matrix.

      print(c(nrow(df.give.to.mice) - n.experiments, "nrows to give to mice"), quote = FALSE)

      mice.test <- tryCatch(mice::mice(df.give.to.mice,
                                       m = 1,
                                       method = "norm",
                                       defaultMethod = "norm",
                                       predictorMatrix = predictorMatrix,
                                       maxit = maxit,
                                       printFlag = FALSE),
                            error = function(mice.err) {
                              return(list())
                            })

      if (length(mice.test) > 0){

        # 11. Turn mice proposals into a new matrix of experiments
        experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names)) # %>% data.frame()
        # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
        experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])
        # And we must also backtransform the logit-transformed values
        experiments[, probability.params] <- exp(experiments[, probability.params]) / (1 + exp(experiments[, probability.params]))
        wave <- wave + 1
      } else {
        wave <- maxwaves + 1
      }
    }
  }

  # 15. Stop clock and return calibration list
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock
  return(calibration.list)
}
