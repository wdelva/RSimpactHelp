#' HPC-version of MICE-assisted calibration
#'
#' Produces a csv file of proposed experiments (input parameter values) to match
#' a vector of target features.
#'
#' This version starts from a dataframe of completed experiments and their
#' associated model features.
#'
#' @param targets.empirical The vector of target features
#' @param previous.experiments Dataframe that stores completed experiments
#'   and their associated model features
#' @param RMSD.tol.max Tolerance for the root mean squared distance between
#'   target features and model output features
#' @param min.givetomice Minimal number of observations in the training dataset
#'   to which MICE is applied
#' @param n.experiments Number of proposed experiments to be produced
#' @param strict.positive.params Vector of indices that indicate which of the
#'   input parameters are strictly positive. Set to zero if there are no such
#'   parameters.
#' @param probability.params Vector of indices that indicate which of the input
#'   parameters are strictly between 0 and 1. Set to zero if there are no such
#'   parameters.
#' @param method Method used by MICE. E.g. "norm" or "rf"
#' @param predictorMatrix Can be "complete", "LASSO", or a user-defined matrix
#'   of indices that indicate which variables are included in the chained
#'   equations in MICE
#' @param maxit The maxit argument used in MICE (number of times that the
#'   chained equations are cycled through)
#' @return A dataframe with proposed input parameter values to match a vector of
#'   target features.
#'
#' @importFrom mice mice
#' @importFrom glmnet cv.glmnet
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom magrittr "%>%"
#' @export

HPC_MaC <- function(targets.empirical = dummy.targets.empirical, #  c(2, 5, 0.8, 2, 2, 16)
                    previous.experiments = input.output.df, # = data("/Users/delvaw/Downloads/input.output.RData"),
                    file.name, # e.g. = "/user/data/gent/vsc400/vsc40070/agemixing/scripts/wave2.csv",
                    RMSD.tol.max = 2,
                    min.givetomice = 12,
                    n.experiments = 48,
                    strict.positive.params = c(1,3),
                    probability.params = 0,
                    method = "norm",
                    predictorMatrix = "complete",
                    maxit = 50){
  sim.results.with.design.df <- mutate_all(previous.experiments, function(x) as.numeric(as.character(x))) %>%
    dplyr::filter(complete.cases(.))
  x.offset <- which.max(names(sim.results.with.design.df) %in% "y.1") - 1
  x.names <- paste0("x", seq(1:(x.offset - 1)))
  y.names <- paste0("y.", seq(1:(length(names(sim.results.with.design.df)) - x.offset)))
  # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
  # targets.diff <- targets.empirical - sim.results.with.design.df.median.features # experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
  candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

  # Initiate n.close.to.targets
  n.close.to.targets <- 0 # This will be overwritten.
  RMSD.tol <- 0 # This will be increased if n.close.to.targets < min.givetomice for this tolerance level

  while (n.close.to.targets < min.givetomice & RMSD.tol <= RMSD.tol.max){
    sum.sq.rel.dist <- rep(0, nrow(sim.results.with.design.df))
    for (i in 1:length(targets.empirical)) { # This for loop can be taken out of the while loop, to increase speed.
      name.dist <- paste0("y.", i, ".sq.rel.dist")
      value.dist <- ((sim.results.with.design.df[,i + x.offset] - targets.empirical[i]) / targets.empirical[i])^2
      assign(name.dist, value.dist)
      sum.sq.rel.dist <- sum.sq.rel.dist + get(name.dist)
    }
    RMSD <- sqrt(sum.sq.rel.dist / length(targets.empirical))
    n.close.to.targets <- sum(RMSD <= RMSD.tol, na.rm = TRUE)
    RMSD.tol <- RMSD.tol + 0.01  # Increasing RMSD.tol
  }
  sim.results.with.design.df$RMSD <- RMSD
  # final.intermediate.features <- candidate.intermediate.features

  # 5. Select n.close.to.targets shortest distances
  dist.order <- order(RMSD) # Ordering the squared distances from small to big.
  selected.distances <- dist.order[1:n.close.to.targets]
  sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]

  #calibration.list$new.sim.results.with.design.df[[wave]] <- new.sim.results.with.design.df

  # 5.aaaa Keeping track of medians
  #calibration.list$sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
  # The median of the simulations in the lastest wave
  #calibration.list$new.sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(new.sim.results.with.design.df, contains("y.")))
  # The median of the simulations to give to mice
  #calibration.list$sim.results.with.design.df.selected.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df.selected, contains("y.")))

  # 5.b. Record highest RMSD value for that the selected experiments
  #calibration.list$max.RMSD[[wave]] <- max(sim.results.with.design.df.selected$RMSD)
  # 5.c. Record n.close.target
  #calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets

  # 6. Record selected experiments to give to mice for this wave
  #calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

  # 7. Put intermediate features in dataframe format
  targets.empirical.df <- as.data.frame(matrix(targets.empirical, ncol = length(targets.empirical)))
  names(targets.empirical.df) <- y.names

  # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
  df.give.to.mice <- dplyr::full_join(dplyr::select(sim.results.with.design.df.selected,
                                                    -one_of(c("RMSD", "seed"))), # adding target to training dataset
                                      targets.empirical.df[rep(1:nrow(targets.empirical.df),
                                                               each = n.experiments), ],
                                      by = names(targets.empirical.df)) # "by" statement added to avoid printing message of the variables were used for joining


  #print(df.give.to.mice)
  if (!identical(strict.positive.params, 0)){
    df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
  }
  # probability.params <- 14
  if (!identical(probability.params, 0)){
    df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation
  }

  # 9. Override default predictorMatrix with a sparser matrix
  # Let's think a bit more carefully about which variables should be allowed as input for which input parameters.
  # IN THE FUTURE THIS COULD BE AUTOMATED WITH VARIABLE SELECTION ALGORITHMS.
  # predictorMatrix <- (1 - diag(1, ncol(df.give.to.mice))) # This is the default matrix.

  #### NEW: Using LASSO to create predictorMatrix (and ignoring the one that was given as a function argument)

  if (is.numeric(predictorMatrix)){
    predictorMatrix.give.to.mice <- predictorMatrix
  }

  if (identical(predictorMatrix, "LASSO")){
    predictorMatrix.LASSO <- diag(0, ncol = ncol(df.give.to.mice), nrow = ncol(df.give.to.mice))
    all.names <- names(df.give.to.mice)

    nrows.training.df <- dplyr::select(sim.results.with.design.df.selected,
                                       -one_of(c("RMSD", "seed", "wave"))) %>% nrow()

    for(y.index in 1:ncol(df.give.to.mice)){
      x4lasso <- as.matrix(df.give.to.mice[1:nrows.training.df, -y.index])
      y4lasso <- as.numeric(df.give.to.mice[1:nrows.training.df, y.index])
      alpha <- 1
      cvfit <- glmnet::cv.glmnet(x = x4lasso,
                                 y = y4lasso,
                                 family = "gaussian",
                                 alpha = alpha,
                                 nlambda = 20)
      remaining.indices <- coef(cvfit, s = "lambda.1se")@i
      nonzero.names <- names(df.give.to.mice[-nrow(df.give.to.mice), -y.index])[remaining.indices] # These are the columns with non-zero coefficients
      col.indices <- all.names %in% nonzero.names
      predictorMatrix.LASSO[y.index, col.indices] <- 1
    }
    predictorMatrix.give.to.mice <- predictorMatrix.LASSO
  }

  if (identical(predictorMatrix, "complete")){
    predictorMatrix.give.to.mice <- (1 - diag(1, ncol(df.give.to.mice)))
  }

  # print(c(nrow(df.give.to.mice) - n.experiments, "nrows to give to mice"), quote = FALSE)
  # do imputation
  mice.test <- tryCatch(mice::mice(df.give.to.mice,
                                   m = 1,
                                   method = method,
                                   defaultMethod = method,
                                   predictorMatrix = predictorMatrix.give.to.mice,
                                   maxit = maxit,
                                   printFlag = FALSE),
                        error = function(mice.err) {
                          return(list())
                        })

  # print(c(length(mice.test), "this is length of mice.test", quote = FALSE))
  if (length(mice.test) > 0){

    # 11. Turn mice proposals into a new matrix of experiments

    experiments <- unlist(mice.test$imp) %>% matrix(., byrow = FALSE, ncol = length(x.names)) # %>% data.frame()
    #colnames(mice.guesses3.df) <- imputed.params.names

    # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
    if (!identical(strict.positive.params, 0)){
      experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])
    }
    # And we must also backtransform the logit-transformed values
    if (!identical(probability.params, 0)){
      experiments[, probability.params] <- exp(experiments[, probability.params]) / (1 + exp(experiments[, probability.params]))
    }
    experiments <- as.data.frame(cbind(1:n.experiments, experiments))
    experiments.col.names <- c("seed", x.names)
    names(experiments) <- experiments.col.names
    write.table(experiments,
                sep = ",",
                dec = ".",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE, #experiments.col.names,
                file = file.name)
  }
  return(experiments)
}
