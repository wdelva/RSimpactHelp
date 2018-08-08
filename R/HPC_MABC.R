#' HPC-version of MABC (weighted MICE ABC)
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
#' @param file.name Path to csv file that will store the proposed experiments
#' @param RMSD.tol.max Tolerance for the root mean squared distance between
#'   target features and model output features
#' @param min.givetomice Minimal number of observations in the training dataset
#'   to which MICE is applied
#' @param n.experiments Number of proposed experiments to be produced
#' @param lls Vector of lower limits of the prior distribution of input
#'   parameter values
#' @param uls Vector of upper limits of the prior distribution of input
#'   parameter values
#' @param strict.positive.params Vector of indices that indicate which of the
#'   input parameters are strictly positive. Set to zero if there are no such
#'   parameters.
#' @param probability.params Vector of indices that indicate which of the input
#'   parameters are strictly between 0 and 1. Set to zero if there are no such
#'   parameters.
#' @param inside_prior TRUE by default. If FALSE, parameter sampling is not
#'   restricted to the initial ranges of the prior distribution during the
#'   subsequent algorithm steps.
#' @param method Method used by MICE. E.g. "norm" or "rf"
#' @param predictorMatrix Can be "complete", "LASSO", or a user-defined matrix
#'   of indices that indicate which variables are included in the chained
#'   equations in MICE
#' @param maxit The maxit argument used in MICE (number of times that the
#'   chained equations are cycled through)
#' @return A dataframe with proposed input parameter values to match a vector of
#'   target features.
#' @examples
#' features.pop.growth <- exp(0.015)
#' features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
#' features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
#' features.art.cov <- c(0.33, 0.38, 0.45, 0.51, 0.61, 0.7, 0.8)
#' features.vl.suppr <- 0.68
#' target.features.EAAA <- c(features.pop.growth, features.hiv.prev, features.hiv.inc, features.art.cov, features.vl.suppr)
#' priors.EAAA <-list(c("unif", 1.0, 1.5), # hivtransmission.param.f1 = 1.1
#' c("unif", 0.0, 0.5), # formation.hazard.agegapry.gap_agescale_man and ~_woman = 0.25
#' c("unif", -3.0, 3.0), # person.agegap.man.dist.normal.mu and ~.woman.~ = 0
#' c("unif", 1.0, 5.0), # person.agegap.man.dist.normal.sigma and ~.woman.~ = 2.5
#' c("unif", 0.1, 1.5), # person.eagerness.man.dist.gamma.a = 0.23
#' c("unif", 0.1, 1.5), # person.eagerness.woman.dist.gamma.a = 0.23
#' c("unif", 10.0, 80), # person.eagerness.man.dist.gamma.b = 45
#' c("unif", 10.0, 80), # person.eagerness.woman.dist.gamma.b = 45
#' c("unif", -2, -0.2), # formation.hazard.agegapry.gap_factor_man_exp and ~_woman_~ = -0.7
#' c("unif", -1.0, 5.0), # formation.hazard.agegapry.baseline = 2.8
#' c("unif", -1, -0.05), # formation.hazard.agegapry.numrel_man = -0.5
#' c("unif", -1, -0.05), # formation.hazard.agegapry.numrel_woman = -0.5
#' c("unif", -5, -1.0), # conception.alpha_base = -2.7
#' c("unif", -2, 1)) # dissolution.alpha_0 = -0.52
#' prior.boundaries.booleans <- !unlist(priors.EAAA) %in% "unif"
#' boundaries.matrix <- unlist(priors.EAAA)[prior.boundaries.booleans] %>% as.numeric() %>% matrix(byrow = TRUE, ncol = 2)
#' nextwave.df <- HPC_MABC(targets.empirical = target.features.EAAA,
#' previous.experiments = data("input.output.EAAA.RData"),
#' file.name = "/temp/wave2.csv",
#' RMSD.tol.max = 2,
#' min.givetomice = 40,
#' n.experiments = 800,
#' lls = boundaries.matrix[, 1],
#' uls = boundaries.matrix[, 2],
#' strict.positive.params = 0,
#' probability.params = 0,
#' inside_prior = TRUE,
#' method = "norm",
#' predictorMatrix = "LASSO",
#' maxit = 20)
#'
#' @import mice
#' @importFrom glmnet cv.glmnet
#' @importFrom dplyr mutate_all
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom magrittr "%>%"
#' @importFrom gtools smartbind
#' @export

HPC_MABC <- function(targets.empirical = dummy.targets.empirical,
                     previous.experiments,
                     file.name,
                     RMSD.tol.max = 2,
                     min.givetomice = 12,
                     n.experiments = 48,
                     lls,
                     uls,
                     strict.positive.params = 0,
                     probability.params = 0,
                     inside_prior = TRUE,
                     method = "norm",
                     predictorMatrix = "complete",
                     maxit = 20){
  sim.results.with.design.df <- mutate_all(previous.experiments, function(x) as.numeric(as.character(x))) %>%
    dplyr::filter(complete.cases(.))
  x.offset <- which.max(names(sim.results.with.design.df) %in% "y.1") - 1
  x.names <- paste0("x", seq(1:length(lls)))
  y.names <- paste0("y.", seq(1:(length(names(sim.results.with.design.df)) - x.offset)))
  # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
  # targets.diff <- targets.empirical - sim.results.with.design.df.median.features # experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
  candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

  # Initiate n.close.to.targets
  n.close.to.targets <- 0 # This will be overwritten.
  RMSD.tol <- 0 # This will be increased if n.close.to.targets < min.givetomice for this tolerance level

  while (n.close.to.targets < min.givetomice & RMSD.tol <= RMSD.tol.max){
    diff.matrix <- sweep(x = sim.results.with.design.df[ , ((1 + x.offset):(x.offset + length(targets.empirical)))], MARGIN = 2, targets.empirical)
    rel.diff.matrix <- sweep(diff.matrix, MARGIN = 2, targets.empirical, FUN = "/")
    squared.rel.diff.matrix <- rel.diff.matrix^2
    sum.squared.rel.diff <- rowSums(squared.rel.diff.matrix)
    RMSD <- sqrt(sum.squared.rel.diff / length(targets.empirical))
    n.close.to.targets <- sum(RMSD <= RMSD.tol, na.rm = TRUE)
    RMSD.tol <- RMSD.tol + 0.01  # Increasing RMSD.tol
  }
  sim.results.with.design.df$RMSD <- RMSD
  # final.intermediate.features <- candidate.intermediate.features

  # 5. Select n.close.to.targets shortest distances
  dist.order <- order(RMSD) # Ordering the squared distances from small to big.
  selected.distances <- dist.order[1:n.close.to.targets]
  sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]

  # 7. Put target features in dataframe format
  targets.empirical.df <- as.data.frame(matrix(targets.empirical, ncol = length(targets.empirical)))
  names(targets.empirical.df) <- y.names

  # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
  df.give.to.mice <- gtools::smartbind(dplyr::select(sim.results.with.design.df.selected,
                                                     -one_of(c("RMSD", "seed"))), # adding target to training dataset
                                       targets.empirical.df[rep(1:nrow(targets.empirical.df),
                                                                each = n.experiments), ])

  if (!identical(strict.positive.params, 0)){
    df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
  }
  if (!identical(probability.params, 0)){
    df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation
  }

  # 9. Override default predictorMatrix with a sparser matrix
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
      # The results of cv.glmnet are random. Perhaps we can solve this by setting a seed?
      set.seed(0) # for reproducibility
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
                                   printFlag = FALSE,
                                   seed = 0),
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
    experiments.df <- data.frame(experiments)

    # We could add an argument to the function to force the new experiments to respect the boundaries of the prior distributions.
    within.prior.limits <- rep(TRUE, n.experiments)
    if (inside_prior == TRUE){ # experiments.df is a dataframe with n.experiments rows and length(lls) columns
      params.above.lls <- sweep(x = experiments.df, MARGIN = 2, lls) %>% sign() %>% rowSums()
      params.below.uls <- sweep(x = -experiments.df, MARGIN = 2, -uls) %>% sign() %>% rowSums()
      within.prior.limits <- params.above.lls %in% length(lls) & params.below.uls %in% length(uls)
      experiments.df <- experiments.df[within.prior.limits, ]
    }

    set.seed(0) # for reproducibility
    experiments <- dplyr::sample_n(experiments.df,
                                   size = n.experiments,
                                   replace = TRUE,
                                   weight = rep(1, nrow(experiments.df))) %>%  #mice.test$imp.rnorm.values.weights[within.prior.limits]) %>%
      unlist %>%
      matrix(., byrow = FALSE, ncol = length(x.names))



    # From here is old code of HPC_MaC
    experiments <- as.data.frame(cbind(1:n.experiments, experiments))
    experiments.col.names <- c("seed", x.names)
    names(experiments) <- experiments.col.names
    write.table(experiments,
                sep = ",",
                dec = ".",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                file = file.name)
  }
  return(experiments)
}
