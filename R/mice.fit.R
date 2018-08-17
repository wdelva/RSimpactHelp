#' MICE with weights
#'
#' Produces a list with imputations with inverse probability weights.
#'
#' @param data training dataframe
#' @param m number of imputations
#' @param method string (vector)
#' @param predictorMatrix predictorMatrix
#' @param blocks blocks
#' @param visitSequence defaults is 1 to ncol(data)
#' @param formulas formulas
#' @param blots blots
#' @param post post
#' @param defaultMethod c("pmm", "logreg", "polyreg", "polr")
#' @param maxit How many cycles through chained equations
#' @param printFlag Do you want stuff printed
#' @param seed seed value to make it reproducible
#' @param data.init Same
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice
#' @export

mice.fit <- function(data, m = 5, method = "norm", predictorMatrix, where = NULL,
          blocks, visitSequence = NULL, formulas, blots = NULL, post = NULL,
          defaultMethod = "norm", maxit = 5,
          printFlag = TRUE, seed = NA, data.init = NULL, ...)
{
  mice.impute.norm <- mice::mice.impute.norm
  check.dataform <- mice:::check.dataform
  check.m <- mice:::check.m
  make.blocks <- mice:::make.blocks
  check.blocks <- mice:::check.blocks
  make.predictorMatrix <- mice:::make.predictorMatrix
  check.predictorMatrix <- mice:::check.predictorMatrix
  check.formulas <- mice:::check.formulas
  make.formulas <- mice:::make.formulas
  construct.blocks <- mice:::construct.blocks
  check.cluster <- mice:::check.cluster
  check.where <- mice:::check.where
  check.visitSequence <- mice:::check.visitSequence
  check.method <- mice:::check.method
  check.post <- mice:::check.post
  check.blots <- mice:::check.blots
  edit.setup <- mice:::edit.setup
  initialize.imp <- mice:::initialize.imp
  check.data <- mice:::check.data
  is.passive <- mice:::is.passive
  check.df <- mice:::check.df
  remove.lindep <- mice:::remove.lindep

  call <- match.call()
  if (!is.na(seed))
    set.seed(seed)
  data <- check.dataform(data)
  m <- check.m(m)
  mp <- missing(predictorMatrix)
  mb <- missing(blocks)
  mf <- missing(formulas)
  if (mp & mb & mf) {
    blocks <- make.blocks(colnames(data))
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  if (!mp & mb & mf) {
    predictorMatrix <- check.predictorMatrix(predictorMatrix,
                                             data)
    blocks <- make.blocks(colnames(predictorMatrix), partition = "scatter")
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (mp & !mb & mf) {
    blocks <- check.blocks(blocks, data)
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  if (mp & mb & !mf) {
    formulas <- check.formulas(formulas, data)
    blocks <- construct.blocks(formulas)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & mf) {
    blocks <- check.blocks(blocks, data)
    z <- check.predictorMatrix(predictorMatrix, data, blocks)
    predictorMatrix <- z$predictorMatrix
    blocks <- z$blocks
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (!mp & mb & !mf) {
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix,
                                             data)
    blocks <- construct.blocks(formulas, predictorMatrix)
  }
  if (mp & !mb & !mf) {
    blocks <- check.blocks(blocks, data, calltype = "formula")
    formulas <- check.formulas(formulas, blocks)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & !mf) {
    blocks <- check.blocks(blocks, data)
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix,
                                             data, blocks)
  }
  chk <- check.cluster(data, predictorMatrix)
  where <- check.where(where, data, blocks)
  visitSequence <- check.visitSequence(visitSequence, data = data,
                                       where = where, blocks = blocks)
  method <- check.method(method = method, data = data, where = where,
                         blocks = blocks, defaultMethod = defaultMethod)
  post <- check.post(post, data)
  blots <- check.blots(blots, data, blocks)
  state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "",
                             out = "")
  setup <- list(method = method, predictorMatrix = predictorMatrix,
                visitSequence = visitSequence, post = post)
  setup <- edit.setup(data, setup, ...)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  nmis <- apply(is.na(data), 2, sum)
  imp <- initialize.imp(data, m, where, blocks, visitSequence,
                        method, nmis, data.init)
  from <- 1
  to <- from + maxit - 1
  # replacing sampler with sampler.fit
  q <- sampler.fit(data, m, where, imp, blocks, method, visitSequence,
               predictorMatrix, formulas, blots, post, c(from, to),
               printFlag, ...)
  imp.fit <- q$imp.fit # New object with regression models
  imp.rnorm.values <- q$imp.rnorm.values # New object with rnorm.values (for which the density can be calculated)
  # Let's turn this list of rnorm.values with one single-column dataframe per missing variable into a single vector or
  # dataframe with one value per imputation. The values are simply the product of the r.norm.values, divided by
  # the sum of all r.norm.values
  imp.rnorm.values.matrix <- matrix(dnorm(as.numeric(unlist(imp.rnorm.values)), log = TRUE),
                                    nrow = max(apply(is.na(data), 2, sum)))
  # Instead of taking the product of the densities, we take the sum of the logdensities
  imp.rnorm.values.sumlog <- apply(imp.rnorm.values.matrix, 1, FUN = sum)
  imp.rnorm.values.logweights <- imp.rnorm.values.sumlog - min(imp.rnorm.values.sumlog)
  imp.rnorm.values.weights <- exp(-imp.rnorm.values.logweights)  # The weights are inverse to their density
  if (!state$log)
    loggedEvents <- NULL
  if (state$log)
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  midsobj <- list(data = data, imp = q$imp,
                  imp.fit = imp.fit,
                  imp.rnorm.values = imp.rnorm.values,
                  imp.rnorm.values.weights = imp.rnorm.values.weights,
                  m = m, where = where,
                  blocks = blocks, call = call, nmis = nmis, method = method,
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                  formulas = formulas, post = post, blots = blots, seed = seed,
                  iteration = q$iteration, lastSeedValue = .Random.seed,
                  chainMean = q$chainMean, chainVar = q$chainVar, loggedEvents = loggedEvents,
                  version = packageVersion("mice"), date = Sys.Date())
  oldClass(midsobj) <- "mids"
  if (!is.null(midsobj$loggedEvents))
    warning("Number of logged events: ", nrow(midsobj$loggedEvents),
            call. = FALSE)
  return(midsobj)
}
