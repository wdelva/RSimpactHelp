#' adjusted mice sampler
#'
#' Augmented to also return details of the regression model and residual error sampled.
#'
#' @param p padded object
#' @param data dataframe with missing data in it
#' @param where points to missing data
#' @param m number of cycles I think
#' @param imp number of imputations
#' @param r Forgot
#' @param visitSequence Default is 1 to ncol(data)
#' @param fromto Forgot
#' @param printFlag Want a flag
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice
#' @export

sampler.fit <- function(p, data, where, m, imp, r, visitSequence, fromto, printFlag,
          ...)
{

  is.passive <- mice:::is.passive
  check.df <- mice:::check.df
  remove.lindep <- mice:::remove.lindep


  imp.fit <- imp # making a copy of imp, to store parameters of fitted regression model
  imp.rnorm.values <- imp
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  chainVar <- chainMean <- NULL
  if (maxit > 0)
    chainVar <- chainMean <- array(0, dim = c(length(visitSequence),
                                              maxit, m), dimnames = list(dimnames(data)[[2]][visitSequence],
                                                                         seq_len(maxit), paste("Chain", seq_len(m))))
  if (maxit < 1)
    iteration <- 0
  if (maxit >= 1) {
    if (printFlag)
      cat("\n iter imp variable")
    for (k in from:to) {
      iteration <- k
      for (i in seq_len(m)) {
        if (printFlag)
          cat("\n ", iteration, " ", i)
        for (j in visitSequence) {
          wy <- where[, j]
          ry <- r[, j]
          p$data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                            i]
        }
        for (j in setdiff(p$visitSequence, visitSequence)) {
          cat.columns <- p$data[, p$categories[j, 4]]
          mm <- model.matrix(~cat.columns - 1, model.frame(~cat.columns,
                                                           na.action = na.pass))[, -1]
          p$data[, (j:(j + p$categories[p$categories[j,
                                                     4], 2] - 1))] <- mm
        }
        for (j in p$visitSequence) {
          theMethod <- p$method[j]
          vname <- dimnames(p$data)[[2]][j]
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(it = k, im = i, co = j, dep = vname,
                           meth = theMethod, log = oldstate$log)
          assign("state", newstate, pos = parent.frame(),
                 inherits = TRUE)
          if (printFlag && theMethod != "dummy")
            cat(" ", vname)
          empt <- theMethod == ""
          elem <- !empt && !is.passive(theMethod) &&
            theMethod != "dummy"
          flat <- elem && substring(theMethod, 1, 2) !=
            "2l"
          pass <- !empt && is.passive(theMethod)
          dumm <- theMethod == "dummy"
          if (elem) {
            if (flat) {
              predictors <- p$predictorMatrix[j, ] ==
                1
            }
            else {
              predictors <- p$predictorMatrix[j, ] !=
                0
            }
            if (!is.null(p$form) && nchar(p$form[j]) >
                0) {
              myform <- paste(p$form[j], "0", sep = "+")
              x <- model.matrix(formula(myform), p$data)
            }
            else {
              x <- p$data[, predictors, drop = FALSE]
            }
            y <- p$data[, j]
            ry <- complete.cases(x, y) & r[, j]
            wy <- complete.cases(x) & where[, j]
            cc <- wy[where[, j]]
            type <- p$predictorMatrix[j, predictors]
            if (k == 1)
              check.df(x, y, ry)
            keep <- remove.lindep(x, y, ry, ...)
            x <- x[, keep, drop = FALSE]
            type <- type[keep]
            f <- paste("mice.impute", theMethod, "fit", sep = ".")
            imputes <- p$data[wy, j]
            imputes[!cc] <- NA
            imputation.list <- do.call(f, args = list(y,
                                                  ry, x, wy = wy, type = type, ...))
            imputes[cc] <- imputation.list$imputation
            imp[[j]][, i] <- imputes
            imp.fit[[j]] <- imputation.list$parm
            imp.rnorm.values[[j]][, i] <- imputation.list$rnorm.values
            p$data[(!r[, j]) & where[, j], j] <- imp[[j]][(!r[,
                                                              j])[where[, j]], i]

          }
          if (pass) {
            wy <- where[, j]
            imp[[j]][, i] <- model.frame(as.formula(theMethod),
                                         p$data[wy, ], na.action = na.pass)
            p$data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                              i]
          }
          if (dumm) {
            cat.columns <- p$data[, p$categories[j, 4]]
            mm <- model.matrix(~cat.columns - 1, model.frame(~cat.columns,
                                                             na.action = na.pass))[, -1]
            p$data[, (j:(j + p$categories[p$categories[j,
                                                       4], 2] - 1))] <- mm
            remove("cat.columns")
          }
          cmd <- p$post[j]
          if (cmd != "") {
            eval(parse(text = cmd))
            p$data[where[, j], j] <- imp[[j]][, i]
          }
        }
      }
      k2 <- k - from + 1
      if (length(visitSequence) > 0) {
        for (j in seq_along(visitSequence)) {
          jj <- visitSequence[j]
          if (!is.factor(data[, jj])) {
            chainVar[j, k2, ] <- apply(imp[[jj]], 2,
                                       var, na.rm = TRUE)
            chainMean[j, k2, ] <- colMeans(as.matrix(imp[[jj]]),
                                           na.rm = TRUE)
          }
          if (is.factor(data[, jj])) {
            for (mm in seq_len(m)) {
              nc <- as.integer(factor(imp[[jj]][, mm],
                                      levels = levels(data[, jj])))
              chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
              chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
            }
          }
        }
      }
    }
    if (printFlag)
      cat("\n")
  }
  return(list(iteration = maxit, imp = imp, imp.fit = imp.fit, imp.rnorm.values = imp.rnorm.values, chainMean = chainMean,
              chainVar = chainVar))
}
