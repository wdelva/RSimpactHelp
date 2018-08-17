#' adjusted mice sampler
#'
#' Augmented to also return details of the regression model and residual error sampled.
#'
#' @param data dataframe with missing data in it
#' @param m number of cycles I think
#' @param where points to missing data
#' @param imp number of imputations
#' @param blocks blocks
#' @param method method
#' @param visitSequence Default is 1 to ncol(data)
#' @param predictorMatrix predictorMatrix
#' @param formulas formulas
#' @param blots blots
#' @param post post
#' @param fromto fromto
#' @param printFlag Want a flag
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice
#' @export

sampler.fit <- function(data, m, where, imp, blocks, method, visitSequence,
          predictorMatrix, formulas, blots, post, fromto, printFlag,
          ...)
{
  initialize.chain <- mice:::initialize.chain
  handles.format <- mice:::handles.format
  is.passive <- mice:::is.passive
  sampler.univ <- mice:::sampler.univ
  check.df <- mice:::check.df
  remove.lindep <- mice:::remove.lindep

  imp.fit <- imp # making a copy of imp, to store parameters of fitted regression model
  imp.rnorm.values <- imp


  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  r <- !is.na(data)
  chainMean <- chainVar <- initialize.chain(blocks, maxit,
                                            m)
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
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            y <- data[, j]
            ry <- r[, j]
            wy <- where[, j]
            data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                            i]
          }
        }
        for (h in visitSequence) {
          ct <- attr(blocks, "calltype")
          calltype <- ifelse(length(ct) == 1, ct[1],
                             ct[h])
          b <- blocks[[h]]
          if (calltype == "formula")
            ff <- formulas[[h]]
          else ff <- NULL
          if (calltype == "type")
            type <- predictorMatrix[h, ]
          else type <- NULL
          user <- blots[[h]]
          theMethod <- method[h]
          empt <- theMethod == "" || (!is.null(type) &&
                                        all(type == 0))
          univ <- !empt && !is.passive(theMethod) &&
            !handles.format(paste0("mice.impute.", theMethod, ".fit"))
          mult <- !empt && !is.passive(theMethod) &&
            handles.format(paste0("mice.impute.", theMethod, ".fit"))
          pass <- !empt && is.passive(theMethod) && length(blocks[[h]]) ==
            1
          if (printFlag & !empt)
            cat(" ", b)
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(it = k, im = i, dep = h, meth = theMethod,
                           log = oldstate$log)
          assign("state", newstate, pos = parent.frame(),
                 inherits = TRUE)
          if (univ) {
            for (j in b) {
              imputation.list <- sampler.univ.fit(data = data,
                                              r = r, where = where, type = type, formula = ff,
                                              method = theMethod, yname = j, k = k,
                                              calltype = calltype, user = user, ...)
              imp[[j]][, i] <- imputation.list$imputes
              imp.fit[[j]] <- imputation.list$parm
              imp.rnorm.values[[j]][, i] <- imputation.list$rnorm.values
              data[(!r[, j]) & where[, j], j] <- imp[[j]][(!r[,
                                                              j])[where[, j]], i]
              cmd <- post[j]
              if (cmd != "") {
                eval(parse(text = cmd))
                data[where[, j], j] <- imp[[j]][, i]
              }
            }
          }
          if (mult) {
            mis <- !r
            mis[, setdiff(b, colnames(data))] <- FALSE
            data[mis] <- NA
            fm <- paste("mice.impute", theMethod, "fit", sep = ".")
            if (calltype == "formula")
              imputes <- do.call(fm, args = list(data = data,
                                                 formula = ff, ...))

            else if (calltype == "type"){
              imputation.list <- do.call(fm, args = list(data = data,
                                                         type = type, ...))
              imputes <- imputation.list$imputation
            }
            else stop("Cannot call function of type ",
                      calltype, call. = FALSE)
            if (is.null(imputes))
              stop("No imputations from ", theMethod,
                   h, call. = FALSE)
            for (j in names(imputes)) {
              imp[[j]][, i] <- imputes[[j]]
              imp.fit[[j]] <- imputation.list$parm
              imp.rnorm.values[[j]][, i] <- imputation.list$rnorm.values
              data[!r[, j], j] <- imp[[j]][, i]
            }
          }
          if (pass) {
            for (j in b) {
              wy <- where[, j]
              ry <- r[, j]
              imp[[j]][, i] <- model.frame(as.formula(theMethod),
                                           data[wy, ], na.action = na.pass)
              data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                              i]
            }
          }
        }
      }
      k2 <- k - from + 1L
      if (length(visitSequence) > 0L) {
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            if (!is.factor(data[, j])) {
              chainVar[j, k2, ] <- apply(imp[[j]], 2L,
                                         var, na.rm = TRUE)
              chainMean[j, k2, ] <- colMeans(as.matrix(imp[[j]]),
                                             na.rm = TRUE)
            }
            if (is.factor(data[, j])) {
              for (mm in seq_len(m)) {
                nc <- as.integer(factor(imp[[j]][, mm],
                                        levels = levels(data[, j])))
                chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
                chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
              }
            }
          }
        }
      }
    }
    if (printFlag) {
      r <- get("loggedEvents", parent.frame(1))
      ridge.used <- any(grepl("A ridge penalty", r$out))
      if (ridge.used) {
        cat("\n * Please inspect the loggedEvents \n")
      }
      else {
        cat("\n")
      }
    }
  }
  return(list(iteration = maxit, imp = imp,
              imp.fit = imp.fit,
              imp.rnorm.values = imp.rnorm.values,
              chainMean = chainMean,
              chainVar = chainVar))
}
