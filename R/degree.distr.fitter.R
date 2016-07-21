#' Fitting parametric distributions to the degree data
#'
#' The user fits various parametric distributions to the degree data and
#' compare their goodness of fit.To choose on the best fitting distribution,
#' a table is returned giving the empirical estimate of the three considered
#' goodness-of-fit statistics.
#'
#' @param degree.df The dataframe that is produced by \code{\link{degree.df.maker()}}
#'
#'
#' @return A table of the empirical estimate of the three considered goodness-of-fit statistics.
#'
#' @examples
#' load(dataframe)
#' bestfit <- degree.distr.fitter(dataframe = degree.df)
#'
#' @import fitdistrplus


# Fitting the parametric distributions to the Degree data by
# maximum likelihood estimation and comparison of goodness of fit statistics.

degree.distr.fitter <- function(dataframe = degree.df){
  # Normal distribution
  (fitn<-fitdist(dataframe$Degree, "norm"))
  gofstat(fitn)

  # Poisson distribution
  (fitpois<-fitdist(dataframe$Degree, "pois"))
  gofstat(fitpois)

  # Negative Bionomial distribution
  (fitnb<-fitdist(dataframe$Degree, "nbinom"))
  gofstat(fitnb)

  # Weibull distribution
  (fitw<-fitdist(dataframe$Degree, "weibull"))
  gofstat(fitw)

  # Gamma distribution
  (fitg<-fitdist(dataframe$Degree, "gamma"))
  gofstat(fitg)

  # Log-normal distribution
  (fitln<-fitdist(dataframe$Degree, "lnorm"))
  gofstat(fitln)

  # Beta distribution

  # (fitb<-fitdist(dataframe$degree, "beta"))
  # gofstat(fitb) # values must be in [0-1]

  (fit<-gofstat(list(fitn,fitpois,fitnb,fitw,fitg,fitln)))
  return(fit)
}

