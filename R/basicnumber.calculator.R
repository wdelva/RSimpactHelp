#' Basci reproductive number
#'
#' Compute the  basic reproductive number of a dynamic sexual network simulated with Simpact
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param beta The probability that an infected individual will infect a suceptible partner over the duration of
#' their relationship
#' @param endpoint Only transmission events that took place before this point in simulation time,
#' are captured in the output
#' @return A value of the basic reproduction number
#'
#' @example ro <- basicnumber.calculator(datalist = datalist, beta = 0.1508)
#'
#'

basicnumber.calculator <- function(datalist = datalist, beta = 0.1508){

  # Ro = beta * C * D

  # C: effective mean of the average of annual number of new partners C = m + (sigma^2)/m
  # where m is the arithmetic mean of the annual number of new partners

  # D: average duration of infectiousness
  # Note: C < 1/(beta * D)

  source("/home/david/RSimpactHelp/R/transmission.rate.calculator.R")

  trans.rate.int <- transmission.rate.calculator(datalist = datalist, timewindow = c(0,40), int = TRUE, by = 1)

  m <- mean(trans.rate.int)

  sigma <- sd(trans.rate.int)

  C <- (m + (sigma^2)/m)



}






