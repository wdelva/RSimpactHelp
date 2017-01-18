#!/usr/bin/env/ Rscript
#' These are the preset boundaries for most of the parameters that are likely to be calibrated in simpact.
#' Most of these parameters will be standard in \code{\link{simpact.config.inputs()}} stage.
#' The user can set the simulation specific parameters as needed
#'
#' @param default is to return all the parameters used during the \code{\link{simpact.config.inputs()}} setup.
#' @return a list of boundaries of parameters to be calibrated
#' @examples
#' varied.params.bound <- simpact.params.boundaries()[[1]]

simpact.params.boundaries <- function(...){

  arg.list <- list(person.eagerness.man.dist.gamma.a =  c(0.6, 1.9), #c(0.1, 2),
                   person.eagerness.man.dist.gamma.b =  c(30, 60), #c(5, 60),
                   person.eagerness.woman.dist.gamma.a =  c(0.6, 1.9), #c(0.5, 2),
                   person.eagerness.woman.dist.gamma.b =  c(30, 60), #c(5, 60),
                   conception.alpha_base =  c(-3.2, -1.2), #c(-3.6, -1.2),
                   formation.hazard.agegapry.numrel_man = c(-1.5, -0.1),
                   formation.hazard.agegapry.numrel_woman = c(-1.5, -0.1),
                   formation.hazard.agegapry.eagerness_diff =  c(-0.06, 0), #c(-0.1, 0),
                   formation.hazard.agegapry.gap_factor_man_exp = c(-1.5, -0.4),
                   formation.hazard.agegapry.gap_factor_woman_exp = c(-1.5, -0.4),
                   person.agegap.man.dist.normal.mu =   c(2, 6), #c(0, 6),
                   person.agegap.woman.dist.normal.mu =  c(2 ,6), #c(0, 6),
                   person.agegap.man.dist.normal.sigma =  c(1, 2), #c(0.5, 2),
                   person.agegap.woman.dist.normal.sigma = c(1, 2), #c(0.5, 2)
                   person.art.accept.threshold.dist.fixed.value = c(0.4, 0.9),
                   ...)

  return(arg.list)
}

