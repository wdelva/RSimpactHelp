#' Create a default parameter configuration list.
#'
#' Create a default parameter configuration list as a starting point for a parameter optimisation algorithm.
#' At this stage, only the most relevant parameters are modifiable.
#'
#' @param conception.alpha_base Baseline parameter for conception rate. (-3)
#' @param conception.alpha_agewoman The Effect of the age of the woman on the conception rate. (0)
#' @param dissolution.alpha_0 Baseline parameter for relationship dissolution rate. (0.1)
#' @param dissolution.alpha_4 Effect of increasing mean age of the  couple on  relationship dissolution rate (-0.05)
#' @param formation.hazard.type Type of hazard function for relationship formation. Choose between "simple", "agegap" and "agegapry".
#' @param formation.hazard.agegapry.gap_factor_man_const Baseline parameter for men (0).
#' @param formation.hazard.agegapry.gap_factor_woman_const Baseline parameter for women (0).
#' @param person.agegap.man.dist.type Distribution of preferred age differences for men ("normal")
#' @param person.agegap.woman.dist.type Distribution of preferred age differences for women ("normal")
#' @param person.agegap.man.dist.normal.mu Mean of preferred age differences distribution for men (-4)
#' @param person.agegap.woman.dist.normal.mu Mean of preferred age differences distribution for women (-4)
#' @param person.agegap.man.dist.normal.sigma Standard deviation of preferred age differences distribution for men (3)
#' @param person.agegap.woman.dist.normal.sigma Standard deviation of preferred age differences distribution for women (3)
#' @param formation.hazard.agegapry.baseline Baseline parameter for relationship formation rate (0)
#' @param formation.hazard.agegapry.numrel_scale_man Effect of the number of relationships on the penalty for deviation from preferred age difference for men (0)
#' @param formation.hazard.agegapry.numrel_scale_woman Effect of the number of relationships on the penalty for deviation from preferred age difference for women (0)
#' @param formation.hazard.agegapry.gap_agescale_man Effect of male age on preferred age difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
#' @param formation.hazard.agegapry.gap_agescale_woman Effect of male age on preferred age difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
#' @param person.eagerness.dist.type Type of distribution for "eagerness" (~ sexual activity) in the population ("gamma")
#' @param formation.hazard.agegapry.eagerness_sum Effect of the sum of eagerness values in both partners on the relationship formation rate (1)
#' @param transmission.param.a Baseline parameter for HIV transmission rate in serodiscordant couples (-1.0352239)
#' @param transmission.param.b Parameter "b" for the linear component of the effect of viral load on the HIV transmission rate in serodiscordant couples (-89.339994)
#' @param transmission.param.c Parameter "c" for the exponential component of the effect of viral load on the HIV transmission rate in serodiscordant couples (-0.4948478)
#' @param transmission.param.d1 Parameter "d1" for effect of number of ongoing relationships the the male partner (coital dilution) (0)
#' @param transmission.param.d2 Parameter "d2" for effect of number of ongoing relationships the the female partner (coital dilution) (0)
#' @param transmission.param.f1 Effect of youngest age on HIV susceptibility (log(5) ~1.6 such that the hazard is x 5 in 15 year olds)
#' @param transmission.param.f2 Effect of female age on HIV susceptibility (log(log(2.5) / log(5)) / 5 ~-0.11 such that the hazard is x 2.5 in 20 year olds, compared to the reference (>>25 year olds)
#' @param formation.hazard.agegapry.numrel_man Effect of number of ongoing relationships on the relationship formation rate for men (-0.2)
#' @param formation.hazard.agegapry.numrel_woman Effect of number of ongoing relationships on the relationship formation rate for women (-0.2)
#' @param formation.hazard.agegapry.numrel_diff Effect of absolute difference in number of ongoing relationships on the relationship formation rate (-0.1)
#' @param person.eagerness.dist.gamma.a Shape parameter (kappa parameter on wiki page) of Gamma distribution of "eagerness" (0.125)
#' @param person.eagerness.dist.gamma.b Scale parameter (theta parameter on wiki page) of Gamma distribution of "eagerness" (8)
#' @param formation.hazard.agegapry.eagerness_diff Effect of absolute difference in eagerness on the relationship formation rate (-1)
#' @param formation.hazard.agegapry.gap_factor_man_exp Effect of deviation from preferred age difference on the relationship formation rate for men (-0.2)
#' @param formation.hazard.agegapry.gap_factor_woman_exp Effect of deviation from preferred age difference on the relationship formation rate for women (-0.2)
#' @param formation.hazard.agegapry.gap_factor_man_age Funnel effect in men (0.05)
#' @param formation.hazard.agegapry.gap_factor_woman_age Funnel effect in women (0.05)
#' @param formation.hazard.agegapry.meanage Effect of increasing mean age of the candidate couple on the relationship formation rate (-0.1)

#' @return a list of model parameters that can be used as input for simpact.run()
#' @examples
#' cfg.list <- input.params.creator(conception.alpha_base = -3,
#' formation.hazard.type = "agegapry",
#' formation.hazard.agegapry.numrel_man = -1,
#' formation.hazard.agegapry.numrel_woman = -1, ...)

input.params.creator <- function(population.eyecap.fraction = 0.5,
                                 population.simtime = 40,
                                 population.numwomen = 200,
                                 population.nummen = 200,
                                 periodiclogging.interval = 1,
                                 periodiclogging.starttime = 0,
                                 syncrefyear.interval = 1,
                                 syncpopstats.interval = 1,
                                 hivseed.time = 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 10,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 30,
                                 conception.alpha_base = -3,
                                 conception.alpha_agewoman = 0,
                                 dissolution.alpha_0 = 0.1,
                                 dissolution.alpha_4 = -0.05,
                                 formation.hazard.type = "agegapry",
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 person.agegap.man.dist.type = "normal",
                                 person.agegap.woman.dist.type = "normal",
                                 person.agegap.man.dist.normal.mu = -4,
                                 person.agegap.woman.dist.normal.mu = -4,
                                 person.agegap.man.dist.normal.sigma = 3,
                                 person.agegap.woman.dist.normal.sigma = 3,
                                 formation.hazard.agegapry.baseline = 0,
                                 formation.hazard.agegapry.numrel_scale_man = 0,
                                 formation.hazard.agegapry.numrel_scale_woman = 0,
                                 formation.hazard.agegapry.gap_agescale_man = 0.3,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.3,
                                 person.eagerness.dist.type = "gamma",
                                 formation.hazard.agegapry.eagerness_sum = 1,
                                 transmission.param.a = -1.0352239,
                                 transmission.param.b = -89.339994,
                                 transmission.param.c = 0.4948478,
                                 transmission.param.d1 = 0,
                                 transmission.param.d2 = 0,
                                 transmission.param.f1 = log(5),
                                 transmission.param.f2 = log(log(2.5) / log(5)) / 5,
                                 formation.hazard.agegapry.numrel_man = -0.2,
                                 formation.hazard.agegapry.numrel_woman = -0.2,
                                 formation.hazard.agegapry.numrel_diff = -0.1,
                                 person.eagerness.dist.gamma.a = 0.125,
                                 person.eagerness.dist.gamma.b = 8,
                                 formation.hazard.agegapry.eagerness_diff = -1,
                                 formation.hazard.agegapry.gap_factor_man_exp = -0.2,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -0.2,
                                 formation.hazard.agegapry.gap_factor_man_age = 0.05,
                                 formation.hazard.agegapry.gap_factor_woman_age = 0.05,
                                 formation.hazard.agegapry.meanage = -0.1
                                 ){
  input.params.list <- list()
  input.params.list$population.eyecap.fraction <- population.eyecap.fraction
  input.params.list$population.simtime <- population.simtime
  input.params.list$population.numwomen <- population.numwomen
  input.params.list$population.nummen <- population.nummen
  input.params.list$periodiclogging.interval <- periodiclogging.interval
  input.params.list$periodiclogging.starttime <- periodiclogging.starttime
  input.params.list$syncrefyear.interval <- syncrefyear.interval
  input.params.list$syncpopstats.interval <- syncpopstats.interval
  input.params.list$hivseed.time <- hivseed.time
  input.params.list$hivseed.type <- hivseed.type
  input.params.list$hivseed.amount <- hivseed.amount
  input.params.list$hivseed.age.min <- hivseed.age.min
  input.params.list$hivseed.age.max <- hivseed.age.max
  input.params.list$conception.alpha_base <- conception.alpha_base
  input.params.list$conception.alpha_agewoman <- conception.alpha_agewoman
  input.params.list$dissolution.alpha_0 <- dissolution.alpha_0
  input.params.list$dissolution.alpha_4 <- dissolution.alpha_4
  input.params.list$formation.hazard.type <- formation.hazard.type
  input.params.list$formation.hazard.agegapry.gap_factor_man_const <- formation.hazard.agegapry.gap_factor_man_const
  input.params.list$formation.hazard.agegapry.gap_factor_woman_const <- formation.hazard.agegapry.gap_factor_woman_const
  input.params.list$person.agegap.man.dist.type <- person.agegap.man.dist.type
  input.params.list$person.agegap.woman.dist.type <- person.agegap.woman.dist.type
  input.params.list$person.agegap.man.dist.normal.mu <- person.agegap.man.dist.normal.mu
  input.params.list$person.agegap.woman.dist.normal.mu <- person.agegap.woman.dist.normal.mu
  input.params.list$person.agegap.man.dist.normal.sigma <- person.agegap.man.dist.normal.sigma
  input.params.list$person.agegap.woman.dist.normal.sigma <- person.agegap.woman.dist.normal.sigma
  input.params.list$formation.hazard.agegapry.baseline <- formation.hazard.agegapry.baseline
  input.params.list$formation.hazard.agegapry.numrel_scale_man <- formation.hazard.agegapry.numrel_scale_man
  input.params.list$formation.hazard.agegapry.numrel_scale_woman <- formation.hazard.agegapry.numrel_scale_woman
  input.params.list$formation.hazard.agegapry.gap_agescale_man <- formation.hazard.agegapry.gap_agescale_man
  input.params.list$formation.hazard.agegapry.gap_agescale_woman <- formation.hazard.agegapry.gap_agescale_woman
  input.params.list$person.eagerness.dist.type <- person.eagerness.dist.type
  input.params.list$formation.hazard.agegapry.eagerness_sum <- formation.hazard.agegapry.eagerness_sum
  input.params.list$transmission.param.a <- transmission.param.a
  input.params.list$transmission.param.b <- transmission.param.b
  input.params.list$transmission.param.c <- transmission.param.c
  input.params.list$transmission.param.d1 <- transmission.param.d1
  input.params.list$transmission.param.d2 <- transmission.param.d2
  input.params.list$transmission.param.f1 <- transmission.param.f1
  input.params.list$transmission.param.f2 <- transmission.param.f2
  input.params.list$formation.hazard.agegapry.numrel_man <- formation.hazard.agegapry.numrel_man
  input.params.list$formation.hazard.agegapry.numrel_woman <- formation.hazard.agegapry.numrel_woman
  input.params.list$formation.hazard.agegapry.numrel_diff <- formation.hazard.agegapry.numrel_diff
  input.params.list$person.eagerness.dist.gamma.a <- person.eagerness.dist.gamma.a
  input.params.list$person.eagerness.dist.gamma.b <- person.eagerness.dist.gamma.b
  input.params.list$formation.hazard.agegapry.eagerness_diff <- formation.hazard.agegapry.eagerness_diff
  input.params.list$formation.hazard.agegapry.gap_factor_man_exp <- formation.hazard.agegapry.gap_factor_man_exp
  input.params.list$formation.hazard.agegapry.gap_factor_woman_exp <- formation.hazard.agegapry.gap_factor_woman_exp
  input.params.list$formation.hazard.agegapry.gap_factor_man_age <- formation.hazard.agegapry.gap_factor_man_age
  input.params.list$formation.hazard.agegapry.gap_factor_woman_age <- formation.hazard.agegapry.gap_factor_woman_age
  input.params.list$formation.hazard.agegapry.meanage <- formation.hazard.agegapry.meanage
  return(input.params.list)
}
