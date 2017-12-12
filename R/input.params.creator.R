#' Create a default parameter configuration list.
#'
#' Create a default parameter configuration list as a starting point for a parameter
#' optimisation algorithm.
#' At this stage, only the most relevant parameters are modifiable.
#'
#' @param conception.alpha_base Baseline parameter for conception rate. (-3)
#' @param dissolution.alpha_0 Baseline parameter for relationship dissolution rate. (0.1)
#' @param dissolution.alpha_4 Effect of increasing mean age of the  couple on  relationship
#' dissolution rate (-0.05)
#' @param formation.hazard.type Type of hazard function for relationship formation.
#' Choose between "simple", "agegap" and "agegapry".
#' @param person.art.accept.threshold.dist.type This specfies the ART acceptance distribution
#' for each person
#' @param person.art.accept.threshold.dist.fixed.value The value for the fixed distribution
#' @param formation.hazard.agegapry.gap_factor_man_const Baseline parameter for men (0).
#' @param formation.hazard.agegapry.gap_factor_woman_const Baseline parameter for women (0).
#' @param person.agegap.man.dist.type Distribution of preferred age differences
#' for men ("normal")
#' @param person.agegap.woman.dist.type Distribution of preferred age differences
#' for women ("normal")
#' @param person.agegap.man.dist.normal.mu Mean of preferred age differences distribution
#' for men (-4)
#' @param person.agegap.woman.dist.normal.mu Mean of preferred age differences distribution
#' for women (-4)
#' @param birth.pregnancyduration.dist.type duration of pregnancy
#' @param birth.pregnancyduration.dist.fixed.value the option is fixed and set to 9 months
#' @param person.agegap.man.dist.normal.sigma Standard deviation of preferred age differences
#'  distribution for men (3)
#' @param person.agegap.woman.dist.normal.sigma Standard deviation of preferred age
#' differences distribution for women (3)
#' @param formation.hazard.agegapry.gap_agescale_man Effect of male age on preferred age
#' difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
#' @param formation.hazard.agegapry.gap_agescale_woman Effect of male age on preferred
#' age difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
#' @param person.eagerness.man.type Type of man distribution for "eagerness" (~ sexual activity)
#' in the population ("independent")
#' @param person.eagerness.woman.type Type of man distribution for "eagerness" (~ sexual activity)
#'  in the population ("independent")
#' @param person.eagerness.man.dist.type Type of man distribution for "eagerness"
#' (~ sexual activity) in the population ("gamma")
#' @param person.eagerness.woman.dist.type Type of women distribution for "eagerness"
#' (~ sexual activity) in the population ("gamma")
#' @param formation.hazard.agegapry.eagerness_sum Effect of the sum of eagerness values
#' in both partners on the relationship formation rate (1)
#' @param hivtransmission.param.a Baseline parameter for HIV transmission rate in serodiscordant
#'  couples (-1.0352239)
#' @param hivtransmission.param.b Parameter "b" for the linear component of the effect of viral
#' load on the HIV transmission rate in serodiscordant couples (-89.339994)
#' @param hivtransmission.param.c Parameter "c" for the exponential component of the effect of
#' viral load on the HIV transmission rate in serodiscordant couples (0.4948478)
#' @param hivtransmission.param.f1 Effect of youngest age on HIV susceptibility (log(5) ~1.6
#' such that the hazard is x 5 in 15 year olds)
#' @param hivtransmission.param.f2 Effect of female age on HIV susceptibility
#' (log(log(2.5) / log(5)) / 5 ~-0.11 such that the hazard is x 2.5 in 20 year olds,
#' compared to the reference (>>25 year olds)
#' @param person.vsp.toaids.x Effect of "initial" AIDS stage versus chronic HIV infection
#'  on infectiousness (7)
#' @param person.vsp.tofinalaids.x Effect of "final" AIDS stage versus chronic HIV
#' infection on infectiousness (12)
#' @param formation.hazard.agegapry.numrel_man Effect of number of ongoing relationships on
#' the relationship formation rate for men (-0.2)
#' @param formation.hazard.agegapry.numrel_woman Effect of number of ongoing relationships
#' on the relationship formation rate for women (-0.2)
#' @param formation.hazard.agegapry.numrel_diff Effect of absolute difference in number of
#' ongoing relationships on the relationship formation rate (-0.1)
#' @param person.eagerness.man.dist.gamma.a Shape parameter (kappa parameter on wiki page) of
#'  Gamma distribution of "eagerness" (0.125)
#' @param person.eagerness.man.dist.gamma.b Scale parameter (theta parameter on wiki page) of
#' Gamma distribution of "eagerness" (8)
#' @param person.eagerness.woman.dist.gamma.a Shape parameter (kappa parameter on wiki page) of
#' Gamma distribution of "eagerness" (0.125)
#' @param person.eagerness.woman.dist.gamma.b Scale parameter (theta parameter on wiki page) of
#'  Gamma distribution of "eagerness" (8)
#' @param formation.hazard.agegapry.eagerness_diff Effect of absolute difference in eagerness
#'  on the relationship formation rate (-1)
#' @param formation.hazard.agegapry.gap_factor_man_exp Effect of deviation from preferred age
#' difference on the relationship formation rate for men (-0.2)
#' @param formation.hazard.agegapry.gap_factor_woman_exp Effect of deviation from preferred
#' age difference on the relationship formation rate for women (-0.2)
#' @param formation.hazard.agegapry.gap_factor_man_age Funnel effect in men (0.05)
#' @param formation.hazard.agegapry.gap_factor_woman_age Funnel effect in women (0.05)
#' @param formation.hazard.agegapry.meanage Effect of increasing mean age of the candidate
#' couple on the relationship formation rate (-0.1)
#' @param mortality.normal.weibull.shape Specifies the shape of a Weibull distribution to make
#'  sure everyone has a limited lifespan (5)
#' @param mortality.normal.weibull.scale Specifies the scale parameter to base the Weibull
#' distribution on (65)
#' @param mortality.normal.weibull.genderdiff Specifies if the gender needs to be taken into
#'  account for the scale parameter of the Weibull distribution (0)
#' @param periodiclogging.interval Allow the periodic log of the events and write teh output
#' to a csv file (1)
#' @param syncrefyear.interval Specifies the time interval with which a reference simulation
#'  time should be saved (1)
#' @param debut.debutage The age a person must have to become sexually active.(14)
#' @param population.simtime Indication of how many year should the simulation run for.(40)
#'  e.g 40 years
#' @param population.nummen Indication of how many men should be in the initial seeded
#' population (500) e.g 500 men
#' @param population.numwomen Indication of how many women should be in the initial seeded
#'  population (500) e.g 500 women
#' @param population.eyecap.fraction Allow for the indication of how many people can a person
#'  possible engage in a relation with. (0.2)
#' @param hivseed.type Allows for a selection of method on seeding HIV. (amount) e.g if amount
#'  number of people else if fraction then a fraction will be set to HIV+
#' @param hivseed.amount These number of people will be seeded to be HIV+ (20)
#' @param hivseed.age.min This allows for the minimal age to be set seeded as HIV+ (20)
#' @param hivseed.age.max Maximum age at which HIV seeding will be infected (30)
#' @param hivseed.time The time within the simulation when the HIV seeding is to take
#'  place (10) e.g after 10 years in the simulation
#' @param diagnosis.baseline Baseline parameter for the diagnosis event (-100)
#' @param monitoring.cd4.threshold The value of cd4 that triggers eligibility for treatment
#' if below this value
#' @param simulation.type To distinguish parameters that vary from the type of simulation
#' being perfomed (simpact-cyan) e.g you can use maxart as well
#' @param mortality.aids.survtime.C HIV-based time of death parameter
#' @param mortality.aids.survtime.k HIV-based time of death parameter
#' @param monitoring.fraction.log_viralload the setpoint viral load will be dropped by this fraction
#' @param birth.boygirlratio the probability of a newly born to be boy
#' @param monitoring.cd4.threshold.prestudy cd4 values prestudy
#' @param monitoring.cd4.threshold.instudy.controlstage cd4 values during the control phase
#'# @param monitoring.cd4.threshold.instudy.transitionstage cd4 threshold when site is transitioning
#'# @param monitoring.cd4.threshold.instudy.interventionstage cd threshold not needed during the intervention stage
#'
#'
#' @return a list of model parameters that can be used as input for simpact.run()
#' @examples
#' cfg.list <- input.params.creator(conception.alpha_base = -3,
#' formation.hazard.type = "agegapry",
#' formation.hazard.agegapry.numrel_man = -1,
#' formation.hazard.agegapry.numrel_woman = -1)
#' @export

input.params.creator <- function(mortality.normal.weibull.shape = 5,
                                 mortality.normal.weibull.scale = 65,
                                 mortality.normal.weibull.genderdiff = 0,
                                 periodiclogging.interval = 1,
                                 syncrefyear.interval = 1,
                                 formation.hazard.type = "agegapry",
                                 person.art.accept.threshold.dist.type = "fixed",
                                 person.art.accept.threshold.dist.fixed.value = 0.5,
                                 person.eagerness.man.type = "independent",
                                 person.eagerness.woman.type = "independent",
                                 person.eagerness.man.dist.type = "gamma",
                                 person.eagerness.woman.dist.type = "gamma",
                                 person.eagerness.man.dist.gamma.a = 0.231989836885,#0.15#0.85#0.1
                                 person.eagerness.man.dist.gamma.b = 45,#70#100#3.5#5#10#20 #170
                                 person.eagerness.woman.dist.gamma.a = 0.231989836885,#0.15#0.1
                                 person.eagerness.woman.dist.gamma.b = 45,#70#100#3.5#5#10#20#170
                                 person.agegap.man.dist.type = "normal",
                                 person.agegap.woman.dist.type = "normal",
                                 person.agegap.man.dist.normal.mu = 0, #-5
                                 person.agegap.woman.dist.normal.mu = 0, #2.5
                                 person.agegap.man.dist.normal.sigma = 1,
                                 person.agegap.woman.dist.normal.sigma = 1,
                                 formation.hazard.agegapry.numrel_man = -0.5,
                                 formation.hazard.agegapry.numrel_woman = -0.5,
                                 formation.hazard.agegapry.gap_factor_man_exp = -0.35,#-0.15#-0.5
                                 formation.hazard.agegapry.gap_factor_woman_exp = -0.35,#-0.15
                                 formation.hazard.agegapry.gap_factor_man_age = 0.05,
                                 formation.hazard.agegapry.gap_factor_woman_age = 0.05,
                                 formation.hazard.agegapry.meanage = -0.1,
                                 formation.hazard.agegapry.numrel_diff = -0.1,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_agescale_man = 0.23,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.1,#0.23
                                 formation.hazard.agegapry.eagerness_sum = 0.1,
                                 person.vsp.tofinalaids.x = 12,
                                 person.vsp.toaids.x = 7,
                                 formation.hazard.agegapry.eagerness_diff = -0.048,#-0.110975
                                 dissolution.alpha_0 = -0.52,#-0.1 # baseline
                                 dissolution.alpha_4 = -0.05,
                                 debut.debutage = 15,
                                 population.simtime = 40,
                                 population.nummen = 500, #1000 #2000
                                 population.numwomen = 500, #1000 #2000
                                 population.eyecap.fraction = 0.2,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 30,
                                 hivseed.time = 10,
                                 hivtransmission.param.a = -1.0352239,
                                 hivtransmission.param.b = -89.339994,
                                 hivtransmission.param.c = 0.4948478,
                                 hivtransmission.param.f1 = log(5), #~1.6 =>hazard is x5 in 15yo
                                 hivtransmission.param.f2 = log(log(2.5) / log(5)) / 5,
                                 mortality.aids.survtime.C = 62,
                                 mortality.aids.survtime.k = -0.2,
                                 conception.alpha_base = -2.35, #-3
                                 diagnosis.baseline = -100,
                                 monitoring.cd4.threshold = 0.1,#Treatment will not start before schedule
                                 monitoring.fraction.log_viralload = 0.3,
                                 population.msm = "yes",
                                 person.eagerness.man.msm.dist.type = "fixed",
                                 person.eagerness.man.msm.dist.fixed.value = 0,
                                 formationmsm.hazard.type = "simple",
                                 formationmsm.hazard.simple.alpha_0 = 2,
                                 formationmsm.hazard.simple.alpha_12 = -0.4,
                                 formationmsm.hazard.simple.alpha_5 = -0.2, # The factor (alpha) 5 controls the relative importance of the age gap between the partners.
                                 formationmsm.hazard.simple.alpha_6 = 0,  # weight for sum of eagerness parameters
                                 formationmsm.hazard.simple.alpha_7 = 0,
                                 birth.pregnancyduration.dist.type = "fixed",
                                 birth.pregnancyduration.dist.fixed.value = 268/365, # just over 38 weeks
                                 birth.boygirlratio = 1.0/2.01, #0.5024876, #101:100
                                 simulation.type = "simpact-cyan",
                                 monitoring.cd4.threshold.prestudy = 350,
                                 monitoring.cd4.threshold.instudy.controlstage = 350
                                 #monitoring.cd4.threshold.instudy.transitionstage = Inf,
                                 #monitoring.cd4.threshold.instudy.interventionstage = Inf
                                 ){
  input.params <- list()
  input.params$mortality.normal.weibull.shape <- mortality.normal.weibull.shape
  input.params$mortality.normal.weibull.scale <- mortality.normal.weibull.scale
  input.params$mortality.normal.weibull.genderdiff <- mortality.normal.weibull.genderdiff
  input.params$periodiclogging.interval <- periodiclogging.interval
  input.params$syncrefyear.interval <- syncrefyear.interval
  input.params$formation.hazard.type <- formation.hazard.type
  input.params$person.art.accept.threshold.dist.type <- person.art.accept.threshold.dist.type
  input.params$person.art.accept.threshold.dist.fixed.value <- person.art.accept.threshold.dist.fixed.value
  input.params$person.eagerness.man.type <- person.eagerness.man.type
  input.params$person.eagerness.woman.type <- person.eagerness.woman.type
  input.params$person.eagerness.man.dist.type <- person.eagerness.man.dist.type
  input.params$person.eagerness.woman.dist.type <- person.eagerness.woman.dist.type
  input.params$person.eagerness.man.dist.gamma.a <- person.eagerness.man.dist.gamma.a
  input.params$person.eagerness.man.dist.gamma.b <- person.eagerness.man.dist.gamma.b
  input.params$person.eagerness.woman.dist.gamma.a <- person.eagerness.woman.dist.gamma.a
  input.params$person.eagerness.woman.dist.gamma.b <- person.eagerness.woman.dist.gamma.b
  input.params$person.agegap.man.dist.type <- person.agegap.man.dist.type
  input.params$person.agegap.woman.dist.type <- person.agegap.woman.dist.type
  input.params$person.agegap.man.dist.normal.mu <- person.agegap.man.dist.normal.mu
  input.params$person.agegap.woman.dist.normal.mu <- person.agegap.woman.dist.normal.mu
  input.params$person.agegap.man.dist.normal.sigma <- person.agegap.man.dist.normal.sigma
  input.params$person.agegap.woman.dist.normal.sigma <- person.agegap.woman.dist.normal.sigma
  input.params$formation.hazard.agegapry.numrel_man <- formation.hazard.agegapry.numrel_man
  input.params$formation.hazard.agegapry.numrel_diff <- formation.hazard.agegapry.numrel_diff
  input.params$formation.hazard.agegapry.numrel_woman <- formation.hazard.agegapry.numrel_woman
  input.params$formation.hazard.agegapry.gap_factor_man_age <- formation.hazard.agegapry.gap_factor_man_age
  input.params$formation.hazard.agegapry.gap_factor_woman_age <- formation.hazard.agegapry.gap_factor_woman_age
  input.params$formation.hazard.agegapry.meanage <- formation.hazard.agegapry.meanage
  input.params$formation.hazard.agegapry.gap_factor_man_exp <- formation.hazard.agegapry.gap_factor_man_exp
  input.params$formation.hazard.agegapry.gap_factor_woman_exp <- formation.hazard.agegapry.gap_factor_woman_exp
  input.params$formation.hazard.agegapry.gap_factor_man_const <- formation.hazard.agegapry.gap_factor_man_const
  input.params$formation.hazard.agegapry.gap_factor_woman_const <- formation.hazard.agegapry.gap_factor_woman_const
  input.params$formation.hazard.agegapry.gap_agescale_man <- formation.hazard.agegapry.gap_agescale_man
  input.params$formation.hazard.agegapry.gap_agescale_woman <- formation.hazard.agegapry.gap_agescale_woman
  input.params$formation.hazard.agegapry.eagerness_sum <- formation.hazard.agegapry.eagerness_sum
  input.params$formation.hazard.agegapry.eagerness_diff <- formation.hazard.agegapry.eagerness_diff
  input.params$dissolution.alpha_0 <- dissolution.alpha_0
  input.params$dissolution.alpha_4 <- dissolution.alpha_4
  input.params$debut.debutage <- debut.debutage
  input.params$population.simtime <- population.simtime
  input.params$population.nummen <- population.nummen
  input.params$population.numwomen <- population.numwomen
  input.params$population.maxevents <-  population.simtime * population.nummen * 5 #4events/person
  input.params$population.eyecap.fraction <- population.eyecap.fraction
  input.params$hivseed.type <- hivseed.type
  input.params$hivseed.amount <- hivseed.amount
  input.params$hivseed.age.min <- hivseed.age.min
  input.params$hivseed.age.max <- hivseed.age.max
  input.params$hivseed.time <- hivseed.time
  input.params$person.vsp.tofinalaids.x <- person.vsp.tofinalaids.x
  input.params$person.vsp.toaids.x <- person.vsp.toaids.x
  input.params$hivtransmission.param.a <- hivtransmission.param.a
  input.params$hivtransmission.param.b <- hivtransmission.param.b
  input.params$hivtransmission.param.c <- hivtransmission.param.c
  input.params$hivtransmission.param.f1 <- hivtransmission.param.f1
  input.params$hivtransmission.param.f2 <- hivtransmission.param.f2
  input.params$mortality.aids.survtime.k <- mortality.aids.survtime.k
  input.params$mortality.aids.survtime.C <- mortality.aids.survtime.C
  input.params$conception.alpha_base <- conception.alpha_base
  input.params$monitoring.fraction.log_viralload <- monitoring.fraction.log_viralload
  input.params$diagnosis.baseline <- diagnosis.baseline
  input.params$birth.boygirlratio <- birth.boygirlratio
  input.params$population.msm <- population.msm
  input.params$person.eagerness.man.msm.dist.type <- person.eagerness.man.msm.dist.type
  input.params$person.eagerness.man.msm.dist.fixed.value <- person.eagerness.man.msm.dist.fixed.value
  input.params$formationmsm.hazard.type <- formationmsm.hazard.type
  input.params$formationmsm.hazard.simple.alpha_0 <- formationmsm.hazard.simple.alpha_0
  input.params$formationmsm.hazard.simple.alpha_12 <- formationmsm.hazard.simple.alpha_12
  input.params$formationmsm.hazard.simple.alpha_5 <- formationmsm.hazard.simple.alpha_5 # The factor (alpha)5 controls the relative importance of the age gap between the partners.
  input.params$formationmsm.hazard.simple.alpha_6 <- formationmsm.hazard.simple.alpha_6  # weight for sum of eagerness parameters
  input.params$formationmsm.hazard.simple.alpha_7 <- formationmsm.hazard.simple.alpha_7
  input.params$birth.pregnancyduration.dist.type <- birth.pregnancyduration.dist.type
  input.params$birth.pregnancyduration.dist.fixed.value <- birth.pregnancyduration.dist.fixed.value

  if(simulation.type == "simpact-cyan"){
    input.params$monitoring.cd4.threshold <- monitoring.cd4.threshold
  }else{
    input.params$facilities.outfile.facilityxypos <- "${SIMPACT_OUTPUT_PREFIX}facilitypositions.csv"
    input.params$monitoring.cd4.threshold.prestudy <- 350
    input.params$monitoring.cd4.threshold.instudy.controlstage <- 350
    input.params$monitoring.cd4.threshold.poststudy <- 20000 #post study everyone was now on test and study
    #input.params$monitoring.cd4.threshold.instudy.transitionstage <- Inf
    #input.params$monitoring.cd4.threshold.instudy.interventionstage <- Inf
  }

  return(input.params)
}
