#I presume this will be called through source()
#Idea is to give the simpact parameters and their ranges
# to be calibrated during the simpact run


source.simpact.parameters <- function(init.design.points = 3, resample = 1){


###### Generate the input parameters for the simulation #################################
  simpact.input.par <- simpact.config.inputs(design.points = init.design.points, resample.count = 1,
                                       conception.alpha_base = c(-5, -0.1),
                                       person.art.accept.threshold.dist.fixed.value = c(0.4, 1),
                                       person.eagerness.man.dist.gamma.a = c(0.3, 1.5),
                                       person.eagerness.man.dist.gamma.b = c(10, 60),
                                       person.eagerness.woman.dist.gamma.a = c(0.3, 1.5),
                                       formation.hazard.agegapry.eagerness_diff =  c(-0.1, 0),
                                       person.eagerness.woman.dist.gamma.b = c(10, 60),
                                       formation.hazard.agegapry.numrel_man = c(-2, -0.1),
                                       formation.hazard.agegapry.numrel_woman = c(-2, -0.1),
                                       formation.hazard.agegapry.gap_factor_man_exp = c(-2, -0.01),
                                       formation.hazard.agegapry.gap_factor_woman_exp = c(-2, -0.01),
                                       person.agegap.man.dist.normal.mu = c(1, 5),
                                       person.agegap.woman.dist.normal.mu = c(1, 5),
                                       person.agegap.man.dist.normal.sigma = c(0.5, 3.5),
                                       person.agegap.woman.dist.normal.sigma = c(0.5, 3.5),
                                       hivtransmission.param.f1 = c(log(2), log(3.5))
  )

    return(simpact.input.par)

}

############### Argument sample above ##################################################

