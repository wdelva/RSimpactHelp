#I presume this will be called through source()
#Idea is to give the simpact parameters and their ranges
# to be calibrated during the simpact run

source.simpact.parameters <- function(init.design.points = 3, resample = 1){

  ###### Generate the input parameters for the simulation #################################
  simpact.input.par <- simpact.config.inputs(design.points = init.design.points, resample.count = 1,
                                             conception.alpha_base = c(-5, -0.1),
                                             person.art.accept.threshold.dist.fixed.value = c(0.4, 0.9),
                                             person.agegap.man.dist.normal.sigma = c(0.5, 3.5),
                                             person.agegap.woman.dist.normal.sigma = c(0.5, 3.5),
                                             hivtransmission.param.f1 = c(log(2), log(3.5))
  )
  return(simpact.input.par)

}

############### Argument sample above ##################################################

