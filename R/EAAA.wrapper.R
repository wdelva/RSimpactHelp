#' Wrapper function for running simpact simulations for the MaxART EAAA simulation study
#'
#' A short description here...
#'
#' @param inputvector Vector of random seed and parameter values
#' @return A vector of model features (summary statistics of simulation output)
#' @import RSimpactCyan
#' @import dplyr
#' @importFrom magrittr %>%
#' @export

EAAA.wrapper <- function(inputvector = input.vector){
  age.distr <- agedistr.creator(shape = 5, scale = 65)


  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 36.75, #40, #Until 1 October 2016
                                   population.nummen = 500,
                                   population.numwomen = 500,
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 20, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   hivtransmission.param.a = -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(2),
                                   hivtransmission.param.f2 = log(log(sqrt(2)) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1,
                                   formation.hazard.agegapry.gap_agescale_man = 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = 0.25,
                                   dissolution.alpha_4 = -0.05,
                                   debut.debutage = 15,
                                   conception.alpha_base = -2.7)

  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3 # 0.3
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000

  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

  cfg.list["person.agegap.man.dist.type"] <- "normal"
  cfg.list["person.agegap.woman.dist.type"] <- "normal"

  cfg.list["monitoring.cd4.threshold"] <- 0
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
  cfg.list["diagnosis.baseline"] <- -2



  # Introducing ART
  art.intro <- list()
  art.intro["time"] <- 20     # ~2000
  art.intro["diagnosis.baseline"] <- -2
  art.intro["monitoring.cd4.threshold"] <- 100

  art.intro1 <- list()
  art.intro1["time"] <- 22     # ~2002
  art.intro1["diagnosis.baseline"] <- -1.8
  art.intro1["monitoring.cd4.threshold"] <- 150

  art.intro2 <- list()
  art.intro2["time"] <- 25     # ~2005
  art.intro2["diagnosis.baseline"] <- -1.5
  art.intro2["monitoring.cd4.threshold"] <- 200

  art.intro3 <- list()
  art.intro3["time"] <- 30     # ~2010
  art.intro2["diagnosis.baseline"] <- -1
  art.intro3["monitoring.cd4.threshold"] <- 350

  art.intro4 <- list()
  art.intro4["time"] <- 33     # ~2013
  art.intro4["monitoring.cd4.threshold"] <- 500

  art.intro5 <- list()
  art.intro5["time"] <- 36.75     # October ~2016
  art.intro5["monitoring.cd4.threshold"] <- 6000


  ART.factual <- list(art.intro,art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  ART.counterfactual <- list(art.intro,art.intro1, art.intro2, art.intro3)


  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3

  cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

  seedid <- inputvector[1]
  cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]
  cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10]
  cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10]
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  cfg.list["conception.alpha_base"] <- inputvector[14]
  cfg.list["dissolution.alpha_0"] <- inputvector[15]


  seedid <- inputvector[1]
  identifier <- paste0(seedid)
  destDir <- paste0("/Users/delvaw/Documents/temp/",
                    identifier) # on laptop

  #destDir <- paste0("/user/scratch/gent/vsc400/vsc40070/agemixing/temp/", # for VSC
  #                  identifier)

  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  intervention = ART.factual,
                                  seed = seedid,
                                  identifierFormat = identifier),
                      error = simpact.errFunction)
  if (length(results) == 0){
    outputvector <- rep(NA, 37)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 37)
    } else {
      datalist.EAAA <- readthedata(results)

      ####
      # Population growth rate
      ####
      growthrate <- pop.growth.calculator(datalist = datalist.EAAA,
                                          timewindow = c(20, datalist.EAAA$itable$population.simtime[1])) # Between January 2000 and October 2016
      ####
      # HIV prevalence. To be compared to SHIMS I estimates (point estimate at March 2011 ~ t = 31.25)
      ####
      #f.18.20
      prev.f.18.19 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(18, 20),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.18.19 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(18, 20),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.20.24 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(20, 25),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.20.24 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(20, 25),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.25.29 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(25, 30),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.25.29 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(25, 30),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.30.34 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(30, 35),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.30.34 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(30, 35),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.35.39 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(35, 40),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.35.39 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(35, 40),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.40.44 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(40, 45),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.40.44 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(40, 45),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      prev.f.45.49 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(45, 50),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      prev.m.45.49 <- prevalence.calculator(datalist = datalist.EAAA,
                                            agegroup = c(45, 50),
                                            timepoint = 31.25) %>%
        dplyr::select(pointprevalence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      ####
      # HIV incidence. Average follow-up period March 2011 until mid Sept 2011 (0.55 years)
      ####
      inc.f.18.19 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(18, 20),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.18.19 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(18, 20),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.20.24 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(20, 25),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.20.24 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(20, 25),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.25.29 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(25, 30),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.25.29 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(25, 30),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.30.34 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(30, 35),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.30.34 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(30, 35),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.35.39 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(35, 40),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.35.39 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(35, 40),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.40.44 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(40, 45),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.40.44 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(40, 45),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      inc.f.45.49 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(45, 50),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(2) %>%
        as.numeric()
      inc.m.45.49 <- incidence.calculator(datalist = datalist.EAAA,
                                          agegroup = c(45, 50),
                                          timewindow = c(31.25, 31.8),
                                          only.active = "No") %>%
        dplyr::select(incidence) %>%
        dplyr::slice(1) %>%
        as.numeric()
      ####
      # ART coverage among adults 15+ years old from UNAIDS (2010 - 2016 estimates)
      ####
      ART.cov.eval.timepoints <- seq(from = 30.5,
                                     to = 36.5)
      ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
      for (art.cov.index in 1:length(ART.cov.vector)){
        ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist.EAAA,
                                                                     agegroup = c(15, 150),
                                                                     timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
          sum(ART.coverage.calculator(datalist = datalist.EAAA,
                                      agegroup = c(15, 150),
                                      timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
      }
      ####
      # VL suppression fraction (all ages in 2016 ~ >= 15 yo) 0.68
      ####
      VL.suppression.fraction <- VL.suppression.calculator(datalist = datalist.EAAA,
                                                           agegroup = c(15, 300),
                                                           timepoint = 36.5,
                                                           vl.cutoff = 1000,
                                                           site="All") %>%
        dplyr::select(vl.suppr.frac) %>%
        dplyr::slice(3) %>%
        as.numeric()

      outputvector <- c(growthrate,
                        prev.f.18.19,
                        prev.m.18.19,
                        prev.f.20.24,
                        prev.m.20.24,
                        prev.f.25.29,
                        prev.m.25.29,
                        prev.f.30.34,
                        prev.m.30.34,
                        prev.f.35.39,
                        prev.m.35.39,
                        prev.f.40.44,
                        prev.m.40.44,
                        prev.f.45.49,
                        prev.m.45.49,
                        inc.f.18.19,
                        inc.m.18.19,
                        inc.f.20.24,
                        inc.m.20.24,
                        inc.f.25.29,
                        inc.m.25.29,
                        inc.f.30.34,
                        inc.m.30.34,
                        inc.f.35.39,
                        inc.m.35.39,
                        inc.f.40.44,
                        inc.m.40.44,
                        inc.f.45.49,
                        inc.m.45.49,
                        ART.cov.vector,
                        VL.suppression.fraction)
    }
  }
  unlink(paste0(destDir, "/*"), recursive = FALSE)
  return(outputvector)
}
