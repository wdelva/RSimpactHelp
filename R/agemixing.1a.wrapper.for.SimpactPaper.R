#' Wrapper function for age mixing example of Simpact paper.
#'
#' This function simulates an HIV epidemic for 20 years and produces model
#' features for the agemixing pattern, the sexual behaviour, the population
#' growth rate, HIV prevalence, HIV incidence, and ART coverage.
#'
#' @param inputvector Vector of model input parameter values
#' @return A vector of model features
#' @import RSimpactCyan
#' @importFrom fitdistrplus fitdist
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export

agemixing.1a.wrapper.for.SimpactPaper <- function(inputvector = input.vector){
  #destDir <- "/Users/delvaw/Documents/temp"                     # on laptop
  destDir <- "/user/data/gent/vsc400/vsc40070/agemixing/temp" # for VSC
  age.distr <- agedistr.creator(shape = 5, scale = 65)


  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 10,
                                   population.nummen = 2500,
                                   population.numwomen = 2500,
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 25, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 25,
                                   hivtransmission.param.a = -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
                                   hivtransmission.param.f2 = log(log(sqrt(2)) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1,
                                   formation.hazard.agegapry.gap_agescale_man = 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = 0.25,
                                   debut.debutage = 15,
                                   conception.alpha_base = -2.5
  )


  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["dropout.interval.dist.uniform.min"] <- 100
  cfg.list["dropout.interval.dist.uniform.max"] <- 200

  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

  cfg.list["person.agegap.man.dist.type"] <- "normal"
  cfg.list["person.agegap.woman.dist.type"] <- "normal"

  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0
  cfg.list["diagnosis.baseline"] <- -99
  cfg.list["population.msm"] <- "no"
  cfg.list["hsv2seed.time"] <- -1 # No HVS2 seeding event
  cfg.list["hsv2seed.time"] <- -1
  cfg.list["hsv2seed.type"] <- "amount"



  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[2]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[2]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[3]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[3]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[4]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[5]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[7]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[8]

  cfg <- cfg.list

  cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 8
  cfg["monitoring.fraction.log_viralload"] <- 0.3
  cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

  seedid <- inputvector[1]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[9]
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[9]
  cfg["formation.hazard.agegapry.baseline"] <- inputvector[10]

  cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[11]
  cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[12]
  cfg["conception.alpha_base"] <- inputvector[13]
  cfg["dissolution.alpha_0"] <- inputvector[14]
  cfg["dissolution.alpha_4"] <- 0



  # # Here we insert the ART acceptability paramter and the ART interventions, so that we can calibrate them
  # cfg["person.art.accept.threshold.dist.fixed.value"] <- 0.3
  #
  # art.intro <- list()
  # art.intro["time"] <- 25
  # art.intro["diagnosis.baseline"] <- -1.2
  # art.intro["monitoring.cd4.threshold"] <- 100
  #
  # # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  # art.intro2 <- list()
  # art.intro2["time"] <- 25 + 5
  # art.intro2["monitoring.cd4.threshold"] <- 200
  #
  # art.intro3 <- list()
  # art.intro3["time"] <- 25 + 8
  # art.intro3["monitoring.cd4.threshold"] <- 350
  #
  # art.intro4 <- list()
  # art.intro4["time"] <- 25 + 11
  # art.intro4["monitoring.cd4.threshold"] <- 500
  #
  # art.intro5 <- list()
  # art.intro5["time"] <- 25 + 13 #25 + 13
  # art.intro5["monitoring.cd4.threshold"] <- 600
  #
  # intervention <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)


  results <- tryCatch(simpact.run(configParams = cfg,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  seed = seedid),
                      error = simpact.errFunction)

  if (length(results) == 0){
    outputvector <- rep(NA, 11)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 11)
    } else {
      datalist.agemix <- readthedata(results)
      agemix.episodes.df <- agemix.episodes.df.maker(datalist.agemix)
      agemix.rels.like.SHIMS.df <- agemix.rels.df.maker(dataframe = agemix.episodes.df,
                                                        agegroup = c(18, 50),
                                                        timepoint = datalist.agemix$itable$population.simtime[1],
                                                        timewindow = 0.5, # Ongoing the 3 most recent relationships in the past 6 months
                                                        start = FALSE) %>%
        dplyr::group_by(ID) %>%
        dplyr::top_n(3, DisTime)

      bignumber <- NA

      # agemix.model <- tryCatch(amp.modeller(dataframe = agemix.episodes.df,
      #                                       agegroup = c(18, 50),
      #                                       timepoint = datalist.agemix$itable$population.simtime[1],
      #                                       timewindow = 0.5,
      #                                       start = FALSE,
      #                                       SHIMS = TRUE,
      #                                       method = "lmer"),
      #                          error = function(agemixing.err) {
      #                            return(list()) # Returns an empty list if the lme(r) models can't be fitted
      #                          })
      #
      #
      # AAD.male <- ifelse(length(agemix.model) > 0, mean(dplyr::filter(agemix.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      # SDAD.male <- ifelse(length(agemix.model) > 0, sd(dplyr::filter(agemix.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      # slope.male <- ifelse(length(agemix.model) > 0, summary(agemix.model$male.model)$coefficients[2, 1], bignumber)
      # WSD.male <- ifelse(length(agemix.model) > 0, summary(agemix.model$male.model)$sigma, bignumber)
      # BSD.male <- ifelse(length(agemix.model) > 0, bsd.extractor(agemix.model, gender = "male"), bignumber)
      # intercept.male <- ifelse(length(agemix.model) > 0, summary(agemix.model$male.model)$coefficients[1,1], bignumber)

      # The above agemixing metrics are the same as those computed in Step 1a.
      ###########
      # We also want to calculate them for the period around the introduction of HIV (t9.5 - t10).
      agemix.10.rels.like.SHIMS.df <- agemix.rels.df.maker(dataframe = agemix.episodes.df,
                                                           agegroup = c(18, 50),
                                                           timepoint = datalist.agemix$itable$hivseed.time[1],
                                                           timewindow = 0.5, # Ongoing the 3 most recent relationships in the past 6 months
                                                           start = FALSE) %>%
        dplyr::group_by(ID) %>%
        dplyr::top_n(3, DisTime)
      agemix.10.model <- tryCatch(amp.modeller(dataframe = agemix.episodes.df,
                                               agegroup = c(18, 50),
                                               timepoint = datalist.agemix$itable$hivseed.time[1],
                                               timewindow = 0.5,
                                               start = FALSE,
                                               SHIMS = TRUE,
                                               method = "lmer"),
                                  error = function(agemixing.err) {
                                    return(list()) # Returns an empty list if the lme(r) models can't be fitted
                                  })
      AAD.10.male <- ifelse(length(agemix.10.model) > 0, mean(dplyr::filter(agemix.10.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      SDAD.10.male <- ifelse(length(agemix.10.model) > 0, sd(dplyr::filter(agemix.10.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      slope.10.male <- ifelse(length(agemix.10.model) > 0, summary(agemix.10.model$male.model)$coefficients[2, 1], bignumber)
      WSD.10.male <- ifelse(length(agemix.10.model) > 0, summary(agemix.10.model$male.model)$sigma, bignumber)
      BSD.10.male <- ifelse(length(agemix.10.model) > 0, bsd.extractor(agemix.10.model, gender = "male"), bignumber)
      intercept.10.male <- ifelse(length(agemix.10.model) > 0, summary(agemix.10.model$male.model)$coefficients[1,1], bignumber)

      # Sexual behaviour metrics
      # num.rels <- agemix.rels.df.maker(dataframe = agemix.episodes.df,
      #                                  agegroup = c(18, 50),
      #                                  timepoint = datalist.agemix$itable$population.simtime[1],
      #                                  timewindow = 1,
      #                                  start = TRUE) %>%
      #   dplyr::filter(Gender == "male") %>%
      #   nrow()
      # num.men <- alive.infected(datalist = datalist.agemix,
      #                           timepoint = datalist.agemix$itable$population.simtime[1],
      #                           site = "All") %>%
      #   dplyr::filter(Gender == 1,
      #                 TOB <= datalist.agemix$itable$population.simtime[1] - 18,
      #                 TOB > datalist.agemix$itable$population.simtime[1] - 50) %>%
      #   nrow()
      # meandegree.male <- num.rels/num.men
      # # # Concurrency point prevalence 6 months before a survey, among men
      # pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
      #                                                    agegroup = c(15, 50),
      #                                                    timepoint = datalist.agemix$itable$population.simtime[1],
      #                                                    hivstatus = 2)[1,2] %>% as.numeric()

      # The above metrics are the same as the one computed in Step 1a.
      ###########
      # We also want to calculate them for the period around the introduction of HIV.
      num.rels.10 <- agemix.rels.df.maker(dataframe = agemix.episodes.df,
                                          agegroup = c(18, 50),
                                          timepoint = datalist.agemix$itable$hivseed.time[1],
                                          timewindow = 1,
                                          start = TRUE) %>%
        dplyr::filter(Gender == "male") %>%
        nrow()
      num.men.10 <- alive.infected(datalist = datalist.agemix,
                                   timepoint = datalist.agemix$itable$hivseed.time[1],
                                   site = "All") %>%
        dplyr::filter(Gender == 1,
                      TOB <= datalist.agemix$itable$hivseed.time[1] - 18,
                      TOB > datalist.agemix$itable$hivseed.time[1] - 50) %>%
        nrow()
      meandegree.male.10 <- num.rels.10/num.men.10
      # # Concurrency point prevalence 6 months before a survey, among men
      pp.cp.6months.male.10 <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                            agegroup = c(15, 50),
                                                            timepoint = datalist.agemix$itable$hivseed.time[1],
                                                            hivstatus = 2)[1,2] %>% as.numeric()

      # AMOUNT OF TIME SPENT IN RELATIONSHIPS in the first 10 years of THE SIMULATION: A NEW METRIC IN THIS EXPLORATION
      rel.start.vector.10 <- pmin(datalist.agemix$rtable$FormTime, 10)
      rel.end.vector.10 <- pmin(datalist.agemix$rtable$DisTime, 10)

      total.rel.time.10 <- sum(rel.end.vector.10 - rel.start.vector.10)

      # Also the ratio of relationship time in the last year of the simulation,
      # divided by number of people alive at the timepoint (population.simtime - 1),
      # to get a sense for the rate at which people for relationships (because the rate at which they break up is the same in all simulations)

      rel.end.minus.one.vector.10 <- pmin(datalist.agemix$rtable$DisTime, 9)#(datalist.agemix$itable$population.simtime[1] - 1))
      rel.start.minus.one.vector.10 <- pmin(datalist.agemix$rtable$FormTime, 9)#(datalist.agemix$itable$population.simtime[1] - 1))
      total.rel.time.minus.one.10 <- sum(rel.end.minus.one.vector.10 - rel.start.minus.one.vector.10)
      rel.time.minus.one.10 <- total.rel.time.10 - total.rel.time.minus.one.10 # reltime in the last year before HIV is introduced
      pop.size.minus.one.10 <- datalist.agemix$ltable$PopSize[(datalist.agemix$itable$hivseed.time[1] - 1)]
      rel.pp.minus.one.10 <- rel.time.minus.one.10 / pop.size.minus.one.10


      # #########
      # # The HIV prevalence, over time
      # #########
      # prev.vector <- tryCatch(prevalence.vector.creator(datalist = datalist.agemix, agegroup = c(15, 50)) %>%
      #                           dplyr::filter(pointestimate == TRUE) %>%
      #                           dplyr::select(prevalence) %>%
      #                           unlist() %>%
      #                           as.numeric(),
      #                         error = function(prev.err) {
      #                           return(rep(NA, (1 + datalist.agemix$itable$population.simtime[1])))
      #                         })
      # ########
      # # The HIV incidence, over time
      # #########
      # inc.vector <- tryCatch(incidence.vector.creator(datalist = datalist.agemix, agegroup = c(15, 50)),
      #                        error = function(inc.err) {
      #                          return(rep(NA, datalist.agemix$itable$population.simtime[1]))
      #                        })
      #
      # ########
      # # ART coverage, over time
      # #########
      # art.cov.vector <- tryCatch(ART.coverage.vector.creator(datalist = datalist.agemix, agegroup = c(15, 50)),
      #                            error = function(art.cov.err) {
      #                              return(rep(NA, (1 + datalist.agemix$itable$population.simtime[1])))
      #                            })
      #
      # #########
      # # The average age of HIV-positive people, over time
      # #########
      # mean.age.pos.vector <- tryCatch(mean.age.pos.vector.creator(datalist = datalist.agemix,
      #                                                             timevector = seq(from = datalist.agemix$itable$hivseed.time[1],
      #                                                                              to = datalist.agemix$itable$population.simtime[1])),
      #                                 error = function(mean.pos.err) {
      #                                   return(rep(NA, (1 + datalist.agemix$itable$population.simtime[1] - datalist.agemix$itable$hivseed.time[1])))
      #                                 })
      #
      # # Population growth rate
      # growthrate <- pop.growth.calculator(datalist = datalist.agemix,
      #                                     timewindow = c(0, datalist.agemix$itable$population.simtime[1]))
      growthrate.10 <- pop.growth.calculator(datalist = datalist.agemix,
                                             timewindow = c(0, datalist.agemix$itable$hivseed.time[1]))

      outputvector <- c(#AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                        AAD.10.male, SDAD.10.male, slope.10.male, WSD.10.male, BSD.10.male, intercept.10.male,
                        #meandegree.male,
                        meandegree.male.10,
                        #pp.cp.6months.male,
                        pp.cp.6months.male.10,
                        total.rel.time.10,
                        rel.pp.minus.one.10,
                        #exp(growthrate),
                        exp(growthrate.10))#,
                        #prev.vector,
                        #inc.vector,
                        #art.cov.vector,
                        #mean.age.pos.vector)
    }
  }
  return(outputvector)
}
