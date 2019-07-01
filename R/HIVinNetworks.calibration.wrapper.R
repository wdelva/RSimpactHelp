#' Wrapper function for running simpact simulations for the calibration step in the HIV in Networks project
#'
#'
#' @param index index
#' @param list_param list of vectors of random number seed and parameter values
#' @return A vector of model features (summary statistics of simulation output)
#' @import RSimpactCyan
#' @import dplyr
#' @importFrom MASS glm.nb
#' @importFrom magrittr %>%
#' @export

HIVinNetworks.calibration.wrapper <- function(index, list_param){
  inputvector <- list_param[[index]]

  age.distr <- agedistr.creator(shape = 5, scale = 65)


  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 10,
                                   population.nummen = 2500,
                                   population.numwomen = 2500,
                                   population.msm = "no",
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 20,
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
                                   conception.alpha_base = -2.7,
                                   dropout.interval.dist.type = "exponential")

  #standard deviation of 200 CD4 cells
  #mu = ln(mean / sqrt(1 + variance/mean^2))
  #sigma^2 = ln(1 + variance/mean^2)
  #Here, we say mean = 825 and variance = 200^2
  mu.cd4 <- 800
  var.cd4 <- 200^2
  mu.cd4.end <- 20
  var.cd4.end <- 5
  cfg.list["person.cd4.start.dist.type"] <- "lognormal"
  cfg.list["person.cd4.start.dist.lognormal.zeta"] <- log(mu.cd4/sqrt(1+var.cd4/mu.cd4^2))
  cfg.list["person.cd4.start.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4/mu.cd4^2))
  cfg.list["person.cd4.end.dist.type"] <- "lognormal"
  cfg.list["person.cd4.end.dist.lognormal.zeta"] <- log(mu.cd4.end/sqrt(1+var.cd4.end/mu.cd4.end^2))
  cfg.list["person.cd4.end.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4.end/mu.cd4.end^2))

  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3 # 0.3

  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

  cfg.list["person.agegap.man.dist.type"] <- "normal"
  cfg.list["person.agegap.woman.dist.type"] <- "normal"

  cfg.list["monitoring.cd4.threshold"] <- 1 # 0
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
  cfg.list["diagnosis.baseline"] <- -99999 # -2
  cfg.list["periodiclogging.interval"] <- 0.25
  cfg.list["dropout.interval.dist.exponential.lambda"] <- 0.1

  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 6

  cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

  seedid <- inputvector[1]
  cfg.list["hivtransmission.param.f1"] = log(1.1)
  cfg.list["hivtransmission.param.f2"] = log(log(sqrt(1.1)) / log(1.1)) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = 0.25
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = 0.25
  cfg.list["person.agegap.man.dist.normal.mu"] <- -1
  cfg.list["person.agegap.woman.dist.normal.mu"] <- -1
  cfg.list["person.agegap.man.dist.normal.sigma"] <- 4
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- 4
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[2]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[3]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[4]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[5]
  cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- -0.4
  cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- -0.4
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[6]
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[7]
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8]
  cfg.list["conception.alpha_base"] <- inputvector[9]
  cfg.list["dissolution.alpha_0"] <- inputvector[10]



  identifier <- paste0(seedid)
  rootDir <- "/tmp" # "/user/scratch/gent/vsc400/vsc40070/EAAA/Fa/temp"

  destDir <- paste0(rootDir, "/", identifier)


  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  seed = seedid,
                                  identifierFormat = identifier),
                      error = simpact.errFunction)
  if (length(results) == 0){
    outputvector <- rep(NA, 4)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 4)
    } else {
      datalist.networks <- readthedata(results)
      agemix.df <- agemix.df.maker(datalist.networks)

      agegroup <- c(18, 50)
      ####
      # Population growth rate
      ####
      growthrate <- pop.growth.calculator(datalist = datalist.networks,
                                          timewindow = c(0, 10))
      ####
      # Point prevalence of concurrent relationships
      ####
      ppconc.tibble <- concurr.pointprev.calculator(datalist = datalist.networks,
                                             agegroup = agegroup,
                                             timepoint = 10,
                                             hivstatus = 2)
      ppconc <- ppconc.tibble$concurr.pointprev[3]

      ####
      # Partner turnover rate
      ####
      degrees.df <- degree.df.maker(df = agemix.df,
                                    agegroup = agegroup,
                                    hivstatus = 2,
                                    survey.time = 10,
                                    window.width = 1,
                                    gender.degree = "male",
                                    only.new = TRUE)


      # If we want to know % of men who had 0 partners in the past 12 months,
      # We need to compare nrow(degree.df) with the number of 18-50 year old men
      # that were alive at the time of the survey
      allsurveymen <- dplyr::filter(datalist.networks$ptable,
                                    Gender == 1, # Male
                                    TOD > 10, # Still alive at the time of the survey
                                    TOB <= 10 - agegroup[1], # Not younger than 18 at the time of the survey
                                    TOB > 10 - agegroup[2])#, # Not older than 50 at the time of the survey

      # Creating the vector of degrees
      degree.vector <- tryCatch(c(rep(0, times = (nrow(allsurveymen) - nrow(degrees.df))), degrees.df$Degree),
                                error = function(degree.vector.err) {
                                  return(NA)
                                })
      meandegree.male <- mean(degree.vector)
      # Fitting the negative binomial distribution to this vector
      mean.nb.male <- NA # Will be overwritten if a negbin distribution can be fitted to the data
      size.nb.male <- NA # Will be overwritten if a negbin distribution can be fitted to the data
      allpartnerships.simulated <- tryCatch(glm.nb(degree.vector ~ 1),
                                                                        error = function(glm.err) {
                                                                          return(NA)
                                                                        })
      if ("negbin" %in% class(allpartnerships.simulated)){
        output.simulated <- summary(allpartnerships.simulated)
        mean.nb.male <- exp(output.simulated$coef[1]) #(mu in dnbinom)
        size.nb.male <- output.simulated$theta  #(size in dnbinom: number of successful trials)
      }

      outputvector <- c(exp(growthrate),
                        exp(ppconc),
                        mean.nb.male,
                        size.nb.male)
    }
  }

  unlink(paste0(rootDir, "/", identifier), recursive = TRUE)
  return(outputvector)
}
