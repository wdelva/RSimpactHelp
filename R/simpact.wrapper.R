simpact.wrapper <- function(inputvector = input.vector){
  age.distr <- agedistr.creator(shape = 5, scale = 65)


  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 20,
                                   population.nummen = 250,
                                   population.numwomen = 250,
                                   hivseed.time = 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 25, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 25,
                                   hivtransmission.param.a = -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5)
  seedid <- inputvector[1]
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[2]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[2]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[3]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[3]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[4]
  identifier <- paste0(seedid)
  destDir <- paste0("$VSC_SCRATCH_VO_USER/agemixing/temp/", # for VSC
                    identifier)

  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = destDir,
                                  agedist = age.distr,
                                  seed = seedid,
                                  identifierFormat = identifier),
                      error = simpact.errFunction)
  if (length(results) == 0){
    outputvector <- rep(NA, 6)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 6)
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
                                               method = "lme"),
                                  error = function(agemixing.err) {
                                    return(list()) # Returns an empty list if the lme(r) models can't be fitted
                                  })
      AAD.10.male <- ifelse(length(agemix.10.model) > 0, mean(dplyr::filter(agemix.10.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      SDAD.10.male <- ifelse(length(agemix.10.model) > 0, sd(dplyr::filter(agemix.10.rels.like.SHIMS.df, Gender =="male")$AgeGap), bignumber)
      slope.10.male <- ifelse(length(agemix.10.model) > 0, as.numeric(summary(agemix.10.model$male.model)$coefficients$fixed[2]), bignumber)
      WSD.10.male <- ifelse(length(agemix.10.model) > 0, summary(agemix.10.model$male.model)$sigma, bignumber)
      BSD.10.male <- ifelse(length(agemix.10.model) > 0, bsd.extractor(agemix.10.model, gender = "male"), bignumber)
      intercept.10.male <- ifelse(length(agemix.10.model) > 0, as.numeric(summary(agemix.10.model$male.model)$coefficients$fixed[1]), bignumber)
      outputvector <- c(AAD.10.male, SDAD.10.male, slope.10.male, WSD.10.male, BSD.10.male, intercept.10.male)
    }
  }
  unlink(paste0(destDir, "/*"), recursive = FALSE)
  return(outputvector)
}


# av AD, sd AD, slope, wsd, bsd, intercept, mean degree, ccp, relationtime, fractimeinrels, exp(pop.growth)
#am.features.1a <- c(4.66, 5.39, 0.39, 2.94, 2.19, 16.42)

#lls = c(0.1, -3, 0.1)
#uls = c(0.9, 0, 5.0)
