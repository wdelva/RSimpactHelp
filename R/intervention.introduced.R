#' We need to decide if we need to include all the intervention parameters here? log viralload, etc
#'
#' Introduce the intervention years and monitoring time for the cd4. Note that
#' art.intro list should include the time when the diagnosis to take effect.
#'
#'
#' The function take in a list of paramenters set at intervals of cd4 count threshhold used to
#' decide if a person can be offerd treatment
#' @param simulation.type An indication of which of simpact simulation to perform.
#' @return a list of lists containing intervention at different times. The first list sets the diagnosis.baseline value.
#' @examples
#' intervention.introduced <- intervention.introduced(list(27,0,100,2),list(30,200,1.5), list(33,350,1),list(36,500,0.5))

intervention.introduced <- function(simulation.type = "simpact-cyan"){

  if(simulation.type == "simpact-cyan"){
    # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
    art.intro <- list()
    art.intro["time"] <- 27
    art.intro["diagnosis.baseline"] <- 0 # Reset to zero, from its original value of -100 at the start of the simulatio
    art.intro["monitoring.cd4.threshold"] <- 100
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

    art.intro2 <- list()
    art.intro2["time"] <- 30
    art.intro2["monitoring.cd4.threshold"] <- 200
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.intro3 <- list()
    art.intro3["time"] <- 33
    art.intro3["monitoring.cd4.threshold"] <- 350
    art.intro3["diagnosis.genderfactor"] <- 1

    art.intro4 <- list()
    art.intro4["time"] <- 36
    art.intro4["monitoring.cd4.threshold"] <- 500
    art.intro4["diagnosis.genderfactor"] <- 0.5

    # person.art.accept.threshold.dist.fixed.value

    iv <- list(art.intro, art.intro2, art.intro3, art.intro4)

  }else if(simulation.type == "maxart"){
    # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
    art.intro <- list()
    art.intro["time"] <- 27
    art.intro["diagnosis.baseline"] <- 0 # Reset to zero, from its original value of -100 at the start of the simulatio
    art.intro["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2014:500

    art.intro2 <- list()
    art.intro2["time"] <- 30
    art.intro2["monitoring.cd4.threshold.prestudy"] <- 350
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.intro3 <- list()
    art.intro3["time"] <- 33
    art.intro3["monitoring.cd4.threshold.prestudy"] <- 350
    art.intro3["diagnosis.genderfactor"] <- 1

    art.intro4 <- list()
    art.intro4["time"] <- 36
    art.intro4["monitoring.cd4.threshold.prestudy"] <- 500
    art.intro4["diagnosis.genderfactor"] <- 0.5

    # person.art.accept.threshold.dist.fixed.value

    iv <- list(art.intro, art.intro2, art.intro3, art.intro4)
  }else{
    iv <-list()
  }

  return(iv)
}
