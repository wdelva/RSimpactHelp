#' We need to decide if we need to include all the intervention parameters here?
#' log viralload, etc
#'
#' Introduce the intervention years and monitoring time for the cd4. Note that
#' art.intro list should include the time when the diagnosis to take effect.
#'
#' For Uganda 1991-2008, ART introduced from 2004
#'
#' The function take in a list of paramenters set at intervals of cd4 count threshhold used to
#' decide if a person can be offerd treatment
#' @param simulation.type An indication of which of simpact simulation to perform.
#' @return a list of lists containing intervention at different times.
#' The first list sets the diagnosis.baseline value.
#' @examples
#' intervention.introduced <- intervention.introduced(simulation.type = "simpact-cyan")
#' @export

Ug.intervention.introduced <- function(simulation.type = "simpact-cyan"){

  if(simulation.type == "simpact-cyan"){
    # Simulation starts in 1970. After 34 years (in 2004), ART is introduced.
    art.intro <- list()
    art.intro["time"] <- 34
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original value @ sim start
    art.intro["monitoring.cd4.threshold"] <- 100
    art.intro["monitoring.interval.piecewise.cd4s"] <- "100,250"
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
    art.intro2 <- list()
    art.intro2["time"] <- 37
    art.intro2["monitoring.cd4.threshold"] <- 200
    art.intro["monitoring.interval.piecewise.cd4s"] <- "200,350"
    art.intro2["diagnosis.genderfactor"] <- 1.5

    # person.art.accept.threshold.dist.fixed.value

    iv <- list(art.intro, art.intro2)
  }

  return(iv)
}
