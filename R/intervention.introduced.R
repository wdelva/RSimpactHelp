#' We need to decide if we need to include all the intervention parameters here? log viralload, etc
#'
#' Introduce the intervention years and monitoring time for the cd4. Note that
#' art.intro list should include the time when the diagnosis to take effect.
#'
#'
#' The function take in a list of paramenters set at intervals of cd4 count threshhold used to
#' decide if a person can be offerd treatment
#'
#' @param time is the simulation time at which the modified setting should be introduced
#' @param diagnosis.baseline is the initial setting at which the initial cd4 threshhold is set as a baseline value.
#' @param monitoring.cd4.threshold is the cd4 value used to decide if teh person is offered treatment
#' @param diagnosis.genderfactor this allows the hazard to take into account gender of the person
#' @return a list of lists containing intervention at different times. The first list sets the diagnosis.baseline value.
#' @examples
#' intervention.introduced <- intervention.introduced(list(27,0,100,2),list(30,200,1.5), list(33,350,1),list(36,500,0.5))


intervention.introduced <- function(art.intro = list(time = 27,diagnosis.baseline = 0,
                                                     monitoring.cd4.threshold = 100,
                                                     diagnosis.genderfactor = 2),
                                    art.intro2 = list(time = 30, monitoring.cd4.threshold = 200,
                                                      diagnosis.genderfactor = 1.5), ...){

  arg.list <- list(art.intro, art.intro2,...)

  iv <- list()

  try(
    if(length(arg.list[[1]])!=4) {
      stop("Check if Diagnosis.baseline is set")
      }
      else{
    for (i in 1:length(arg.list)){
      #print(length(arg.list[[1]]))
      if(i==1){
        art.intro <- list()
        art.intro["time"] <- arg.list[[i]][1]
        art.intro["diagnosis.baseline"] <- arg.list[[i]][2]
        art.intro["monitoring.cd4.threshold"] <- arg.list[[i]][3]
        art.intro["diagnosis.genderfactor"] <- arg.list[[i]][4]
      }else{
      art.intro <- list()
      art.intro["time"] <- arg.list[[i]][1]
      art.intro["monitoring.cd4.threshold"] <- arg.list[[i]][2]
      art.intro["diagnosis.genderfactor"] <- arg.list[[i]][3]
      }

      iv[[length(iv)+1]] <- art.intro

    }

    return(iv)
  })
}
