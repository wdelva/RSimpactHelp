#' Basci reproductive number
#'
#' Compute HIV basic reproductive number from a dynamic sexual network simulated with Simpact
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param beta The probability that an infected individual will infect a suceptible partner over the duration of
#' their relationship
#' @param trans.rate.int A vector of new relationship per time interval produced by \code{\link{transmission.rate.calculator}}
#' @return A value of the basic reproduction number
#'
#' @example ro <- basicnumber.calculator(datalist = datalist, beta = 0.1508, trans.rate.int = trans.rate.int)
#'
#' @export

basicnumber.calculator <- function(datalist = datalist, beta = 0.1508, trans.rate.int = trans.rate.int){

  # Ro = beta * C * D

  # C: effective mean of the average of annual number of new partners C = m + (sigma^2)/m
  # where m is the arithmetic mean of the annual number of new partners

  # D: average duration of infectiousness
  # Ideally Note: C < 1/(beta * D)

  # source("/home/david/RSimpactHelp/R/transmission.rate.calculator.R")
  #
  # trans.rate.int <- transmission.rate.calculator(datalist = datalist, timewindow = c(0,40), int = TRUE, by = 1)

  m <- mean(trans.rate.int)

  sigma <- sd(trans.rate.int)

  C <- (m + (sigma^2)/m)


  # Search D

  infec.indiv <- datalist$ptable[InfectType!=-1]

  # keep id and infectime

  infec.id.time <-subset(infec.indiv, select = c(ID,InfectTime))

  id.infec <- infec.id.time$ID

  # or
  # keeps <- c("ID","InfectTime")
  # infec.id.time <- as.data.frame(infec.indiv)[keeps]

  # person table of alive individuals
  alive.timepoint <- as.data.frame(alive.infected(datalist = datalist, timepoint = 40,
                                       site = "All"))

  # person table of alive individuals but infected
  alive.timepoint.infec <- subset(alive.timepoint, alive.timepoint$Infected=="TRUE")

  # IDs of infected by alive
  id.infec.alive <- alive.timepoint.infec$ID

  # IDs and time of infection for alive individuals
  id.time.infec.alive <- subset(infec.id.time, id.infec%in%id.infec.alive)

  # removal time for alive individuals : endpoint
  fin.time <- rep(40, nrow(id.time.infec.alive)) # this 40 is for population.simt.time

  # attach removal time to IDs and infection time
  id.time.infec.alive.fin <- cbind(id.time.infec.alive, fin.time)

  # search removal time for those who died >> id.time.infec.die

  deathtime.hiv <- datalist$etable[eventname=="aidsmortality"]

  deathtime.hiv.df <- subset(deathtime.hiv, select = c(p1ID,eventtime))

  deathtime.norm.raw <- datalist$etable[eventname=="normalmortality"] # !!! but death norm for infected individuals

  deathtime.hiv.norm <- subset(deathtime.norm.raw, deathtime.norm.raw$p1ID%in%id.infec)

  deathtime.norm.df <- subset(deathtime.hiv.norm, select = c(p1ID,eventtime))

  deathtime.df <- rbind(deathtime.hiv.df, deathtime.norm.df)

  id <- vector()
  infectime <- vector()
  removtime <- vector()
  for(i in 1:length(id.infec)){
    for(j in 1:nrow(deathtime.df)){
      if(id.infec[i] == deathtime.df$p1ID[j]){
        id <- c(id, id.infec[i])
        infectime <- c(infectime, infec.id.time$InfectTime[i])
        removtime <- c(removtime, deathtime.df$eventtime[j])
      }
    }
  }

  id.time.infec.die.fin.raw <- as.data.frame(cbind(id,infectime,removtime))

  names(id.time.infec.die.fin.raw) <- c("ID", "InfectTime", "fin.time")

  id.time.infec.die.fin <- id.time.infec.die.fin.raw

  DF.id.infec.remov.time <- rbind(as.data.frame(id.time.infec.alive.fin), as.data.frame(id.time.infec.die.fin))

  D <- sum((DF.id.infec.remov.time$fin.time - DF.id.infec.remov.time$InfectTime), na.rm = TRUE)/nrow(DF.id.infec.remov.time)


  Ro = beta * C * D

  return(Ro)

}






