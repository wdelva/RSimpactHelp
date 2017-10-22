
intervention.introduced <- function(simulation.type = "pre.hhohho", simulation.start = "1970-03-31"){

  start.time.sim <- as.Date(simulation.start)

  if(simulation.type == "simpact-cyan"){
    # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.

    art.time1 <- as.Date("2004-03-31")
    art.intro <- list()
    art.intro["time"] <- round(as.numeric(difftime(art.time1,start.time.sim, units = "days")/365.242), 0)
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original value @ sim start
    art.intro["monitoring.cd4.threshold"] <- 100
    art.intro["monitoring.interval.piecewise.cd4s"] <- "100,250"
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
    art.time2 <- as.Date("2007-03-31")
    art.intro2 <- list()
    art.intro2["time"] <- round(as.numeric(difftime(art.time2,start.time.sim, units = "days")/365.242),0)
    art.intro2["monitoring.cd4.threshold"] <- 200
    art.intro2["monitoring.interval.piecewise.cd4s"] <- "200,350"
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.time3 <- as.Date("2010-03-31")
    art.intro3 <- list()
    art.intro3["time"] <- round(as.numeric(difftime(art.time3,start.time.sim, units = "days")/365.242),0)
    art.intro3["monitoring.cd4.threshold"] <- 350
    art.intro3["monitoring.interval.piecewise.cd4s"] <- "350,500"
    art.intro3["diagnosis.genderfactor"] <- 1

    art.time4 <- as.Date("2013-03-31")
    art.intro4 <- list()
    art.intro4["time"] <- round(as.numeric(difftime(art.time4,start.time.sim, units = "days")/365.242),0)
    art.intro4["monitoring.cd4.threshold"] <- 500
    art.intro4["monitoring.interval.piecewise.cd4s"] <- "500,650"
    art.intro4["diagnosis.genderfactor"] <- 0.5

    # person.art.accept.threshold.dist.fixed.value

    iv <- list(art.intro, art.intro2, art.intro3, art.intro4)

  }else if(simulation.type == "pre.hhohho"){
    ### Eligibility increased from <200 to <350 in 2010 and then to <500 in 2015

    # Simulation starts in 1970. After 33 years (in 2003), ART is introduced.

    art.time1 <- as.Date("2003-03-31")
    art.intro <- list()
    art.intro$facilities.outfile.facilityxypos <- "" #reset the writing of the facilities.xy position to none.
    art.intro["time"] <- round(as.numeric(difftime(art.time1,start.time.sim, units = "days")/365.242),0)
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original @ sim start
    art.intro["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro["monitoring.cd4.threshold.poststudy"] <- 20000
    art.intro["monitoring.interval.piecewise.cd4s"] <- "100,250"
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2014:500
    art.time2 <- as.Date("2007-03-31")
    art.intro2 <- list()
    art.intro2["time"] <- round(as.numeric(difftime(art.time2,start.time.sim, units = "days")/365.242),0)
    art.intro2["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro2["monitoring.interval.piecewise.cd4s"] <- "200,350"
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.time3 <- as.Date("2010-03-31")
    art.intro3 <- list()
    art.intro3["time"] <- round(as.numeric(difftime(art.time3,start.time.sim, units = "days")/365.242),0)
    art.intro3["monitoring.cd4.threshold.prestudy"] <- 350
    art.intro3["monitoring.interval.piecewise.cd4s"] <- "350,500"
    art.intro3["diagnosis.genderfactor"] <- 1

    art.time4 <- as.Date("2014-03-31")
    art.intro4 <- list()
    art.intro4["time"] <- round(as.numeric(difftime(art.time4,start.time.sim, units = "days")/365.242),0)
    art.intro4["monitoring.cd4.threshold.prestudy"] <- 500
    art.intro4["monitoring.interval.piecewise.cd4s"] <- "500,650"
    art.intro4["diagnosis.genderfactor"] <- 0.5

    art.time7 <- as.Date("2016-10-01")
    art.intro7 <- list()
    art.intro7["time"] <- round(as.numeric(difftime(art.time7,start.time.sim, units = "days")/365.242),0)
    art.intro7["monitoring.cd4.threshold.prestudy"] <- 20000
    art.intro7["monitoring.interval.piecewise.cd4s"] <- "500,2000"
    art.intro7["diagnosis.genderfactor"] <- 0.5

    #reset stepinterval to last facility and end of study
    art.time5 <- as.Date("2016-09-02")
    stepinterval <- list()
    stepinterval["time"] <- round(as.numeric(difftime(art.time5,start.time.sim, units = "days")/365.242),0)
    stepinterval["maxart.stepinterval"] <- 1/12

    #reset stepinterval to the end of study
    art.time6 <- as.Date("2016-10-02")
    stepinterval1 <- list()
    stepinterval1["time"] <- round(as.numeric(difftime(art.time6,start.time.sim, units = "days")/365.242),0)
    stepinterval1["maxart.stepinterval"] <- 7/12


    iv <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro4, stepinterval, stepinterval1 )
  }else{
    iv <-list()
  }

  return(iv)
}
