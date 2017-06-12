# Compute target statistics

# growth rate: OK
# relationship rate
# transmission rate: OK

##################################
# 1. Growth rate: done by Trust ##
##################################

growth.rate <- pop.growth.calculator(datalist = chunk.datalist.test,
                                     timewindow = c(0, timewindow.max = end.time.wind))

###########################
# 2. Transmission rate ####
###########################

# Transmission rate in time interval
transmission.rate <- function(datalist = datalist,
                                  timewindow = c(0, 40)){

  Infec.pop.table <- datalist$ptable[InfectType==1]

  numb.infec.pop <- nrow(Infec.pop.table %>%
    subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))


  transm.rate <- numb.infec.pop / diff(timewindow)

  return(transm.rate)
}

# Transmission rate in time interval subdivided in different steps

# Transmission rate
transmission.rate.int <- function(datalist = datalist,
                              timewindow = c(0, 40), by=1){

  Infec.pop.table <- datalist$ptable[InfectType==1]

  upper.limit <- ceiling(diff(timewindow)/by)

  interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

  infec.pop.int <- vector()
  trans.rate.int <- vector()


  for(i in 0:(upper.limit-2)){

    timewindow.int <- c(interv.time[1+i], interv.time[2+i])

    infec.pop.numb <- nrow(Infec.pop.table %>%
      subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))


    infec.pop.int <- c(infec.pop.int,  infec.pop.numb)


    trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))

  }

  return(trans.rate.int)
}


# Transmission ratio for men and women

transm.gender.ratio <- function(datalist = datalist,
                              timewindow = c(0, 40)){

  pers.table.hiv <- datalist$ptable[InfectType==1]

  pers.table.hiv.window <- pers.table.hiv %>%
    subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1])

  numb.transm <- nrow(pers.table.hiv.window)

  numb.transm.men <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID))) # Gender 0 men

  numb.transm.women <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID)))# Gender 1 women

  transm.ratio <- vector("list", length(c(0,1)))

  for(i in 1:2){
    transm.ratio$men <- (numb.transm.men/numb.transm)  #/diff(timewindow)
    transm.ratio$women <- (numb.transm.women/numb.transm)  #/diff(timewindow)
  }
  return(transm.ratio)
}

###############################
# 3. Relationship rate ########
###############################

# Women, men and overall relationships rate in time interval
# Overall relationship rate
# Rste of women (men) who have been invlved in relationship udirng a time interval.s

rels.rate <- function(datalist = datalist,
                              timewindow = c(0, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.rate <- vector("list", length(c(0,1,2)))

  for(i in 1:3){
    rels.rate$overall <- (numb.rels)/ diff(timewindow)
    rels.rate$men <- (numb.rels.men)/diff(timewindow)
    rels.rate$women <- (numb.rels.women)/diff(timewindow)
  }
  return(rels.rate)
}


# Relationship rate in time interval subdivided in different steps

# Relationship rate
rels.rate.int <- function(datalist = datalist,
                                  timewindow = c(0, 40), by=1){

  Rels.table <- datalist$rtable

  upper.limit <- ceiling(diff(timewindow)/by)

  interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

  rels.int <- vector()
  rels.rate.int <- vector()


  for(i in 0:(upper.limit-2)){

    timewindow.int <- c(interv.time[1+i], interv.time[2+i])

    rels.numb <- nrow(Rels.table %>%
                             subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))


    rels.int <- c(rels.int,  rels.numb)


    rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))

  }

  return(rels.rate.int)
}


# Relationship ratio for men and women

rels.gender.ratio <- function(datalist = datalist,
                            timewindow = c(0, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.ratio <- vector("list", length(c(0,1)))

  for(i in 1:3){
    rels.ratio$men <- (numb.rels.men/numb.rels)  #/diff(timewindow)
    rels.ratio$women <- (numb.rels.women/numb.rels)  #/diff(timewindow)
  }
  return(rels.ratio)
}



