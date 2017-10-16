#' Build a list of HIV transmission networks data with past demographic of virus populations
#' The variation of viral load is tracked at transmission and sampling times and their rates of changes
#' Each transmission network emerges from one seed individual in the population.
#' Thus we do not consider super-infection in our case.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param endpoint Only transmission events that took place before this point in simulation time, are captured in the output.
#' @return a list with the transmission networks data ( each of them considering different sampling/removal dates), as required by the epi2tree function.
#' infec.times, Transmission time
#' samp.times, Sampling time
#' rec.IDs, recipients IDs
#' don.IDs, donors IDs
#' t.infec, time interval from getting the infection until current transmission for donors
#' t.samp, time interval from getting the infection until this sampling for recipients
#' init.vload, initial viral load for recipient for recipients
#' infec.vload, viral load of the donor at infection time when is going to transmit
#' samp.vload, viral load of the receipient at sampling time for recipients
#' rate.1back.infec.vload, rate of change of the donor from present at infection time to initial viral load
#' rate.2back.infec.vload, rate of change of the donor from present at infection time to previous viral load
#' rate.1back.samp.vload, rate of change of the recipient from present at sampling time to initial viral load
#' rate.2back.samp.vload, rate of change of the recipient from present at sampling time to previous viral load

#' @examples
#' transm.col.ls <- transNetworkPastDemographic(datalist = datalist,endpoint = 40)
#' @note
#' transm.col.ls[[1]] is the transmission network of the first seed
#' @import igraph
#' @import data.table

# Build a transmission network data per seed to be handled by epi2tree function of expotree package



transNetworkPastDemographic <- function(datalist = datalist, endpoint = 40){

  # HIV seed time
  hivseed.time <- datalist$etable[eventname=="HIV seeding"]$eventtime


  # 1. Table of donors and recipients and time of infection

  ## person table of all infected people
  infec.all.raw <- as.data.frame(datalist$ptable[InfectTime !="Inf"])

  infec.all <- infec.all.raw[which(infec.all.raw$InfectTime <= endpoint), ]

  ## recipients
  id.infec.all <- infec.all$ID

  ## donors
  id.infec.origin.all <- infec.all$InfectOrigID

  ## time of infection
  time.infec.all <- infec.all$InfectTime

  ## table of donors, recipients,  and time of infection including the eternal donor -1

  recDonTime <- cbind(id.infec.origin.all, id.infec.all, time.infec.all)
  recDonTime.data <- as.data.frame(recDonTime)

  ## order the table by time of infection (increasing)
  recDonTimeOrd <- recDonTime.data[order(recDonTime.data$time.infec.all),]


  # 2. Among donors and recipients who got infection by seed/transmission event

  #### Consider the infection types

  # -1: non infected people
  # 0: infected due to seed event >> this help to identifie the seeds IDs
  # 1: infected due to a transmission event

  # person table with all infected people (by seed or transmission event),
  # we then remove those who are not infected
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])

  pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= endpoint),]

  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)

  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do

  # person table of infected individuals by transmission event
  pers.table.trans <- subset(pers.infec,pers.infec$InfectType==1)



  # 3. One Network per Seed ID: now automated for all the seeds


  ## Now as we have the table of donors, recipients and time of infection,
  ## as we known the seed IDs, let them trace pathways of each seed
  ## bear in mind that each individual can be infected one time only
  # this means that a recipient is unique in the column of recipients

  ## To get the pathways of each seed, let use igraph function to filter these pathways

  # Take the table of donors, recipients and time of infection
  # which is NewrecDonTimeOrd and make a graph object

  library(igraph)

  ga <- recDonTimeOrd

  ga[,1] <- as.character(recDonTimeOrd[,1]) # donors
  ga[,2] <- as.character(recDonTimeOrd[,2]) # recipients
  gag = as.matrix(ga)
  ga.graph = graph.edgelist(gag[,1:2])


  # Let know extract from the above big network
  # the subcomponent for each seed (transmission pathways due to each seed)

  subcomps <- vector("list", length(seeds.id))
  for (i in 1: length(seeds.id)){
    subcomps[[i]] <- subcomponent(ga.graph, paste(seeds.id[i]), "out")
  }

  # table of donors, recipients,  and time of infection of infected people due to each seed
  net.recdontime.raw <- vector("list", length(seeds.id)) # with eternal infector -1
  init.recdontime <- vector("list", length(seeds.id)) # take the eternal infector -1
  net.recdontime <- vector("list", length(seeds.id)) # remove init.recdontime in net.recdontime.raw
  for(i in 1:length(seeds.id)){
    for(j in 1:length(subcomps)){
      if(seeds.id[i] == as.numeric(names(subcomps[[j]]))[1]){
        net.recdontime.raw[[i]] <- subset(recDonTimeOrd,recDonTimeOrd$id.infec.origin.all%in%as.numeric(names(subcomps[[j]])))
        init.recdontime[[i]] <- c(-1,seeds.id[i], unique(as.numeric(datalist$itable$hivseed.time)))
        net.recdontime[[i]] <- rbind(init.recdontime[[i]],net.recdontime.raw[[i]]) # 1st col donors, 2nd recipients and 3rd time of infection
      }
    }
  }



  ################################################################################
  ### itimes, dtimes, id and parent  for epi2tree function of expoTRee package ###
  ################################################################################

  # Note that for dtimes are times of sampling or removing,
  # they are diagnosis time or death time.
  # Thus, we have
  # diagnosed people nd died >> save death time
  # diagnosed people and alive >> save diagnosed time
  # non diagnosed people and alive >> save endpoint time
  # non diagnosed people and died >> save death time


  # Recode the ids of donors and recipients
  # Recipients help to recorde since we know that they are unique in their column

  # ids (infections ids = recipients ids)

  # table of donors, recipients,  and time of infection of infected people due to see 920
  # augmented with new recorded ids

  id1 <- vector("list", length(seeds.id))
  new.recdontime <- vector("list", length(seeds.id))
  # new.parent <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id1[[i]] <- as.vector(seq(from = 0, to = length(net.recdontime[[i]][,2])-1, by = 1))
    new.recdontime[[i]] <- cbind(id1[[i]],net.recdontime[[i]][,2],net.recdontime[[i]][,1],net.recdontime[[i]][,3])
  }

  # parents ids (donors)
  # the first donor is the universal donor -1
  NewParent <- vector("list", length(seeds.id))
  dat.recdontime <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    NewParent1 <- c(-1)
    for(j in 1:length(new.recdontime[[i]][,3])){
      for(k in 1:length(new.recdontime[[i]][,2])){
        if(new.recdontime[[i]][j,3] == new.recdontime[[i]][k,2]){
          NewParent1 <- c(NewParent1, new.recdontime[[i]][k,1])

        }
        else{
          NewParent1 <- c(NewParent1)
        }
      }
    }
    NewParent[[i]] <- c(NewParent1)

    # per seed we have a contact matrix with
    # c("NewID", "OldID", "OldParent", "NewParent", "InfecTime")
    dat.recdontime[[i]] <- cbind(new.recdontime[[i]][,1],new.recdontime[[i]][,2],new.recdontime[[i]][,3],NewParent[[i]], new.recdontime[[i]][,4])
  }

  # Searching sampling and removal times

  diagnosis.event <- as.data.frame(datalist$etable[eventname == "diagnosis"])
  death.hiv.event <- as.data.frame(datalist$etable[eventname == "aidsmortality"])
  death.norm.event <- as.data.frame(datalist$etable[eventname == "normalmortality"])
  keeps <- c("p1ID","eventtime")

  diagnosis <- diagnosis.event[keeps]
  aids.death <- death.hiv.event[keeps]
  normal.death <- death.norm.event[keeps]
  deaths.all <- as.data.frame(rbind(aids.death,normal.death))
  deaths.all.ord <- deaths.all[order(deaths.all$eventtime),]

  # since people can be diagnosed more than one time remove duplicates!
  diagnosis.unique <-  diagnosis[!duplicated(diagnosis$p1ID), ]
  # or subset(diagnosis, !duplicated(p1ID))

  # people alive at the end of simulation >> use of alive.infected() FUNCTION
  alive.at.endpoint <- alive.infected(datalist = datalist, timepoint = endpoint,
                                      site = "All")

  # person table of all people alive at the end of simulation but are HIV positive
  alive.hiv <- subset(alive.at.endpoint, alive.at.endpoint$InfectType !=-1)

  # person table of people alive at the end of simulation but are HIV positive
  # due to only each seed
  transm.dat <- vector("list", length(seeds.id))
  transm.dat.id <- vector("list", length(seeds.id))
  alive.hiv.pers <- vector("list", length(seeds.id))

  for(i in 1:length(seeds.id)){
    transm.dat[[i]] <- dat.recdontime[[i]]
    transm.dat.id[[i]] <- transm.dat[[i]][,2]
    alive.hiv.pers[[i]] <- subset(alive.hiv, alive.hiv$ID%in%transm.dat.id[[i]]) # person table of people alive at the end of simulation but are HIV positive
  }


  # person table of dead people
  pers <- vector("list", length(seeds.id))
  pers.alive.dat <- vector("list", length(seeds.id))
  pers.alive <- vector("list", length(seeds.id))
  pers.death.dat <- vector("list", length(seeds.id))
  pers.death <- vector("list", length(seeds.id))

  for(i in 1:length(seeds.id)){
    pers.alive.dat[[i]] <- subset(alive.hiv.pers[[i]], alive.hiv.pers[[i]]$ID %in% transm.dat.id[[i]])
    pers.alive[[i]] <- pers.alive.dat[[i]]$ID # alive people in the transmission network

    pers.death.dat[[i]] <- subset(deaths.all.ord, deaths.all.ord$p1ID%in%transm.dat.id[[i]]) # pers replaced by transm.dat.id[[i]]
    pers.death[[i]] <- pers.death.dat[[i]]$p1ID # dead people in the transmission network

  }


  # id of dead people in the transmission network = pers - alive.hiv
  id.dead <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id.dead[[i]] <- subset(transm.dat.id[[i]], !transm.dat.id[[i]]%in%alive.hiv.pers[[i]]$ID)
  }

  # id of alive people

  id.alive <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id.alive[[i]] <- setdiff(transm.dat.id[[i]],id.dead[[i]])
  }


  # split times
  ##############

  # 1. alive and diagnosed people > id & time of diagnosis
  alive.diag <- vector("list", length(seeds.id))
  d1 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    alive.diag[[i]] <- subset(diagnosis.unique, diagnosis.unique$p1ID%in%pers.alive[[i]]) # OK
    d1[[i]] <- as.data.frame(cbind(alive.diag[[i]]$p1ID, alive.diag[[i]]$eventtime))
  }


  # 2. died and diagnosed people > id & time of diagnosis
  died.diag <- vector("list", length(seeds.id))
  d2 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    died.diag[[i]] <- subset(diagnosis.unique, diagnosis.unique$p1ID%in%id.dead[[i]]) # OK
    d2[[i]] <- as.data.frame(cbind(died.diag[[i]]$p1ID, died.diag[[i]]$eventtime))
  }


  # 3. alive non diagnosed people, below person table  > Id & time = endpoint
  id.non.diag.all <- vector("list", length(seeds.id))
  alive.non.diag <- vector("list", length(seeds.id))
  d3 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    # Ids of all no diagnosed
    id.non.diag.all[[i]] <- setdiff(transm.dat.id[[i]],c(alive.diag[[i]]$p1ID, died.diag[[i]]$p1ID))
    alive.non.diag[[i]] <- subset(alive.hiv.pers[[i]], alive.hiv.pers[[i]]$ID%in%id.non.diag.all[[i]])

    d3[[i]] <- as.data.frame(cbind(alive.non.diag[[i]]$ID, rep(endpoint, length(alive.non.diag[[i]]$ID))))
  }


  # 4. died and non diagnosed > Id & death time
  died.non.diag <- vector("list", length(seeds.id))
  d4 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    died.non.diag[[i]] <- subset(deaths.all.ord, deaths.all.ord$p1ID%in%id.non.diag.all[[i]])
    d4[[i]] <- as.data.frame(cbind(died.non.diag[[i]]$p1ID, died.non.diag[[i]]$eventtime))
  }


  # All together: dtimes

  # id and sampling times together

  times.raw <- vector("list", length(seeds.id))
  id.times.raw <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    times.raw[[i]] <- do.call("rbind", list(d1[[i]], d2[[i]], d3[[i]], d4[[i]]))

    id.times.raw[[i]] <- cbind(times.raw[[i]]) # , transm.dat.id[[i]] are infectees for i^th seed
  }


  # rearange the above table according to the order of infection time which is the same
  # order in which we have old ids for each seed
  time.id.dt <- vector("list", length(seeds.id))
  time.dt.all <- vector("list", length(seeds.id))
  time.dt <- vector("list", length(seeds.id))
  time.id <-vector("list", length(seeds.id))


  # define the function to rearange the order anmed times.pers.sort()
  times.pers.sort <- function(id.times.raw,transm.dat.id){
    for(k in 1:length(seeds.id)){
      id.time.ord <- as.data.frame(id.times.raw[[k]]) # old ids and their sampling times
      id.old <- transm.dat.id[[k]] # old ids of infectees for each seed in the infection time order
      time.id.dt <- NULL
      time.dt.init <- vector()
      time.id.init <- vector()
      for(i in 1:nrow(id.time.ord)){
        for(j in 1:length(id.old)){
          if(id.old[i] == id.time.ord$V1[j]){
            # time.id.dt <- cbind(c(time.id, id.times.raw$pers[i]), c(time.dt, id.times.raw$V2[j]))

            time.dt.init <- c(time.dt.init, id.time.ord$V2[j])

            time.id.init <- c(time.id.init, id.old[i])

          }
        }
        time.id.dt <- cbind(time.dt.init, time.id.init)
      }
      time.id[[k]] <- time.id.dt
    }
    return(time.id)
  }

  b = times.pers.sort(id.times.raw, transm.dat.id)

  # to build the epi object handled by epi2tree function
  # to build a transmission tree, we reverse the time for infections time
  # itimes show how old is a given transmitted infection


  # initialize the list of of epi object (one per seed)


  # Viral load and its rates of change and time period of infection

  # From dat.recdontime data frame
  # From b, 1 sampling time and 2 Ids (Rec)
  # 2 Ids (rec), 3 Parents (Don), 5 Infectin time

  vltable <- datalist$vltable
  vldat.init <- vector("list", length(seeds.id)) # viral load at strating of infection for recipient
  vldat.samp <- vector("list", length(seeds.id)) # viral load at sampling for recipients
  vldat.infec <- vector("list", length(seeds.id)) # viral load at infection time for the donor
  t.samp.table <- vector("list", length(seeds.id)) # time interval from infection to sampling for recipient
  t.infec.table <- vector("list", length(seeds.id)) # time interval from infection to transmission for donor
  rate1.rec <- vector("list", length(seeds.id))  # viral load rate between sampling and initial viral for recipients
  rate2.rec  <- vector("list", length(seeds.id)) # viral load rate between sampling and previous change of viral load for recipients
  rate1.don <- vector("list", length(seeds.id)) # viral load rate between transmission and initial viral load for donors
  rate2.don <- vector("list", length(seeds.id)) # viral load rate between transmission and previous change of viral load for donors


  for(i in 1:length(seeds.id)){
    ids <- dat.recdontime[[i]][,2] # recipients
    don.pars <- dat.recdontime[[i]][,3] # donors
    vl.ids <- subset(vltable, vltable$ID%in%ids) # Viral load tale for these recipients individuals
    # vl.pars <- subset(vltable, vltable$ID%in%don.pars) # Viral load tale for donors individuals

    dat <- as.data.table(dat.recdontime[[i]])
    infecTime.all <- dat$V5 # infection time
    sampTime.all  <- b[[i]][,1] # sampling time

    t.samp.table[[i]] <- sampTime.all - infecTime.all

    # Viral load table of everyone
    vl.id <- vector("list",length(ids)) # > REC
    vl.init <- vector() # initial viral load arranged regarind recipients > REC
    vl.init.don <- vector() # initial viral load arranged regarind donors > DON
    vl.samp <- vector() # viral load at sampling time > REC
    vl.infec <- vector() # viral load at infection time > DON
    t.infection <- vector() # time interval from infection until current transmission > DON
    individuals.rec <- vector() # ID of recipients
    individuals.don <- vector() # ID of donors
    r1.rec <- vector() # viral load rate between sampling and initial viral for recipients
    r2.rec <- vector() # viral load rate between sampling and previous change of viral load for recipients
    r1.don <- vector()  # viral load rate between transmission and initial viral for donor
    r2.don <- vector() # viral load rate between transmission and previous change of viral load for donors
    for(j in 1:length(ids)){ # retrieve viral load table for each individual
      id <- ids[j] # rec. individual
      id.don.par <- dat$V3[j] # don. individual
      vl.id[[j]] <- vl.ids[ID==id] # rec. viral load
      vl.dat <- vl.id[[j]]
      infecTime <- dat$V5[j] # infection time
      sampTime  <- b[[i]][j] # sampling time

      vl.init <- c(vl.init,vl.dat$Log10VL[1])

      # vl at samp
      viralLoadSamp <- function(vl.dat){
        if(sampTime >= vl.dat$Time[nrow(vl.dat)]){
          vl.samp.fun <- vl.dat$Log10VL[nrow(vl.dat)]
        }else{
          vlTime <- vl.dat$Time
          index <- which.min(abs(vlTime-sampTime)) # check  nearest time value

          vl.samp.fun <- vl.dat$Log10VL[index] # get viral load at that time interval
        }
        return(vl.samp.fun)
      }
      vl.samp.index <- viralLoadSamp(vl.dat)

      vl.samp <- c(vl.samp,vl.samp.index)

      # vl at infection time for donor
      viralLoadInfec <- function(vl.dat){
        if(infecTime >= vl.dat$Time[nrow(vl.dat)]){
          vl.infec.fun <- vl.dat$Log10VL[nrow(vl.dat)]
        }else{
          vlTime.i <- vl.dat$Time
          index.i <- which.min(abs(vlTime.i-infecTime)) # check  nearest time value

          vl.infec.fun <- vl.dat$Log10VL[index.i]
        }
        return(vl.infec.fun)
      }
      vl.infec.index <- viralLoadInfec(vl.dat)
      vl.infec <- c(vl.infec, vl.infec.index)

      InPutDat <- dat
      # Find parent of the current donor here "h": due to above mentionned argument
      parent.find <- function(h){
        for(l in 1:length(InPutDat$V2)){
          if(InPutDat$V2[l] == h){ # when current donor was a recipient
            par <- InPutDat$V3[l]
          }
        }
        return(par)
      }

      t.infectionFun <- function(InPutDat,j){
        parentDon <- InPutDat$V3[j]
        t.infection.give <- InPutDat$V5[j]
        if(parentDon==-1){
          t.infection.out <- 10 # for the universal infector -1 we set up 1000 yrs
        }else{
          par.infector <- parent.find(parentDon) # parent of current donor
          par.infector.dat <- InPutDat[V2==parentDon] # table of transmission of this donor
          # in the past when he (she) got the infection
          t.infection.got <- par.infector.dat$V5 # the time when s/he got it
          t.infection.out <-  t.infection.give - t.infection.got # time of todays transmission - time when he got the infection
        }
        return(t.infection.out)
      }
      t.infec.pt <- t.infectionFun(InPutDat,j)
      t.infection <- c(t.infection,t.infec.pt)

      # Rates of change of viral load #
      #################################

      ## REC

      # r1 > further, rate of change of viral load between current viral load
      # at sampling time and initial viral load when got the infection
      r1.rec.index <- vl.samp.index/vl.dat$Log10VL[1]
      r1.rec <- c(r1.rec, r1.rec.index)

      # r2 > rate of change in time interval: sampling time and time of previous change of viral load

      # vl at previous viral load change before sampling
      viralLoadSamp.rate <- function(vl.dat){
        if(sampTime >= vl.dat$Time[nrow(vl.dat)]){
          vl.samp.fun.rate <- vl.dat$Log10VL[nrow(vl.dat)]
        }else{
          vlTime <- vl.dat$Time
          index <- which.min(abs(vlTime-sampTime)) # check  nearest time value
          if(index==1){
            vl.samp.fun.rate <- vl.dat$Log10VL[1] # if no changes index==1 and we remain with same amount of viral load in that time interval
          }else{
            vl.samp.fun.rate <- vl.dat$Log10VL[index-1]
          }
        }
        return(vl.samp.fun.rate)
      }
      vl.samp.index.rate <- viralLoadSamp.rate(vl.dat)

      r2.rec.index <- vl.samp.index/vl.samp.index.rate
      r2.rec <- c(r2.rec, r2.rec.index)


      ## DON

      # initial viral load for donors > we do have initial vl for all recipients
      # rearrange them regarding donors recipients

      vl.init.don.FUN <- function(dat, vl.ids){
        id.don.par <- dat$V3 # donor
        for (l in 1:length(id.don.par)){
          parent.don <- id.don.par[l]
          if(parent.don==-1){
            vl.init.index <- 9
          }else{
            vl.donor.index.dat <- vl.ids[ID==parent.don]
            vl.donor.index <- vl.donor.index.dat
            vl.init.index <- vl.donor.index$Log10VL[1]
          }
        }
        return(vl.init.index)
      }
      vl.init.don.index <- vl.init.don.FUN(dat, vl.ids)
      vl.init.don <- c(vl.init.don, vl.init.don.index)

      r1.don.index <- vl.infec.index/vl.init.don.index

      r1.don <- c(r1.don, r1.don.index)


      # vl at infection time for donor
      viralLoadInfec <- function(vl.dat){
        if(infecTime >= vl.dat$Time[nrow(vl.dat)]){
          vl.infec.fun <- vl.dat$Log10VL[nrow(vl.dat)]
        }else{
          vlTime.i <- vl.dat$Time
          index.i <- which.min(abs(vlTime.i-infecTime)) # check  nearest time value

          vl.infec.fun <- vl.dat$Log10VL[index.i]
        }
        return(vl.infec.fun)
      }
      vl.infec.index <- viralLoadInfec(vl.dat)

      # vl at previous viral load change before transmission
      viralLoadInfec.rate <- function(vl.dat){
        if(infecTime >= vl.dat$Time[nrow(vl.dat)]){
          vl.infec.fun.rate <- vl.dat$Log10VL[nrow(vl.dat)]
        }else{
          vlTime <- vl.dat$Time
          index <- which.min(abs(vlTime-infecTime)) # check  nearest time value
          if(index==1){
            vl.infec.fun.rate <- vl.dat$Log10VL[1] # if no changes index==1 and we remain with same amount of viral load in that time interval
          }else{
            vl.infec.fun.rate <- vl.dat$Log10VL[index-1]
          }
        }
        return(vl.samp.fun.rate)
      }
      vl.infec.index.rate <- viralLoadInfec.rate(vl.dat)

      r2.don.index <- vl.infec.index/vl.infec.index.rate
      r2.don <- c(r2.don, r2.don.index)

    }

    vldat.init[[i]] <- vl.init # REC
    vldat.samp[[i]] <- vl.samp # REC
    vldat.infec[[i]] <- vl.infec # DON
    t.infec.table[[i]] <- t.infection # DON
    rate1.rec[[i]] <- r1.rec # REC
    rate2.rec[[i]] <- r2.rec # REC
    rate1.don[[i]] <- r1.don # DON
    rate2.don[[i]] <- r2.don # DON


  }


  transm.Colescent.ls <- vector("list", length(seeds.id))

  transNet <- list()
  for(q in 1:length(seeds.id)){
    transNet$infec.times <- (dat.recdontime[[q]][,5]) # Transmission time
    transNet$samp.times <- (b[[q]][,1]) # Sampling time
    transNet$rec.IDs <- dat.recdontime[[q]][,2] # recipients IDs
    transNet$don.IDs <- dat.recdontime[[q]][,3] # donors IDs
    transNet$t.infec <- t.infec.table[[q]] # time interval from getting the infection until current transmission > DON.
    transNet$t.samp <- t.samp.table[[q]] # time interval from getting the infection until this sampling > REC.
    transNet$init.vload <- vldat.init[[q]] # initial viral load for recipient > REC.
    transNet$infec.vload <- vldat.infec[[q]] # viral load of the donor at infection time > DON.
    transNet$samp.vload <- vldat.samp[[q]] # viral load of the receipient at sampling time > REC.
    transNet$rate.1back.infec.vload <- rate1.don[[q]] # rate of change of the donor from present to previous viral load > DON.
    transNet$rate.2back.infec.vload <- rate2.don[[q]] # rate of change of the donor from present to before previous viral load > DON.
    transNet$rate.1back.samp.vload <- rate1.rec[[q]] # rate of change of the recipient from present to previous viral load > REC.
    transNet$rate.2back.samp.vload <- rate2.rec[[q]] # rate of change of the recipient from present to before previous viral load > REC.

    transm.Colescent.ls[[q]] <- transNet
  }

  return(transm.Colescent.ls)

}
