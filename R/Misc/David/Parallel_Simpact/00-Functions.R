pacman::p_load(igraph,RSimpactHelper,magrittr,dplyr,nlme,lme4,fitdistrplus,
               lmtest,tidyr,data.table,gsubfn,expoTree,ape,readr,phylosim,utils,
               phangorn,apTreeshape,pcaPP)

freq <- c(0.3353293, 0.2035928, 0.2628077, (1 - sum(c(0.3353293, 0.2035928, 0.2628077))))
rate <- list("a"=0.2, "b"=0.6, "c"=0.12,"d"=0.001, "e"=0.25, "f"=0.24)

overall.rate <- 0.00003 * (8760/52)  # 3 x 10^−5 per nucleotide base per replication cycle (which takes 52 hours). There are 8760 hours in 1 year. ~ 0.005 per nucleotide base per year.

agemixing.lme.errFunction <- function(e)
{
  return(list())
}

simpact.errFunction <-function(e){
  if (length(grep("NaN",e$message)) != 0){
    return(list())
  }
}

ampmodel <- function(data = dplyr::filter(agemix.model[[1]], Gender =="male")) {
  lmer(pagerelform ~ agerelform0 + (1 | ID),
       data = data,
       REML = TRUE)
}

bvar <- function(model) {

  # Outputs a df with between-subject variance, upr & lwr limits

  # Must take an merMod  object

  bsd <- as.numeric(as.data.frame(VarCorr(model))[1,5])


}

pattern.modeller <- function(dataframe,
                             agegroup,
                             timepoint,
                             timewindow,
                             start = FALSE) {
  #Warnings
  if (!is.data.frame(dataframe)) {
    stop("dataframe wrong type")
  }

  if (length(agegroup) != 2) {
    stop("Need an upper and lower bound")
  }

  if (timepoint < 1) {
    stop("Time must be at least 1 year")
  }

  if (timewindow < 0) {
    stop("Window must be a whole number")
  }

  #Constants
  time <- timepoint
  window <- timepoint - timewindow
  lwrage <- agegroup[1]
  uprage <- agegroup[2]


  if (start == TRUE) {
    # This only includes relationships that started
    # in the time window

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
                      (FormTime <= time & FormTime >= window) &
                      age >= lwrage &
                      age < uprage &
                      Gender == "male" &
                      TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
                      (FormTime <= time & FormTime >= window) &
                      age >= lwrage &
                      age < uprage &
                      Gender == "female" &
                      TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

  } else {
    # This includes all relationships that were ongoing
    # at somepoint during the time window, but may have
    # started long before the time window.

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
                      DisTime > window &
                      age >= lwrage &
                      age < uprage &
                      Gender == "male" &
                      TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
                      DisTime > window &
                      age >= lwrage &
                      age < uprage &
                      Gender == "female" &
                      TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

  }

  comb <- bind_rows(men, women)

  modoutput <- matrix(nrow = 1, ncol = 14) %>%
    as.data.frame()
  colnames(modoutput) <- c("slopem", "slopew", "interceptm", "interceptw",
                           "powerm",  "lowerpowerm", "upperpowerm", "powerw",
                           "lowerpowerw", "upperpowerw", "bvarm", "bvarw",
                           "wvarm", "wvarw")

  agemix.pieces <- list(comb, modoutput)

  return(agemix.pieces)

}

degree.df.maker <- function(dataframe.df, agegroup = c(15, 30), hivstatus = 0,
                             survey.time = 40, window.width = 1, gender.degree = "female",
                             only.new = TRUE)
{
  dataframe.rels.df <- dataframe.df
  dataframe.rels.df <- dplyr::select(dataframe.rels.df, -episodeorder, -FormTime, -DisTime)
  dataframe.rels.df <- unique.data.frame(dataframe.rels.df)
  rels.form.dis.df <- dplyr::summarise(dplyr::group_by(dataframe.df,
                                                       relid), FormTime = min(FormTime), DisTime = max(DisTime))
  dfnew <- dplyr::left_join(x = dataframe.rels.df, y = rels.form.dis.df,
                            by = "relid")
  dfnew <- dplyr::filter(dfnew, TOD > survey.time)
  {
    if (hivstatus == 0) {
      dfnew <- dplyr::filter(dfnew, InfectTime > survey.time)
    }
    else if (hivstatus == 1) {
      dfnew <- dplyr::filter(dfnew, InfectTime <= survey.time)
    }
  }
  {
    if (only.new) {
      dfnew <- dplyr::filter(dfnew, FormTime >= survey.time -
                               window.width, FormTime < survey.time, DisTime >
                               survey.time - window.width, Gender == gender.degree,
                             survey.time - TOB >= agegroup[1], survey.time -
                               TOB < agegroup[2])
    }
    else {
      dfnew <- dplyr::filter(dfnew, FormTime < survey.time,
                             DisTime > survey.time - window.width, Gender ==
                               gender.degree, survey.time - TOB >= agegroup[1],
                             survey.time - TOB < agegroup[2])
    }
  }
  uniqueOut <- dfnew %>% dplyr::select(ID, relid) %>% distinct %>%
    rename(Degree = relid)
  if (dim(uniqueOut)[1] != 0) {
    degreedata.df <- aggregate(Degree ~ ID, data = uniqueOut,
                               length)
  } else {
    degreedata.df <- data.frame(ID = NA, Degree = NA)
  }
  return(degreedata.df)
}


concurr.pointprev.calculator <- function(datalist,
                                         timepoint = datalist$itable$population.simtime[1] - 0.5){

  output <- data.table()
  DTalive.infected <- alive.infected(datalist = datalist,
                                     timepoint = timepoint, site = "All") # First we only take the data of people who were alive at time_i
  agemix.df <- agemix.df.maker(datalist)
  degrees.df <- degree.df.maker(dataframe.df = agemix.df,
                                agegroup = c(15, 50),
                                hivstatus = 2,
                                survey.time = timepoint,
                                window.width = 0,
                                gender.degree = "male",
                                only.new = FALSE)

  number.people.with.cps <- sum(degrees.df$Degree > 1)
  popsize <- nrow(DTalive.infected)
  concurr.pointprevalence <- number.people.with.cps / popsize

  return(concurr.pointprevalence)
}




cd4.atARTinit <- function(datalist = datalist,
                          agegroup = c(15, 30),
                          timewindow = c(15, 30),
                          cd4count=350, site="All"){

  cd4.atARTinit <- age.group.time.window(datalist = datalist,
                                         agegroup = agegroup,
                                         timewindow = timewindow, site="All")

  cd4.atARTinit <- subset(cd4.atARTinit, TreatTime !=Inf) #HIV positive individuals

  raw.df <- data.frame(cd4.atARTinit)
  art.df <- subset(datalist$ttable, ID %in% cd4.atARTinit$ID &
                     TStart > timewindow[1] & TStart < timewindow[2])

  ##What if the person dropped out and come back again?
  art.df <- data.frame(dplyr::summarise(dplyr::group_by(art.df, ID, Gender),
                                        CD4atARTstart = min(CD4atARTstart)))

  #indicate those who started their treatment when their CD4 count was below a given threshold
  art.df <- art.df %>% dplyr::mutate(ART.start.CD4 = CD4atARTstart < cd4count)

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #provide a summary of those that are on treatment and those that started below a threshold
  cd4count.atARTInit <- data.frame(dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                                    TotalCases = n(),
                                                    LessCD4initThreshold =sum(ART.start.CD4)))

  return(cd4count.atARTInit)
}




transmNetworkBuilder.baseline <- function(datalist = datalist, endpoint = 40){

  # HIV seed time
  hivseed.time <- datalist$etable[eventname=="HIV seeding"]$eventtime

  # Simulation time
  population.simtime <- unique(datalist$itable$population.simtime)

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

 # library(igraph)

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
  net.recdontime.raw <- vector("list", length(seeds.id))
  init.recdontime <- vector("list", length(seeds.id))
  net.recdontime <- vector("list", length(seeds.id))
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

  # Note that for dtimes since we are only interested on building a transmission tree dtype = tranmission event
  # Thus dtimes will be fixed at 0


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




  # to build the epi object handled by epi2tree function
  # to build a transmission tree, we reverse the time for infections time
  # itimes show how old is a given transmitted infection


  # initialize the list of of epi object (one per seed)

  transm.ls <- vector("list", length(seeds.id))

  transNet <- list()
  for(i in 1:length(seeds.id)){
    transNet$itimes <- (endpoint - (dat.recdontime[[i]][,5]))
    transNet$dtimes <- rep(0,length(dat.recdontime[[i]][,5]))
    transNet$id <- dat.recdontime[[i]][,1]
    transNet$parent <- dat.recdontime[[i]][,4]

    transm.ls[[i]] <- transNet
  }

  ## Raw data of the transmission matrix, itimes, dtimes, id & parent

  # transNet <- list()
  # for(i in 1:length(seeds.id)){
  #   transNet$itimes <- (dat.recdontime[[i]][,5])-hivseed.time
  #   transNet$dtimes <- rep(0,length(dat.recdontime[[i]][,5]))
  #   transNet$id <- dat.recdontime[[i]][,1]
  #   transNet$parent <- dat.recdontime[[i]][,4]
  #
  #   transm.ls[[i]] <- transNet
  # }

  return(transm.ls)

}


transmNetworkBuilder.diff <- function(datalist = datalist, endpoint = 40){

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

  transm.ls <- vector("list", length(seeds.id))

  transNet <- list()
  for(i in 1:length(seeds.id)){
    transNet$itimes <- (endpoint - (dat.recdontime[[i]][,5]))#+hivseed.time
    transNet$dtimes <- (endpoint - (b[[i]][,1]))
    transNet$id <- dat.recdontime[[i]][,1]
    transNet$parent <- dat.recdontime[[i]][,4]

    transm.ls[[i]] <- transNet
  }

  return(transm.ls)

}

time.mrca.matrix <- function(tree = phylo.tree){

  Mytree <- tree


  # Use of mrca() and branching.time() functions from ape package

  # Symmetric matrix with tips which where entries are the internal nodes which represents the MRCA between two tips
  # from the phylogenetic tree

  Mytree.mrca.tips <- as.data.frame(mrca(Mytree, full = FALSE)) # MRCA for tips only >> full = FALSE

  # The assignment is to use the internal nodes and replace them by branching time


  # tips which represent individuals (sequences)
  inds <- Mytree$tip.label

  # branching time
  branch.time <- branching.times(Mytree) # each internal_node is associated with a branching time

  # Preparing the matrix of time to the MRCA
  ##########################################
  # Change the Mytree.mrca.tips matrix a bit to be easily handle for as igraph or network obejct
  # we have to make the diagonal 0
  n = as.matrix(Mytree.mrca.tips)
  mat.fun.tree <- function(n){
    diag(n) <- 0
    m <- n
    return(m)
  }

  m <- as.data.frame(mat.fun.tree(n))

  # Get internal nodes of the phyloegentic tree

  i.nodes <- NULL
  int.nodes <- function(m){
    for(i in 1:nrow(m)){
      for(j in 1:nrow(m)){
        i.nodes <- c(i.nodes,m[i,j])
      }
    }
    return(unique(i.nodes))
  }

  h = int.nodes(m)
  # remove the first element which is the 0 on diagonal
  g <- h[(2:length(h))] # these are really the internal nodes which represents the MRCA in the phylogenetic tree

  # MRCA times
  # each internal node has its branching time
  # thus, we have to replace the internal nodes in the matrix with their branching times respectively

  mrca.time.val <- as.data.frame(matrix(,nrow(m),nrow(m))) # initiate an n*n empty matrix (n tips)

  mrca.time.matrix <- function(g,m,branch.time){ # write a function which take as parms the internal nodes (g)
    # the matrix of internal nodes (MRCA) which is m and
    # the branching time which correspond with time where the two seq eveoled separelty
    k <- as.data.frame(branch.time)
    p <- k$branch.time

    for(l in 1:length(g)){
      for(i in 1:nrow(m)){
        for(j in 1:nrow(m)){
          if(m[i,j] == g[l]){
            mrca.time.val[i,j] <- p[l]
          }
        }
      }

    }
    return(mrca.time.val)
  }

  v = mrca.time.matrix(g,m,branch.time) # the diagonal is empty with entries NA

  # remove the NA in the diagonal and put 0 since no MRCA of a given sequence comparing at itself
  matrix.mrca.time.func <- function(v){
    diag(v) <- 0
    times.mrca.fine <- v
    return(times.mrca.fine)
  }

  names.inds <- names(m) # get the names of the tips
  mrca.times.final <- matrix.mrca.time.func(v)
  names(mrca.times.final) <- names.inds # put the names of the tips on the columns

  mrca.times.done <- mrca.times.final

  return(mrca.times.done)
}

sequence.simulation <- function(transtree = tree0, seedSeq = hivSeq,
                                base.freq = freq){

  # define the substitution processes >> package phylosim
  # proc <- GTR(rate.params = rate.list,
  #             base.freqs = base.freq)

  proc <- F81(base.freqs = base.freq)

  # attach process to the nucleotides sequence >> package phylosim
  nucleproc <- NucleotideSequence(string = seedSeq, processes = list(list(proc)))

  # plusGamma(nucleproc,proc,alpha)

  # simulate the sequences >> package phylosim
  simsequence <- Simulate(PhyloSim(root.seq = nucleproc,
                                   phylo = transtree))

  # return the alignement
  return(simsequence)
}

transmNetworkBuilder.baseline <- function(datalist = datalist, endpoint = 40){

  # HIV seed time
  hivseed.time <- datalist$etable[eventname=="HIV seeding"]$eventtime

  # Simulation time
  population.simtime <- unique(datalist$itable$population.simtime)

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
  net.recdontime.raw <- vector("list", length(seeds.id))
  init.recdontime <- vector("list", length(seeds.id))
  net.recdontime <- vector("list", length(seeds.id))
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

  # Note that for dtimes since we are only interested on building a transmission tree dtype = tranmission event
  # Thus dtimes will be fixed at 0


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




  # to build the epi object handled by epi2tree function
  # to build a transmission tree, we reverse the time for infections time
  # itimes show how old is a given transmitted infection


  # initialize the list of of epi object (one per seed)

  transm.ls <- vector("list", length(seeds.id))

  transNet <- list()
  for(i in 1:length(seeds.id)){
    transNet$itimes <- (endpoint - (dat.recdontime[[i]][,5]))#-hivseed.time
    transNet$dtimes <- rep(0,length(dat.recdontime[[i]][,5]))
    transNet$id <- dat.recdontime[[i]][,1]
    transNet$parent <- dat.recdontime[[i]][,4]

    transm.ls[[i]] <- transNet
  }

  ## Raw data of the transmission matrix, itimes, dtimes, id & parent

  # transNet <- list()
  # for(i in 1:length(seeds.id)){
  #   transNet$itimes <- (dat.recdontime[[i]][,5])-hivseed.time
  #   transNet$dtimes <- rep(0,length(dat.recdontime[[i]][,5]))
  #   transNet$id <- dat.recdontime[[i]][,1]
  #   transNet$parent <- dat.recdontime[[i]][,4]
  #
  #   transm.ls[[i]] <- transNet
  # }

  return(transm.ls)

}

epi2tree2 <- function(epi){
  make.internal <- function(epi, node, parent, cur.id) {
    node.i <- which(epi$id == node)
    dt <- epi$dtimes[node.i]
    offspring <- epi$parent == node
    if (sum(offspring) > 0) {
      off.id <- epi$id[offspring]
      inf.t <- epi$itimes[offspring]
      dt <- epi$dtimes[node.i]

      intervals <- c(epi$itimes[node.i], inf.t) - c(inf.t, dt)

      # intervals <- epi$itimes[node.i] - dt


      ages <- c(inf.t, dt)
      new.id <- 1:length(inf.t) + cur.id
      cur.id <- max(new.id)
      parents <- c(parent, new.id)
      edge.list <- cbind(parents, c(new.id, node), intervals,
                         ages)
      for (j in 1:length(off.id)) {
        out <- make.internal(epi = epi, node = off.id[j],
                             parent = new.id[j], cur.id = cur.id)
        edge.list <- rbind(edge.list, out[[1]])
        cur.id <- out[[2]]
      }
      return(list(edge.list, cur.id))
    }
    else {
      edge.list <- c(parent, node, epi$itimes[node.i] -
                       dt, dt)
      return(list(edge.list, cur.id))
    }
  }
  edge.list <- make.internal(epi, 0, -1, max(epi$id))
  edges <- matrix(edge.list[[1]][-1, 1:2] + 1, ncol = 2)
  root.edge <- edge.list[[1]][1, 3]
  edge.lengths <- edge.list[[1]][-1, 3]
  tip.label <- as.character(epi$id)
  coords <- expoTree:::epi.coords(epi)[[2]]
  coords[, 1] <- coords[, 1] + 1
  tree <- list(edge = edges, edge.length = edge.lengths, tip.label = tip.label,
               Nnode = max(edge.list[[1]][,1:2]) - max(epi$id), root.edge = root.edge,
               coords = coords)
  class(tree) <- "phylo"
  return(tree)
}

trans.network2tree <- function(transnetwork = transnetwork){

  transm.phylo <- epi2tree2(transnetwork)

  return(transm.phylo)

}



attributes.trans.network <- function(datalist = datalist, endpoint = 40){


  # Simulation time
  population.simtime <- unique(datalist$itable$population.simtime)

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

  #library(igraph)

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

  # Searching sampling and removal times: diagnosis and death

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

    d3[[i]] <- as.data.frame(cbind(alive.non.diag[[i]]$ID, rep(population.simtime - endpoint, length(alive.non.diag[[i]]$ID))))
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

  # Rearange IDs and Sampling times

  # define the function to rearange the order named times.pers.sort()
  # times.pers.sort <- function(id.times.raw,transm.dat.id){
  #   for(k in 1:length(seeds.id)){
  #     id.time.ord <- as.data.frame(id.times.raw[[k]]) # old ids and their sampling times
  #     id.old <- transm.dat.id[[k]] # old ids of infectees for each seed in the infection time order
  #     time.id.dt <- NULL
  #     time.dt.init <- vector()
  #     time.id.init <- vector()
  #     for(i in 1:nrow(id.time.ord)){
  #       for(j in 1:length(id.old)){
  #         if(id.old[i] == id.time.ord$V1[j]){
  #           # time.id.dt <- cbind(c(time.id, id.times.raw$pers[i]), c(time.dt, id.times.raw$V2[j]))
  #
  #           time.dt.init <- c(time.dt.init, id.time.ord$V2[j])
  #
  #           time.id.init <- c(time.id.init, id.old[i])
  #
  #         }
  #       }
  #       time.id.dt <- cbind(time.dt.init, time.id.init)
  #     }
  #     time.id[[k]] <- time.id.dt
  #   }
  #   return(time.id)
  # }
  #
  # b.dt = times.pers.sort(id.times.raw, transm.dat.id)

  # to build the epi object handled by epi2tree function
  # to build a transmission tree, we reverse the time for infections time
  # itimes show how old is a given transmitted infection




  ## Add other filter attributes for individuals:
  ## Gender, number of relatioship before sampling/removal within a time window, and
  ## age at sampling time [ add soon viral load at sampling and CD4]

  # (i) Gender: OK

  # gender.filter <- function(datalist){
  #   gender.all <- vector("list", length(seeds.id))
  #   pers.table.filter <- vector("list", length(seeds.id))
  #   id <- vector("list", length(seeds.id))
  #   gender <- vector("list", length(seeds.id))
  #   pers.table <- as.data.frame(datalist$ptable)
  #
  #   for(k in 1:length(seeds.id)){
  #     pers.table.filter <- subset(pers.table, pers.table$ID%in%b.dt[[k]][,2])
  #     id.raw <- vector()
  #     gender.raw <- vector()
  #     c <- NULL
  #     for(i in 1:nrow(pers.table.filter)){
  #       for(j in 1:nrow(pers.table.filter)){
  #         if(b.dt[[k]][,2][i] == pers.table.filter$ID[j]){
  #           id.raw <- c(id.raw,b.dt[[k]][,2][i])
  #           # id[[k]] <- c(id, b[[k]][,2][i])
  #
  #           gender.raw <- c(gender.raw,pers.table.filter$Gender[j])
  #           #  gender[[k]] <- c(gender,pers.table.filter$Gender[j])
  #         }
  #       }
  #     }
  #
  #     id[[k]] <-id.raw
  #     gender[[k]] <- gender.raw
  #     c <- cbind(id[[k]],gender[[k]])
  #     gender.all[[k]] <- as.data.frame(c)
  #   }
  #   return(gender.all)
  #
  # }
  #
  # filter.gender <- gender.filter(datalist = datalist)
  #
  #
  # # (ii) Number of relatioships at the sampling times
  #
  # count.rels.filter <- function(b.dt=b.dt, datalist=datalist, window=endpoint){
  #   count.rels <- vector("list", length(seeds.id))
  #   rels.ind.times <- vector("list", length(seeds.id))
  #   z1 <- vector("list", length(seeds.id))
  #   z2 <- vector("list", length(seeds.id))
  #
  #   for(k in 1:length(seeds.id)){
  #
  #     # k = 16
  #
  #     rels.table <- datalist$rtable
  #     rels.infec1 <- subset(rels.table, rels.table$ID1%in%b.dt[[k]][,2]) # some IDs in trans network are rec or don
  #     rels.infec2 <- subset(rels.table, rels.table$ID2%in%b.dt[[k]][,2]) # thus, some are in rels.table$ID1 others in rels.table$ID2
  #
  #     count.relationships.1 <- function(k,rels.infec1){ # recipients rels
  #       ids.rels <- b.dt[[k]][,2]
  #       rel.est.numb.time1 <- NULL
  #       rel.est.numb1 <- vector()
  #       rel.est.time1 <- vector()
  #
  #       for(i in 1:nrow(b.dt[[k]])){
  #         if(nrow(rels.infec1) > 0){
  #           for(j in 1:nrow(rels.infec1)){
  #             if(ids.rels[[i]] == rels.infec1$ID1[j]){ # make diff1
  #               rel.est.numb1 <- c(rel.est.numb1,rels.infec1$ID1[j])
  #               rel.est.time1 <- c(rel.est.time1,rels.infec1$FormTime[j])
  #             }
  #           }
  #         }
  #         rel.est.numb.time1 <- cbind(rel.est.numb1,rel.est.time1)
  #       }
  #       return(rel.est.numb.time1)
  #     }
  #     z1[[k]] <- count.relationships.1(k,rels.infec1)
  #
  #     count.relationships.2 <- function(k,rels.infec2){ # recipients rels
  #       ids.rels <- b.dt[[k]][,2]
  #       rel.est.numb.time1 <- NULL
  #       rel.est.numb2 <- vector()
  #       rel.est.time2 <- vector()
  #
  #       for(i in 1:nrow(b.dt[[k]])){
  #         if(nrow(rels.infec2) > 0){
  #           for(j in 1:nrow(rels.infec2)){
  #             if(ids.rels[[i]] == rels.infec2$ID1[j]){ # make diff1
  #               rel.est.numb2 <- c(rel.est.numb2,rels.infec2$ID1[j])
  #               rel.est.time2 <- c(rel.est.time2,rels.infec2$FormTime[j])
  #             }
  #           }
  #         }
  #         rel.est.numb.time2 <- cbind(rel.est.numb2,rel.est.time2)
  #       }
  #       return(rel.est.numb.time2)
  #     }
  #     z2[[k]] <- count.relationships.2(k,rels.infec2)
  #
  #
  #     ## All relationships of the infected individuals: IDs and time of relationship est.
  #
  #     rels.ind.times[[k]] <- rbind(z1[[k]],z2[[k]])
  #
  #
  #     ## number of relationships per individual
  #
  #     rels.per.ind <- vector("list", (nrow(b.dt[[k]])))
  #     w <- vector("list", (nrow(b.dt[[k]])))
  #     v <- vector("list", (nrow(b.dt[[k]])))
  #
  #     for(i in 1:(nrow(b.dt[[k]]))){
  #
  #       rels.per.ind[[i]] <- subset(rels.ind.times[[k]], rels.ind.times[[k]][,1]%in%b.dt[[k]][,2][i])
  #
  #     }
  #
  #     for(j in 1:(nrow(b.dt[[k]]))){
  #       #pers.infec.raw[which(pers.infec.raw$InfectTime <= endpoint),]
  #       d <- as.data.frame(rels.per.ind[[j]])
  #
  #       v[[j]] <- d[which(d$rel.est.time1 >= b.dt[[k]][,1][j]-window), ] # take into account window time for relationships >> sampling time  window (backward)
  #       # we can miss an ongoing relationship if the latter last for a while and is no longer in the time window
  #       # between samling time and window time (backward)
  #       # e.g how many established relationships in the lst pasy 5 yrs whereas the individual established one or more rels in 7 yrs
  #       # they still ongoing or broken
  #     }
  #
  #     # Count relationships per individual before sampling and ongoing
  #
  #     rels.count <- vector()
  #     rels.identif <- vector()
  #     for(i in 1:(nrow(b.dt[[k]]))){
  #       rels.count <- c(rels.count,nrow(v[[i]]))
  #       rels.identif <- c(rels.identif,b.dt[[k]][,2][i])
  #       ids.rels.count <- as.data.frame(cbind(rels.identif, rels.count))
  #     }
  #     count.rels[[k]] <- ids.rels.count
  #   }
  #   return(count.rels)
  # }
  #
  #
  #
  # filter.rels <- count.rels.filter(b.dt = b.dt, datalist = datalist, window = endpoint)
  #
  #
  # # (iii) Age at sampling/removal time: OK
  #
  # age.filter <- function(b.dt, datalist){
  #   age.count <- vector("list", length(seeds.id))
  #   age.sampling <- vector("list", length(seeds.id))
  #   id.age.sampling <- vector("list", length(seeds.id))
  #   pers.dat <- datalist$ptable
  #
  #   # Get first age at the beginning of simulation
  #
  #   for(k in 1:length(seeds.id)){
  #     baseline.dat <- subset(pers.dat, pers.dat$ID%in%b.dt[[k]][,2])
  #     id.TOB <- cbind(baseline.dat$ID,baseline.dat$TOB,b.dt[[k]][,2])
  #     age.TOB.sort <- function(id.TOB){
  #       id.TOB.sort <- NULL
  #       id.sort <- vector()
  #       TOB.sort <- vector()
  #
  #       for(i in 1:nrow(id.TOB)){
  #         for(j in 1:nrow(id.TOB)){
  #           if(id.TOB[,3][i] == id.TOB[,1][j]){
  #             id.sort <- c(id.sort,id.TOB[,3][i])
  #             TOB.sort <- c(TOB.sort,id.TOB[,2][j])
  #           }
  #         }
  #       }
  #       id.TOB.sort <- cbind(id.sort,TOB.sort)
  #       return(id.TOB.sort)
  #     }
  #
  #     c <- age.TOB.sort(id.TOB) # IDs and their age at the begining of the simulation
  #     # (negatif means individuals was there at the begining)
  #     # -12, this person was 12 yaer old at the begining
  #
  #     # the age at the sampling time = - age at begining of sim + sampling time
  #     #
  #     # ifelse(c[,2] <= 0,
  #     #        age.sampling[[k]] <- (b.dt[[k]][,1] - c[,2]),
  #     #        age.sampling[[k]] <- (b.dt[[k]][,1] + c[,2]))
  #
  #     age.sampling[[k]] <- (abs(c[,2]) + b.dt[[k]][,1]) # c > TOB, b.dt > sampling time
  #
  #     id.age.sampling[[k]] <- as.data.frame(cbind(b.dt[[k]][,2],age.sampling[[k]]))
  #   }
  #   return(id.age.sampling)
  # }
  # filter.age <- age.filter(b.dt, datalist)

  # (iv) Viral load at sampling time: interpolation

  # viralload <- vector("list", length(seeds.id))
  #
  # for(i in 1:length(seeds.id)){
  #   IDs <- b.dt[[i]][,2]
  #   viralload.dat <- datalist$vltable
  #   viralload.sub <- subset(viralload.dat, viralload.dat$ID%in%IDs)
  #   viralload[[i]] <- viralload.sub
  # }
  #
  # unique.vl <- vector("list", length(seeds.id))
  #
  # ids.vl <- vector("list", length(seeds.id))
  #
  # for(i in 1:length(seeds.id)){ # seed
  #   IDs <- b.dt[[i]][,2]
  #   for(j in 1:length(IDs)){
  #     diag.time <- b.dt[[i]][,1][j]
  #     v <- subset(viralload[[i]], viralload[[i]]$ID==IDs[j])
  #     time.vl <- v$Time
  #
  #     if(time.vl < diag.time){ # above max time in vtable > take last value of VL
  #       unique.vl[[i]] <- tail(v$Log10VL, n=1)
  #
  #       ids.vl[[i]] <- b.dt[[i]][,2][j]
  #     }
  #
  #     if(time.vl > diag.time){ # below min time in vtable > take first value of VL
  #       unique.vl[[i]] <- v$Log10VL[1]
  #
  #       ids.vl[[i]] <- b.dt[[i]][,2][j]
  #     }
  #
  #
  #   }
  # }


  filter.ls <- vector("list", length(seeds.id))

  filter.unique <- list()
  for(i in 1:length(seeds.id)){
    filter.unique$itimes <- endpoint - dat.recdontime[[i]][,5] # Infection time > seroconversion time OK
    #filter.unique$dtimes <- endpoint - b.dt[[i]][,1] # Sampling/Removal time > Diagnosis time OK
    filter.unique$id <- dat.recdontime[[i]][,1] # Receipients ID > Infected individuals OK
    filter.unique$parent <- dat.recdontime[[i]][,4] # New Donors ID
    #filter.unique$age <- filter.age[[i]][,2] # age at sampling/removal OK
    #filter.unique$rels <- filter.rels[[i]][,2] # number of relationships OK
    #filter.unique$gender <- filter.gender[[i]][,2] # gender OK
    filter.unique$id.orig.recip <- dat.recdontime[[i]][ , 2] # Original ID of the recipient, as in the datalist$ptable
    filter.unique$id.orig.donor <- dat.recdontime[[i]][ , 3] # Original ID of the donor, as in the datalist$ptable




    filter.ls[[i]] <- filter.unique
  }

  return(filter.ls)


  # diagnosis time, > OK
  # age at sampling, > OK
  # the gender of individual, > OK
  # number of relationships, > Ok
  # viral load level at sampling time, > interpolate
  # type of sexual relationship, > Ok
  # infection stage, > infection stage
  # sero-conversion date, > Ok
  # being on treatment > when he (she) went on treatment after seroconversion > treatment time

}

