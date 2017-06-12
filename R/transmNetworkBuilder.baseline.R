#' Build a list with HIV transmission networks data.
#'
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param endtime Only transmission events that took place before this point in simulation time, are captured in the output.
#' @return a list with the transmission network data, as required by the epi2tree function ( each of them considering same sampling/removing date).
#' @examples
#' transm.ls <- transm.network.builder(datalist = datalist,endpoint = 40)
#' @note
#' transm.ls[[1]] is the transmission network (epi object) of seed 1
#' @import igraph

# Build a transmission network data per seed to be handled by epi2tree function of expotree package

transm.network.builder <- function(datalist = datalist, endpoint = 40){

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
        transNet$itimes <- rev(dat.recdontime[[i]][,5]-10)
        transNet$dtimes <- rep(0,length(dat.recdontime[[i]][,5]))
        transNet$id <- dat.recdontime[[i]][,1]
        transNet$parent <- dat.recdontime[[i]][,4]

        transm.ls[[i]] <- transNet
    }

    return(transm.ls)

}


