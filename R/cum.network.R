#' Build a cumulative relationships network - an igraph object
#'
#' Build a cumulative relationships network from the begining of the simulation untill a given time point.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param endtime Only relationships events that took place before this point in simulation time, are captured in the output.
#' @return A graph object for relationships network data with difference between gender and HIV status
#' men are represented by square, women by sphere, HIV positive with red color and HIV negative with bue clolor
#' @examples
#' c.net <- cum.network(datalist = datalist, endtime = 40)
#' plot.igraph(c.net,layout=layout.fruchterman.reingold)
#' print(c.net, full = igraph_opt("print.full")) # to print the attributes
#'
#'@import igraph

cum.network <- function(datalist = datalist, endtime = 15){

  # NETWOTK

  # Relationships table
  rels <- as.data.frame(datalist$rtable)

  # keep those variables
  newdat <- c("ID1", "ID2", "FormTime", "DisTime", "AgeGap")

  # new relationship table with "ID1", "ID2" and "FormTime"
  newrels <- rels[newdat]

  # relationships before time endtime

  newrels.t <- newrels[which(newrels$FormTime <= endtime),]

  id <- c("ID1", "ID2")
  id.net.cont <- newrels.t[id]
  id.net <- id.net.cont  # for network


  id.net[,1] <- as.character(id.net[,1])
  id.net[,2] <- as.character(id.net[,2])

  matgraph <- as.matrix(id.net)

  rels.graph <- graph.edgelist(matgraph[,1:2])

  # Weight of the edges: time of a relationship

  # lasttime = (newrels$DisTime - newrels$FormTime)

  # We then add the edge weights to this network by assigning an edge attribute called 'weight'.
  # E(rels.graph)$weight=as.numeric(lasttime)



  # NODE ATTRIBUTES : gender and hiv status

  # 1. id and gender

  id.all <- unique(c(id.net.cont$ID1,id.net.cont$ID2)) # use id.net.cont otherwise id.net has become matrix of characters

  pers <- as.data.frame(datalist$ptable) # person table

  id.gen.p <- vector()
  Gender.p <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers$ID)){
      if(id.all[i] == pers$ID[j]){
        id.gen.p <- c(id.gen.p,id.all[i])
        Gender.p <- c(Gender.p,pers$Gender[j])
      }
    }
  }

  per.Gender <- cbind(id.gen.p,Gender.p) # 0 man & 1 woman



  # 2. id and hiv status

  #  person table of infected individuals
  pers.hiv.raw <- as.data.frame(datalist$ptable[InfectType !="-1"])

  #  person table of infected individuals before time endtime
  pers.hiv <- pers.hiv.raw[which(pers.hiv.raw$InfectTime <= endtime),]

  id.hiv.p <- vector()
  hivStatus.p <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers.hiv$ID)){
      if(id.all[i] == pers.hiv$ID[j]){
        id.hiv.p <- c(id.hiv.p,id.all[i])
        hivStatus.p <- c(hivStatus.p,4) # 4 means that the person is HIV positive
      }
    }
  }

  per.hivstatus <- cbind(id.hiv.p, hivStatus.p)


  # V(rels.graph)$name

  # Add attributes

  # This code says to create a vertex attribute called "Gender" by extracting the value
  # of the column per.Gender.t[,2] of individuals gender and match it with appropriate ID of individual

 V(rels.graph)$Gender=as.character(per.Gender[,2][match(V(rels.graph)$name, per.Gender[,1])])
 V(rels.graph)$Gender

 # This code says to create a vertex attribute called "HIVstatus" by extracting the value
 # of the column per.hivstatus.t[,2] of individuals HIV status and match it with appropriate ID of individual

 V(rels.graph)$HIVstatus = as.character(per.hivstatus[,2][match(V(rels.graph)$name, per.hivstatus[,1])])
 V(rels.graph)$HIVstatus

 V(rels.graph)$color=V(rels.graph)$HIVstatus
 V(rels.graph)$color=gsub("4","red",V(rels.graph)$color) # hiv + individuals
 V(rels.graph)$color=gsub("NA","blue",V(rels.graph)$color) # hiv - individuals

 V(rels.graph)$shape=V(rels.graph)$Gender
 V(rels.graph)$shape=gsub("0","sphere",V(rels.graph)$shape) # women individuals
 V(rels.graph)$shape=gsub("1","square",V(rels.graph)$shape) # men individuals


 # plot.igraph(rels.graph,layout=layout.fruchterman.reingold)

 return(rels.graph)


}

