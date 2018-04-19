#' Build a cumulative relationships network in a given time interval
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param duration Only relationships events that took place between that time interval.
#' @return A graph object for relationships network data.
#' @examples
#' tw.net <- time.wind.network(datalist = datalist, duration = c(1,15))
#' plot.igraph(tw.net,layout=layout.fruchterman.reingold)
#' print(tw.net, full = igraph_opt("print.full"))
#' vertex_attr(tw.net, index = V(tp.net)) # to print the attributes
#' @import igraph
#'

time.wind.network <- function(datalist = datalist, duration = c(1,15)){

  # NETWOTK

  # Relationships table
  rels <- as.data.frame(datalist$rtable)

  # keep those variables
  newdat <- c("ID1", "ID2", "FormTime", "DisTime", "AgeGap")

  # new relationship table with "ID1", "ID2" and "FormTime"
  newrels <- rels[newdat]

  # relationships between c(1,15)


  newrels.tt <- newrels[which(newrels$FormTime <= duration[2] & newrels$FormTime >= duration[1]),]

  id <- c("ID1", "ID2")
  id.net.contt <- newrels.tt[id]
  id.nett <- id.net.contt  # for network



  library(igraph)

  id.nett[,1] <- as.character(id.nett[,1])
  id.nett[,2] <- as.character(id.nett[,2])

  matgraph.t <- as.matrix(id.nett)

  rels.graph.t <- graph.edgelist(matgraph.t[,1:2])

  # Weight of the edges: time of a relationship

  # lasttime = (newrels$DisTime - newrels$FormTime)

  # We then add the edge weights to this network by assigning an edge attribute called 'weight'.
  # E(rels.graph.t)$weight=as.numeric(lasttime)


  # NODE ATTRIBUTES : gender and hiv status

  # 1. id and gender

  id.all <- unique(c(id.net.contt$ID1,id.net.contt$ID2)) # use id.net.contt otherwise id.nett has become matrix of characters

  pers <- as.data.frame(datalist$ptable) # person table

  # Search for gender

  id.gen.t <- vector()
  Gender.t <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers$ID)){
      if(id.all[i] == pers$ID[j]){
        id.gen.t <- c(id.gen.t,id.all[i])
        Gender.t <- c(Gender.t,pers$Gender[j])
      }
    }
  }

  per.Gender.t <- cbind(id.gen.t,Gender.t) # 0 man & 1 woman



  # 2. id and hiv status

  #  person table of infected individuals
  pers.hiv.raw <- as.data.frame(datalist$ptable[InfectType !="-1"])

  #  person table of infected individuals before time endtime
  pers.hiv.t <- pers.hiv.raw[which(pers.hiv.raw$InfectTime <= duration[2]),]


  # Search for HIV status

  id.hiv.t <- vector()
  hivStatus.t <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers.hiv.t$ID)){
      if(id.all[i] == pers.hiv.t$ID[j]){
        id.hiv.t <- c(id.hiv.t,id.all[i])
        hivStatus.t <- c(hivStatus.t,1) # 1 means that the person is HIV positive
      }
    }
  }

  per.hivstatus.t <- cbind(id.hiv.t, hivStatus.t)



  # V(rels.graph.t)$name


  # Add attributes

  # This code says to create a vertex attribute called "Gender" by extracting the value
  # of the column per.Gender.t[,2] of individuals gender and match it with appropriate ID of individual

  V(rels.graph.t)$Gender=as.character(per.Gender.t[,2][match(V(rels.graph.t)$name, per.Gender.t[,1])])
  V(rels.graph.t)$Gender

  # This code says to create a vertex attribute called "HIVstatus" by extracting the value
  # of the column per.hivstatus.t[,2] of individuals HIV status and match it with appropriate ID of individual

  V(rels.graph.t)$HIVstatus = as.character(per.hivstatus.t[,2][match(V(rels.graph.t)$name, per.hivstatus.t[,1])])
  V(rels.graph.t)$HIVstatus

  V(rels.graph.t)$color=V(rels.graph.t)$HIVstatus
  V(rels.graph.t)$color=gsub("1","red",V(rels.graph.t)$color) # hiv + individuals
  V(rels.graph.t)$color=gsub("NA","yellow",V(rels.graph.t)$color) # hiv - individuals

  V(rels.graph.t)$shape=V(rels.graph.t)$Gender
  V(rels.graph.t)$shape=gsub("0","circle",V(rels.graph.t)$shape) # women individuals
  V(rels.graph.t)$shape=gsub("1","square",V(rels.graph.t)$shape) # men individuals


  return(rels.graph.t)


}
