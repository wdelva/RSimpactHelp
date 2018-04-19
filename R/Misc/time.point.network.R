#' Build current relationships network at a given point of time
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param time Only relationships present at that time point.
#' @return A graph object for relationships network data.
#' Men are represented by square, women by sphere, HIV positive with red color and HIV negative with bue clolor
#' @examples
#' tp.net <- time.point.network(datalist = datalist, time = 30)
#' plot.igraph(tp.net,layout=layout.fruchterman.reingold)
#' print(tp.net, full = igraph_opt("print.full"))
#' vertex_attr(tp.net, index = V(tp.net)) # to print the attributes
#'
#'
#'@import igraph
#'

time.point.network <- function(datalist = datalist, time = 30){

  # NETWOTK

  # Relationships table
  rels <- as.data.frame(datalist$rtable)

  # keep those variables
  newdat <- c("ID1", "ID2", "FormTime", "DisTime", "AgeGap")

  # new relationship table with "ID1", "ID2" and "FormTime"
  newrels <- rels[newdat]

  # How long last a relationship

  lasttime <- (newrels$DisTime - newrels$FormTime)

  # current relationships at time = 30, those with DistTime >= time (ongoing relationships)
  # FormTime <= time and DisTime >= time

  condform.rels <- rels[which(rels$FormTime <= time),]
  conddiss.rels <- condform.rels[which(condform.rels$DisTime >= time),]

  newrels.tp <- conddiss.rels


  id <- c("ID1", "ID2")
  id.net.conttp <- newrels.tp[id]
  id.nettp <- id.net.conttp  # for network



  library(igraph)

  id.nettp[,1] <- as.character(id.nettp[,1]) # remember to add -1 the universal donor which is in contact with seeds
  id.nettp[,2] <- as.character(id.nettp[,2])

  matgraph.tp <- as.matrix(id.nettp)

  rels.graph.tp <- graph.edgelist(matgraph.tp[,1:2])

  # Weight of the edges: time of a relationship

  # lasttime = (newrels$DisTime - newrels$FormTime)

  # We then add the edge weights to this network by assigning an edge attribute called 'weight'.
  # E(rels.graph.tp)$weight=as.numeric(lasttime)



  # NODE ATTRIBUTES : gender and hiv status

  # 1. id and gender

  id.all <- unique(c(id.net.conttp$ID1,id.net.conttp$ID2)) # use id.net.cont otherwise id.net has become matrix of characters

  pers <- as.data.frame(datalist$ptable) # person table

  id.gen.tp <- vector()
  Gender.tp <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers$ID)){
      if(id.all[i] == pers$ID[j]){
        id.gen.tp <- c(id.gen.tp,id.all[i])
        Gender.tp <- c(Gender.tp,pers$Gender[j])
      }
    }
  }

  per.Gender.tp <- cbind(id.gen.tp,Gender.tp) # 0 man & 1 woman



  # 2. id and hiv status
  #  person table of infected individuals
  pers.hiv.raw <- as.data.frame(datalist$ptable[InfectType !="-1"])

  #  person table of infected individuals before time endtime
  pers.hiv.tp <- pers.hiv.raw[which(pers.hiv.raw$InfectTime <= time),]

  id.hiv.tp <- vector()
  hivStatus.tp <- vector()
  for(i in 1:length(id.all)){
    for(j in 1:length(pers.hiv.tp$ID)){
      if(id.all[i] == pers.hiv.tp$ID[j]){
        id.hiv.tp <- c(id.hiv.tp,id.all[i])
        hivStatus.tp <- c(hivStatus.tp,1) # 1 means that the person is HIV positive
      }
    }
  }

  per.hivstatus.tp <- cbind(id.hiv.tp, hivStatus.tp)


  # Add attributes

  # This code says to create a vertex attribute called "Gender" by extracting the value
  # of the column per.Gender.t[,2] of individuals gender and match it with appropriate ID of individual

  V(rels.graph.tp)$Gender=as.character(per.Gender.tp[,2][match(V(rels.graph.tp)$name, per.Gender.tp[,1])])
  V(rels.graph.tp)$Gender

  # This code says to create a vertex attribute called "HIVstatus" by extracting the value
  # of the column per.hivstatus.t[,2] of individuals HIV status and match it with appropriate ID of individual

  V(rels.graph.tp)$HIVstatus = as.character(per.hivstatus.tp[,2][match(V(rels.graph.tp)$name, per.hivstatus.tp[,1])])
  V(rels.graph.tp)$HIVstatus

  V(rels.graph.tp)$color=V(rels.graph.tp)$HIVstatus
  V(rels.graph.tp)$color=gsub("1","red",V(rels.graph.tp)$color) # hiv + individuals
  V(rels.graph.tp)$color=gsub("NA","blue",V(rels.graph.tp)$color) # hiv - individuals

  V(rels.graph.tp)$shape=V(rels.graph.tp)$Gender
  V(rels.graph.tp)$shape=gsub("0","circle",V(rels.graph.tp)$shape) # women individuals
  V(rels.graph.tp)$shape=gsub("1","square",V(rels.graph.tp)$shape) # men individuals


  return(rels.graph.tp)

}
