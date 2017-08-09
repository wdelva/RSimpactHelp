#' Summarise the properties of a network

properties_network <- function(graph = graph){

  test.graph <- graph

  if(length(E(test.graph))){

    # Additionally I will use the NetIndices package,
    # since its function "GenInd()" outputs several network properties
    library(NetIndices)

    # The "GenInd()" function requires an input of an adjacency matrix
    test.graph.adj<-get.adjacency(test.graph,sparse=F)
    # in older versions of igraph the default was sparse=F,
    # but now you must specify, other wise you get a matrix of 1s and .s

    test.graph.properties<-GenInd(test.graph.adj)

    # The function output consists of 10 network properties.
    # I will consider five of them here:

    test.graph.properties$N            #number of nodes
    test.graph.properties$Ltot # links
    test.graph.properties$LD        #link density (average # of links per node)



    # The degree of a node refers to the number of links associated with a node.
    # Degree can be measured as the links going in ("in degree"), out ("out degree"), or both.
    # The degree() function takes a graph input and gives the degree of specified nodes.
    # With the argument "v=V(graph)" you tell the function to give the degree of all nodes in the graph,
    # while the "mode" argument specifies in, out, or both.

    # plot and fit the power law distribution
    fit_power_law = function(graph) {
      # calculate degree
      d = degree(graph, mode = "all")
      dd = degree.distribution(graph, mode = "all", cumulative = T)
      degree = 1:max(d)
      probability = dd[-1]
      # delete blank values
      nonzero.position = which(probability != 0)
      probability = probability[nonzero.position]
      degree = degree[nonzero.position]
      reg = lm(log(probability) ~ log(degree))
      cozf = coef(reg)
      power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
      alpha = -cozf[[2]]
      beta = -cozf[[1]]
      R.square = summary(reg)$r.squared
      # print(paste("Alpha =", round(alpha, 3)))
      # print(paste("Beta =", round(beta, 3)))
      # print(paste("R square =", round(R.square, 3)))
      # # plot
      # plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
      #      col = 1, main = "Degree Distribution")
      # curve(power.law.fit, col = "red", add = T, n = length(d))
      ret <- list()
      ret$Alpha <- round(alpha, 3)
      ret$Beta <- round(beta, 3)
      ret$RSquare <- round(R.square, 3)

      return(ret)
    }


    fit.degree.net <- fit_power_law(test.graph)

    # Clustering coefficient is the proportion of
    # a nodes neighbors that can be reached by other neighbors
    # in igraph this property is apparently called "transitivity"

    glob.clust <- transitivity(test.graph)
    Nodes <- test.graph.properties$N
    links <- test.graph.properties$Ltot
    linksdensity <- test.graph.properties$LD
    alpha <- fit.degree.net$Alpha
    beta <- fit.degree.net$Beta
    Rsquare <- fit.degree.net$RSquare
    # gives the clustering coefficient of the whole network
    net.properties <- as.data.frame(cbind(Nodes,links,linksdensity,
                                          alpha,beta,Rsquare,
                                          glob.clust))

    return(net.properties)

  }else{
    return("NA")
  }

}
