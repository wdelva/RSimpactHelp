source(file = "~/Documents/MiceABC/R/00-Functions-MSM.R")
destDir <- "/Users/delvaw/Documents/temp"
seedid <- 1

cfg <-list()
cfg$population.msm <- "yes"
cfg$population.numwomen <- 0
cfg$hivtransmission.param.a <- -1
cfg$hivtransmission.param.b <- -90
cfg$hivtransmission.param.c <- 0.5

msm.results <- simpact.run(configParams = cfg,
            destDir = destDir,
            seed = seedid)
datalist.msm <- readthedata(msm.results)


transm.ls <- transmNetworkBuilder.baseline(datalist = datalist.msm,
                                           endpoint = datalist.msm$itable$population.simtime[1])



transm.ls.attrib <- attributes.trans.network(datalist = datalist.msm,
                                             endpoint = datalist.msm$itable$population.simtime[1])

# We start with a transnetwork object (the 3th chain of the 15 chains)
transnetwork <- transm.ls[[3]] # 1 seed considered at a time
# testing with different, unequal dtimes
transnetwork$dtimes <- transnetwork$dtimes + runif(length(transnetwork$dtimes), min = 0, max = 0.68)

transnetw.attr <- transm.ls.attrib[[3]]


# We want to know how many people were infected by each node
# And we want to know when these infections happened
transmlist.df <- dplyr::filter(datalist.msm$etable,
                               eventname == "transmission")
indices.transm <- transmlist.df$p2ID %in% transnetw.attr$id.orig.recip
tranms.times <- c(datalist.msm$itable$hivseed.time[1], transmlist.df$eventtime[indices.transm])
node.ids <- c(transnetw.attr$id.orig.recip[1], transmlist.df$p2ID[indices.transm])

edges.df <- data.frame(from = transnetwork$parent,
                       to = transnetwork$id,
                       inf.time = datalist.msm$itable$population.simtime[1] - transnetwork$itimes)
vertices.df <- data.frame(nodes = unique(c(transnetwork$parent, transnetwork$id)),
                          tranms.times = c(-1, tranms.times))

net <- graph_from_data_frame(d = edges.df, vertices = vertices.df, directed=T)
V(net)$size <- 3

#V(net)$frame.color <- "white"
#V(net)$color <- "orange"
#V(net)$label <- ""
#####E(net)$arrow.size <- .5
#E(net)$label <- round(E(net)$inf.time, 1)


l <- layout_as_tree(net, flip.y = FALSE)
l.horiz <- l[, ncol(l):1]
l.horiz[, 2] <- -l.horiz[, 2]

# layout_with_lgl
# layout_with_kk # Or just a way to initiate the layout
# layout_with_fr
# layout_as_tree
# layout_nicely
plot(net, layout = l)

net.minus.root <- delete_edges(net, 1) %>% delete_vertices(., "-1")
E(net.minus.root)$label <- round(E(net.minus.root)$inf.time, 1)

l <- layout_as_tree(net.minus.root, flip.y = FALSE)
l.horiz <- l[, ncol(l):1]
l.horiz[, 2] <- -l.horiz[, 2]
plot(net.minus.root,
     layout = l.horiz, #layout_as_tree(net, flip.y = FALSE),
     vertex.size = 40,
     vertex.color = "lightblue",
     vertex.label.color = "darkblue",
     rescale = FALSE,
     asp = 1,
     xlim = c(-1, 7),
     ylim = c(-0.5,1.5),
     vertex.label.family = "sans",
     vertex.label.font = 2,
     edge.label.family = "sans",
     edge.label.cex = 0.8,
     edge.label.color = "black",
     edge.label.font = 2,
     edge.arrow.size = 1.5,
     edge.arrow.width = 2,
     edge.width = 20)     #vertex.label = NA)

tranms.times

xlims <- c(-2, datalist.msm$itable$population.simtime[1] - datalist.msm$itable$hivseed.time[1])
ylims <- c(0, nrow(vertices.df) + 1)
l[, 1] <- c(-1, tranms.times)
l[, 2] <- (nrow(vertices.df) + 1) / 2
l[5, 2] <- 5
l[6, 2] <- 3
l[7, 2] <- 2
l[8, 2] <- 11
l[9, 2] <- 5
l[10, 2] <- 9
l[11, 2] <- 3
plot(net,
     rescale = FALSE,
     asp = 0,
     axes = TRUE,
     xlim = xlims,
     ylim = ylims,
     layout = l, vertex.size = 3, vertex.label = NA)

igraph::degree(net, mode = "out")
igraph::degree_distribution(net, mode = "out")


# Now the reconstructed tree

epi.tree <- epi2tree2(transnetwork)
epi.tree$root.edge


tree.dat.full <- simSeq(epi.tree,
                        l = 100,
                        bf = freq,
                        #rootseq = hiv.seq.env.short,
                        type = "DNA",
                        rate = overall.rate)

tree.dat.gtr <-simSeq(epi.tree,
                      l = 1000,
                      Q = c(3.37,14.50,1.44,1.21,14.50,1.00),
                      rate = 10*overall.rate)


fit.ini <- pml(epi.tree, tree.dat.full, k = 4)

fit.ini <- pml(epi.tree, tree.dat.gtr)#, k = 4)



fit <- optim.pml(fit.ini,
                 optNni = FALSE,
                 optBf = FALSE,
                 optQ = FALSE,
                 model = "GTR",
                 control = pml.control(epsilon = 1e-08,
                                       maxit = 10),
                 optGamma = TRUE,
                 optEdge = FALSE,
                 optRate = FALSE,
                 optRooted = FALSE)


fit <- optim.pml(fit.ini,
                 optNni = FALSE,
                 optBf = FALSE,
                 optQ = FALSE,
                 model = "GTR",
                 control = pml.control(epsilon = 1e-09,
                                       maxit = 30),
                 optGamma = FALSE,
                 optEdge = TRUE,
                 optRate = FALSE,
                 optRooted = FALSE)

# fit$tree stores the tree
plot(fit$tree)
nodelabels()

plot(midpoint(fit$tree))

rooted <- root(fit$tree, node = "18", resolve.root = TRUE)
plot(midpoint(rooted))
nodelabels()
test.tree <- rotate(fit$tree, node = c("11"))
test.tree <- rotateConstr(fit$tree, test.tree$tip.label)

plot(midpoint(test.tree, show.node.label=TRUE))

rerooted.fittedtree <- multi2di(fit$tree, random = FALSE)
is.rooted(rerooted.fittedtree)

plot(rerooted.fittedtree)
plot(fit$tree)
all.equal(epi.tree, fit$tree) # We didn't change the tree, only the model for the viral evolution, and the simulated sequences.
tree.dat.full.GTR <- simSeq(fit) # The phyDat object that stores the sequences
plot(fit$tree,
     root.edge = TRUE,
     type = "phylogram", #"fan",
     open.angle = 20,
     rotate.tree = 10,
     show.tip.label = FALSE) #fan
axisPhylo(root.time = tranms.times[2],
          backward = FALSE)#),

#add.scale.bar(x = 0, y = 0, length = 15, lwd = 2, lcol = "blue3", col = "blue3")
tiplabels(text = fit$tree$tip.label,
          frame = "circle",
          bg = "lightblue",
          col = "darkblue")

#pch = 21, bg = gray(1:23/23), cex = 2, adj = 1.4)


plot(epi.tree, root.edge = FALSE) # TRUE

# This tree is the "hyperfull" tree: it's topology is given entirely by the transmission events and
# it does not check whether everybody is still alive at the time of the sequencing.

head(datalist.msm$rtable, 20)
