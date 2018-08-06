#' A function that returns age mixing patterns quantities in transmission clusters
#' in scenarios where individuals are missing at completly at random
#' @param simpact.trans.net a list of transmission networks produced by \code{\link{transm.network.builder}}
#' @param work.dir working directory
#' @param dirfasttree directory where is the fastTree tool
#' @param sub.dir.rename subdurectory required when we have to run more than one simulations
#' @param limitTransmEvents Number of minimum transmission events to be considered in each transmission networks
#' @param timewindow Time interval
#' @param seq.cov Percentage of individulas considered for this transmission pattern scenario
#' @param age.group.15.25 age group between 15 and 25 years old
#' @param age.group.25.40 age group between 25 and 40 years old
#' @param age.group.40.50 age group between 40 and 50 years old
#' @return a vector of number of men and women in different age group, number of transmissions within all age groups, and mean and SD of age different between infectors and infectees
#' @examples
#' w <- phylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#'                                work.dir = work.dir,
#'                                dirfasttree = dirfasttree,
#'                                sub.dir.rename = sub.dir.rename,
#'                                limitTransmEvents = 7,
#'                                timewindow = c(30,40),
#'                                seq.cov = 70,
#'                                age.group.15.25 = c(15,25),
#'                                age.group.25.40 = c(25,40),
#'                                age.group.40.50 = c(40,50))

#'
#' @importFrom magrittr %>%
#' @import dplyr filter
#' @export
#'


# infection with phylo - sampling time

phylo.CAR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
                                        work.dir = work.dir,
                                        dirfasttree = dirfasttree,
                                        sub.dir.rename = sub.dir.rename,
                                        limitTransmEvents = 7,
                                        timewindow = c(30,40),
                                        seq.cov = 70,
                                        age.group.15.25 = c(15,25),
                                        age.group.25.40 = c(25,40),
                                        age.group.40.50 = c(40,50)){



  seeds.id <- length(simpact.trans.net)

  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)

  for(i in 1:seeds.id){

    transm.age.df.ic <- as.data.frame(simpact.trans.net[[i]])

    age.samp.Rec <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBRec
    age.samp.Don <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBDon

    transm.age.i <- cbind(transm.age.df.ic, age.samp.Rec, age.samp.Don)

    new.transm.tab[[i]] <- transm.age.i

  }


  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow,
                             seq.cov = seq.cov,
                             age.limit = age.group.40.50[2])

  if(length(mCAr.IDs)>5){


    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mCAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))


    mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                          sub.dir.rename = sub.dir.rename,
                                                          fasttree.tool = "FastTree",
                                                          calendar.dates = "samplingtimes.all.csv",
                                                          simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                          count.start = 1977,
                                                          endsim = 40,
                                                          clust = FALSE)

    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")))



    # run ClusterPicker

    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))

    # Read clusters' files

    d <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                    all.files = FALSE,
                    full.names = FALSE, recursive = FALSE)





    # id of people who got infection by seed event: seeds.id
    trans.network <- new.transm.tab

    seeds.id <- length(trans.network)


    ID.select <- vector() # ID of selected transmission network
    ID.select.count <- vector() # number of individuals in these networks

    for (i in 1: seeds.id) {


      trans.network.i <- as.data.frame(trans.network[[i]])

      if(nrow(trans.network.i)>=limitTransmEvents){


        ID.select <- c(ID.select, i)
        ID.select.count <- c(ID.select.count, nrow(trans.network.i))

      } # X if

    } # Y for


    infectionTable <- vector("list", length(ID.select))

    for(j in 1:length(ID.select)){

      p <- ID.select[j]

      trans.network.i <- as.data.frame(trans.network[[p]])

      trans.network.i <- trans.network.i[-1,]

      id.lab <- paste0(p,".",trans.network.i$id,".C")

      trans.network.i$id.lab <- id.lab

      infectionTable[[p]] <- trans.network.i
    }


    infecttable <- rbindlist(infectionTable)



    transm.df <- infecttable
    # transm.df <- agemixing.trans.df(trans.network = simpact.trans.net,
    #                                 limitTransmEvents = 7)


    # define to filter pairings
    sort.partners.fun.phylo <- function(partner.table = partner.table){ # for receivers

      # age and gender structured receiver individuals

      num.15.25.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                error=function(e) return(NULL))

      num.15.25.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                  error=function(e) return(NULL))


      num.25.40.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                error=function(e) return(NULL))

      num.25.40.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                  error=function(e) return(NULL))


      num.40.50.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                error=function(e) return(NULL))

      num.40.50.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                  error=function(e) return(NULL))

      comb = function(n, x) {

        if(n>=x){
          va <- factorial(n) / factorial(n-x) / factorial(x)
        }else{
          va <- factorial(x) / factorial(x-n) / factorial(n)
        }
        return(va)
      }


      # Possibles pairings


      pairs.15.25.men.women.15.25 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.15.25.women)), # C(1,n)
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.15.25 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.15.25.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.40.50.men.women.15.25 <- comb(1, nrow(num.40.50.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.15.25.women)),
      #                                       error=function(e) return(NA))
      #

      pairs.15.25.men.women.25.40 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.25.40 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.40.50.men.women.25.40 <- comb(1, nrow(num.40.50.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #


      pairs.15.25.men.women.40.50 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.40.50 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.40.50.men.women.40.50 <-  comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #


      pairings.al <- c(pairs.15.25.men.women.15.25, pairs.15.25.men.women.25.40, pairs.15.25.men.women.40.50,
                       pairs.25.40.men.women.15.25, pairs.25.40.men.women.25.40, pairs.25.40.men.women.40.50,
                       pairs.40.50.men.women.15.25, pairs.40.50.men.women.25.40, pairs.40.50.men.women.40.50)

      men.women  <- c(nrow(num.15.25.men), nrow(num.15.25.women),
                      nrow(num.25.40.men), nrow(num.25.40.women),
                      nrow(num.40.50.men), nrow(num.40.50.women))


      val.names <- c("num.men.15.25", "num.women.15.25",
                     "num.men.25.40", "num.women.25.40",
                     "num.men.40.50", "num.women.40.50",

                     "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                     "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                     "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50")

      num.ind.pairings.al <- c(men.women, pairings.al)


      names(num.ind.pairings.al) <- val.names


      return(num.ind.pairings.al)


    }

    if(length(d)>=1){

      pairings.clust.tab.list <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster


      for (i in 1:length(d)) {

        clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
        #size <- c(size, nrow(clus.read))

        transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster

        pairings.clust.tab <- sort.partners.fun.phylo(partner.table = transm.df.cl)


        # Age difference statistics #
        #############################
        AD <- abs(abs(transm.df.cl$TOBDon) - abs(transm.df.cl$TOBRec))
        mean.AD <- mean(AD)
        med.AD <- median(AD)
        sd.AD <- sd(AD)

        # Mixed effect models #
        #######################
        # fit.agemix.trans.women <- fit.agemix.trans.women(datatable = data.transm.agemix)
        # fit.agemix.trans.men <- fit.agemix.trans.men(datatable = data.transm.agemix)

        AD.stat <- c(mean.AD, med.AD, sd.AD)

        pairings.clust.tab.AD <- c(pairings.clust.tab, AD.stat)

        pairings.clust.tab.AD <- as.numeric(pairings.clust.tab.AD)

        val.names <- c("num.men.15.25", "num.women.15.25",
                       "num.men.25.40", "num.women.25.40",
                       "num.men.40.50", "num.women.40.50",

                       "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                       "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                       "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",

                       "mean.AD", "median.AD", "sd.AD")

        names(pairings.clust.tab.AD) <- val.names


        pairings.clust.tab.list[[i]] <- pairings.clust.tab.AD

      }

    }else{

      val.names <- c("num.men.15.25", "num.women.15.25",
                     "num.men.25.40", "num.women.25.40",
                     "num.men.40.50", "num.women.40.50",

                     "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                     "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                     "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",

                     "mean.AD", "median.AD", "sd.AD")

      clust.stat.table <- rep(NA, length(val.names))

      names(clust.stat.table) <- val.names


    }





    clust.stat.table <- as.data.frame(do.call(rbind, pairings.clust.tab.list)) # data.table & data.frame

  }else{

    val.names <- c("num.men.15.25", "num.women.15.25",
                   "num.men.25.40", "num.women.25.40",
                   "num.men.40.50", "num.women.40.50",

                   "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                   "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                   "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",

                   "mean.AD", "median.AD", "sd.AD")

    clust.stat.table <- rep(NA, length(val.names))

    names(clust.stat.table) <- val.names

  }

  return(clust.stat.table)

}


# d <- phylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#                                  limitTransmEvents = 7,
#                                  timewindow = c(30,40),
#                                  seq.cov = 70,
#                                  # seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
#                                  age.group.15.25 = c(15,25),
#                                  age.group.25.40 = c(25,40),
#                                  age.group.40.50 = c(40,50))





