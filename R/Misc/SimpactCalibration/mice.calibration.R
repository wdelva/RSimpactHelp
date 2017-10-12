#loading packages
pacman::p_load(dplyr, mice)

dirname <- "~/Documents/GIT_Projects/SimulationParameterSetup"

comp <- "lin" #lin #mac #chpc #gent
#dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"
#dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"

#simpact simulation traing dataset
training.df <- "1rep.mean.adj.sse.TrainingSet.csv"

#Read the saved training set
complete.df.wna <- data.frame(read.csv(file = paste0(dirname,"/",training.df), header = TRUE))

mice.calibration <- function(df.wna, method.id = 1, out.sets = 5, iter.num = 5, seed.id = 2525){

  #selected methods used
  cal.list <- c("pmm", "norm", "norm.boot")

  mice.imputted.df <- mice(df.wna, m=out.sets, maxit = iter.num,
                           method = cal.list[method.id], seed = seed.id)

  #collect all the imputed sets
  mice.out.imputted.df <- tail(mice::complete(mice.imputted.df,1),1)

  for (i in 2:out.sets){
    complete.mice.wo.na <- mice::complete(mice.imputted.df,i)
    mice.out.imputted.df <- rbind(mice.out.imputted.df, tail(complete.mice.wo.na,1))
  }

  #remove all the rows with out of range imputed files
  mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.man.dist.normal.mu >0)
  mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.woman.dist.normal.mu >0)
  mice.out.imputted.df <- subset(mice.out.imputted.df, person.art.accept.threshold.dist.fixed.value < 0.9)
  mice.out.imputted.df <- subset(mice.out.imputted.df, person.agegap.man.dist.normal.sigma > 0)
  mice.out.imputted.df <- subset(mice.out.imputted.df, hivtransmission.param.f1 >0)

  mice.out.imputted.df$method <- cal.list[method.id]

  write.csv(mice.out.imputted.df, file = paste0(dirname,"/mice.",cal.list[method.id],".",training.df),
            row.names = FALSE)

  return(mice.out.imputted.df)

}

#check time
mice.startt <- as.numeric(proc.time()[3])
mice.calibration <- mice.calibration(complete.df.wna, method.id = 1)
mice.end.pmm <- as.numeric(proc.time()[3]) - mice.startt

