#load packages
pacman::p_load(dplyr, missForest, mice)

dirname <- "/home/trust/Documents/TrustLander-IBM/analysis_round2"


complete.df <- data.frame(read.csv(file = paste0(dirname,"/inputANDoutput.SSE.df.csv"), header = TRUE))

#select the row that have minimal SSE
complet.df.min.sse <- dplyr::filter(complete.df, mean.match == TRUE)

targets.df <- complet.df.min.sse[,1:14]
param.df <- complet.df.min.sse[,41:52]


target.stats <- c(0.015, #gr
                  0.016, 0.043, 0.031, 0.027,  #inc
                  0.008, 0.143, 0.21, 0.47, 0.47, 0.538, #prev
                  0.33, 0.34, #art
                  5) #ad

missing.params <- rep(NA, ncol(param.df))

targets.df <- rbind(targets.df, target.stats)
param.df <- rbind(param.df, missing.params)

complete.wna.df <- cbind(targets.df, param.df)


miss.forest.out <- missForest(complete.wna.df, ntree=500, verbose = TRUE) #number of tree 500
complet.df.wo.na <- miss.forest.out$ximp #get the cleaned data frame

tail(complet.df.wo.na,1)


input.sets <- 100
mice.imputted.out <- mice(complete.wna.df, m=input.sets, maxit = 100, method = 'pmm', seed = 500)
#summary(mice.imputted.out)
complete.mice.wo.na <- complete(mice.imputted.out,1)
#collect all the imputed sets
mice.df.out <- tail(complete.mice.wo.na,1)
for (i in 2:input.sets){
  complete.mice.wo.na <- complete(mice.imputted.out,i)
  mice.df.out <- rbind(mice.df.out, tail(complete.mice.wo.na,1))
}




