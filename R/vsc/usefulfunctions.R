
#remove columns 

#df <- subset(df, select = -c(x,y))
#drop = c("x","y")
#df <- df[,!(names(df) %in% drop)]

#new.xdesign.ssd.pc.c <- new.xdesign.ssd.pc.c[order(new.xdesign.ssd.pc.c$ssd.c),]
#inputANDoutput.select <- inputANDoutput.select %>% dplyr::group_by(sim.id) %>% dplyr::summarise_each(funs(mean))

#pairs(inputANDoutput.select[, x.variables[10:14]], col = 1+inputANDoutput.select$is.complete, pch = 16, cex = 2)#,

# ICDCodesfreq_rt$Description <- NA
# iterate_freq <-dim(ICDCodesfreq_rt)[1]
# for (i in 1:iterate_freq){
#   if(!is.na(ICDCodesfreq_rt$Var1[i]) & ICDCodesfreq_rt$Var1[i]!="NA" ) {
#     ICDCodesfreq_rt$Description[i] <- ICDCODES_Mapping$Clinical.syndrome[which(ICDCODES_Mapping$ICDcode==ICDCodesfreq_rt$Var1[i])]}
# }



#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}

#qplot(ICD1_TB_General, data = plotme_now, fill = hiv_status, ylab = "Frequency", xlab = "Top 10 causes of death with HIV status") +
#  scale_x_discrete(labels=addline_format(c("A_2" = "Tuberculosis","J_5" = "Pneumonia","I_25" = "Stroke","N_6" = "Kidney Failure","A_8" = "Other Sepsis","G_5" = "Meningitis", "B_4" = "HIV Infection","A_1" = "Intestinal Infectious" ,"I_1" = "Hypertension" ,"B_9" = "Crypto coccosis"))) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))


