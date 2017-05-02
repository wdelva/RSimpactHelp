
#remove columns

#df <- subset(df, select = -c(x,y))
#drop = c("x","y")
#df <- df[,!(names(df) %in% drop)]

#new.xdesign.ssd.pc.c <- new.xdesign.ssd.pc.c[order(new.xdesign.ssd.pc.c$ssd.c),]
#inputANDoutput.select <- inputANDoutput.select %>% dplyr::group_by(sim.id) %>% dplyr::summarise_each(funs(mean))

#pairs(inputANDoutput.select[, x.variables[10:14]], col = 1+inputANDoutput.select$is.complete, pch = 16, cex = 2)#,

#rename col
#names(df)[names(df) == 'old.var.name'] <- 'new.var.name'

#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}

#qplot(ICD1_TB_General, data = plotme_now, fill = hiv_status, ylab = "Frequency", xlab = "Top 10 causes of death with HIV status") +
#  scale_x_discrete(labels=addline_format(c("A_2" = "Tuberculosis","J_5" = "Pneumonia","I_25" = "Stroke","N_6" = "Kidney Failure","A_8" = "Other Sepsis","G_5" = "Meningitis", "B_4" = "HIV Infection","A_1" = "Intestinal Infectious" ,"I_1" = "Hypertension" ,"B_9" = "Crypto coccosis"))) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#visualise the correlation matrix
#The asterisks indicate the significance levels of the correlations.
#Each significance level is associated to a symbol :
#p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
#library(PerformanceAnalytics)
#chart.Correlation(mydata, histogram=TRUE, pch=19)

#rle produces two vectors (lengths and values ). The length of values vector gives you the number of unique values.
#x<-c(1,2,3,1,2,3,4,6)
#length(rle(sort(x))$values)

#values in a list
#vc <- c(1,2,3,10)
#dt[dt$id %in% vc,]
#dplyr::filter(df, fct %in% vc)


