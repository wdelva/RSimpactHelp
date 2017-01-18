

#inputANDoutput.select <- inputANDoutput.select %>% dplyr::group_by(sim.id) %>% dplyr::summarise_each(funs(mean))


#
#pairs(inputANDoutput.select[, x.variables[10:14]], col = 1+inputANDoutput.select$is.complete, pch = 16, cex = 2)#,





