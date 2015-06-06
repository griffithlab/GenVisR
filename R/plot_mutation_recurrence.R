#' plot mutation burden
#' 
#' plot a barchart showing mutations per MB
#' @name plot_mutation_recurrence
#' @param data_frame a data frame in MAF format
#' @param coverage_space an integer specifying the coverage space in base pairs from which a mutation could occur
#' @return a ggplot object

plot_mutation_recurrence <- function(data_frame, coverage_space)
{
  #############################################################################################################################
  #################################### Function to plot the top margin barplot ################################################
  #############################################################################################################################
  
  # Add in mutations per MB calculation
  data_frame$mutation_per_MB <- (data_frame$mutation_total * (37500/1599))/coverage_space * 1000000
  temp <- data_frame[-which(data_frame$sample == 'H_KU-1240-1226244' | data_frame$sample == 'H_KU-1827-1307712' | data_frame$sample == 'H_KU-2035-1226280' | data_frame$sample == 'H_KU-2064-1226283' | data_frame$sample == 'H_KU-2095-1226285' | data_frame$sample == 'H_KU-2616-1226315' | data_frame$sample == 'H_KU-2664-1320548' | data_frame$sample == 'H_KU-2829-1226367' | data_frame$sample == 'H_KU-2859-1226000' | data_frame$sample == 'H_KU-3069-1226489' | data_frame$sample == 'H_KU-3156-1226535' | data_frame$sample == 'H_KU-3163-1226021' | data_frame$sample == 'H_KU-3274-1226545' | data_frame$sample == 'H_KU-3286-1226025' | data_frame$sample == 'H_KU-3381-1226037' | data_frame$sample == 'H_KU-3689-1226056' | data_frame$sample == 'H_KU-3742-1226061' | data_frame$sample == 'H_KU-3767-1226066' | data_frame$sample == 'H_KU-3769-1226067' | data_frame$sample == 'H_KU-3899-1226074' | data_frame$sample == 'H_KU-3909-1226076' | data_frame$sample == 'H_KU-4032-1226086' | data_frame$sample == 'H_KU-4150-1226090' | data_frame$sample == 'H_KU-4156-1226091' | data_frame$sample == 'H_KU-4250-1226142' | data_frame$sample == 'H_KU-4289-1226146' | data_frame$sample == 'H_KU-4306-1226150' | data_frame$sample == 'H_KU-436-1226195' | data_frame$sample == 'H_KU-44-1226160' | data_frame$sample == 'H_KU-4407-1226339' | data_frame$sample == 'H_KU-4514-1226132' | data_frame$sample == 'H_KU-4521-1226134' | data_frame$sample == 'H_KU-4667-1226110' | data_frame$sample == 'H_KU-739-1307702' | data_frame$sample == 'H_LV-2280-1320495'),]
  tempA <- temp[which(temp$trv_type == 'Synonymous'),]
  tempB <- temp[which(temp$trv_type == 'Non Synonymous'),]
  tempC <- tempA$mutation_per_MB + tempB$mutation_per_MB
  print(summary(tempC))
  
  #print(data_frame)
  
  # Alter GGplot2 Theme 
  theme <- theme(axis.title.y=element_text(size=18), panel.border =  element_blank(), axis.line =  element_line(), panel.background=element_rect(fill='white'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_text(size=18), legend.text=element_text(size=16))
  
  # Add Legend
  legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c("Synonymous", "Non Synonymous"), drop=FALSE)
  
  # add y label
  y_label <- ylab('Mutations per MB')
  
  # ggplot2 call
  p1 <- ggplot(data_frame, aes(x=sample, y=mutation_per_MB, fill=trv_type)) + geom_bar(stat='identity', alpha=.75, width=1) + theme + y_label + legend
  
  return(p1)
}