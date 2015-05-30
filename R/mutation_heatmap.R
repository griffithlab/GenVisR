#' Plot a mutation landscape
#' 
#' Plot a mutation landscape plot for a cohort in an annotation file
#' @name mutation_heatmap
#' @param data_frame a data frame in annotation format
#' @param recurrence_cutoff an integer value to remove genes that do not have x number of mutations
#' @param grid a boolean value to overlay a grid on the primary plot
#' @param label_x a boolean value to plot samples on the x axis
#' @param title a character string for the plot title
#' @param gene_label_size an integer specifying the size of labels on Y axis
#' @param coverage_space an integer specifying the size in bp of the genome covered from which mutations could be called
#' @param file_type a character string specifying the file format of the data frame, one of "TGI", "MAF"
#' @param genes a character vector specifying genes to plot
#' @return a grob for plotting
#' @export

mutation_heatmap <- function(data_frame, recurrence_cutoff = 0, grid = TRUE, label_x = FALSE, title ='NULL', gene_label_size=8, coverage_space=63564965, file_type='TGI', genes=NULL)
{
  ############################################################################################
  ######## Function to create a mutation heatmap given a file in TGI annotation format #######
  ############################################################################################
  require(ggplot2)
  
  if(toupper(file_type) == toupper('MAF'))
  {
    data_frame <- MAF_to_anno(data_frame)
  }
  
  # Extract columns from annotation format needed for a mutation heatmap
  data_frame <- data_frame[,c('sample', 'gene', 'trv_type')]
  
  # add in a count of mutations at the sample level before anything is stripped out and save for mutation recurrence plot
  data_frame2 <- add_mutation_counts(data_frame)
  
  # Remove trv_type that are not the most deleterious for a given gene/sample
  data_frame <- hiearchial_remove_trv_type(data_frame, file_type=file_type)
  
  # reorder the genes based on frequency of mutations in the gene
  gene_sorted <- gene_sort(data_frame)
  data_frame$gene <- factor(data_frame$gene, levels=gene_sorted)
  
  # reorder the samples based on hiearchial sort on ordered gene list
  sample_order <- sample_sort(data_frame)
  data_frame$sample <- factor(data_frame$sample, levels=sample_order)
  
  # Subset the data based on the recurrence of mutations at the gene level
  data_frame <- mutation_recurrence_subset(data_frame, recurrence_cutoff)
  
  # Subset the data based on a vector of genes if supplied
  if(!is.null(genes))
  {
    data_frame <- mutation_sample_subset(data_frame, genes)
  }
  
  # Reorder the sample levels in data_frame2 to match the main plot's levels, and then plot the top margin plt
  data_frame2$sample <- factor(data_frame2$sample, levels=sample_order)
  p3 <- plot_mutation_recurrence(data_frame2, coverage_space)
  
  # Plot the Left Bar Chart
  p2 <- plot_bar(data_frame)
  
  # if there are any NA values in the data frame for a gene, give these NA values a gene name so they are plotted properly
  data_frame <- add_gene_to_NA(data_frame)
  
  # Plot the Heatmap
  p1 <- plot_heatmap(data_frame, grid = grid, label_x = label_x, gene_label_size = gene_label_size, file_type = file_type)	
  
  # Align the Plots and return as 1 plot
  pA <- align_y(p2, p1, p3, title)
  
  return(pA)
}