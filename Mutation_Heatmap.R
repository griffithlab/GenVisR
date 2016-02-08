#******************************************************************************************************************************#
################################## R script to create a "mutational heatmap"  ##################################################
#******************************************************************************************************************************#

################################################################################################################################
############################################ preliminary setup/function definitions ############################################
################################################################################################################################

setwd('~/RTools/Mutation_Landscape/')

annotation_header <- function(data_frame)
{
	##################################################################################
	####### function to add column names to file in standard annotation format #######
	##################################################################################
	
	colnames(data_frame) <- c('chromosome_name', 'start', 'stop', 'reference', 'variant', 'type', 'gene', 'transcrip_name', 'transcript_species', 'transcript_source', 'transcript_version', 'strand', 'transcript_status', 'trv_type', 'c_position', 'amino_acid_change', 'ucsc_cons', 'domain', 'all_domains', 'deletion_substructures', 'transcript_error', 'default_gene_name', 'gene_name_source', 'ensembl_gene_id')
	
	return(data_frame)
}

MAF_to_anno <- function(data_frame)
{
	##################################################################################################################
	############## Function to take a MAF file and coerce it into a format recognizable by other functions ###########
	##################################################################################################################
	
	data_frame <- data_frame[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')]
	colnames(data_frame) <- c('sample', 'gene', 'trv_type')
	
	return(data_frame)
}

add_mutation_counts <- function(data_frame)
{
	#############################################################################################################################
	###################### Function to obtain Mutation Frequency counts on a sample level as a data frame #######################
	#############################################################################################################################
	
	require(reshape2)
	
	# Change trv_type calls to either synonymous or non synonymous, for use in the mutation per Mb plot
	data_frame$trv_type <- as.character(data_frame$trv_type)
	data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <- 'Non Synonymous'
	data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <- 'Synonymous'
	data_frame$trv_type <- factor(data_frame$trv_type, levels=c('Synonymous', 'Non Synonymous'))
	
	# Obtain a data frame of mutation counts on the sample level
	mutation_counts <- table(data_frame[,c('sample', 'trv_type')])
	mutation_counts <- as.data.frame(melt(mutation_counts))
	colnames(mutation_counts) <- c('sample', 'trv_type', 'mutation_total')
	
	return(mutation_counts)
}

hiearchial_remove_trv_type <- function(data_frame, file_type)
{
	#############################################################################################################################
	####### Function to hiearchial collapse the data frame on sample/gene leaving only the most deleterious trv_type, ###########
	####### the function expects a dataframe with column names sample, gene, trv_type                                 ###########
	#############################################################################################################################
	
	# reorder the trv_type in terms of deleterious effect and refactor the data frame
	if(toupper(file_type) == toupper('TGI'))
	{
		mutation_order <- c("nonsense", "in_frame_del_stop_gain", "frame_shift_del", "frame_shift_ins", "in_frame_del", "in_frame_ins", "nonstop", "splice_site_del", "splice_site_ins", "splice_site", "missense", "5_prime_flanking_region", "3_prime_flanking_region", "3_prime_untranslated_region", "5_prime_untranslated_region", "rna", "intronic", "silent")
	} else if (toupper(file_type) == toupper('MAF'))
	{
		mutation_order <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Splice_Site", "Missense_Mutation", "5\'Flank", "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron", "IGR", "Silent", "Targeted_Region")
	}
	
	data_frame$trv_type <- factor(data_frame$trv_type, levels=mutation_order)
	
	# sort the data frame so that the duplicated call will remove the proper trv_type
	data_frame <- data_frame[order(data_frame$sample, data_frame$gene, data_frame$trv_type),]
	
	# collapse the data on sample/gene
	data_frame <- data_frame[!duplicated(data_frame[, c("sample", "gene")]), ]
	
	return(data_frame)
}

mutation_recurrence_subset <- function(data_frame, recurrence_cutoff)
{
	##############################################################################################################################
	####### Function to subset a data frame based on a mutation recurrence cutoff at the gene level, expects a data frame  #######
	####### in long format with columns 'gene', 'trv_type'                                                                 #######
	##############################################################################################################################
	
	# convert the dataframe to a table of mutation counts at the gene level
	mutation_recurrence <- table(data_frame[,c('gene', 'trv_type')], useNA="ifany")
	
	# add totals to the data frame (total mutations for each gene and total mutations for each mutation type), remove the last row (total mutations for each mutation)
	mutation_recurrence_total <- as.data.frame.matrix(addmargins(mutation_recurrence, FUN = list(Total = sum), quiet = TRUE))
	mutation_recurrence_total <- mutation_recurrence_total[-nrow(mutation_recurrence_total),]
	
	# If recurrence cutoff specified exceeds upper limit such that no plot useful would be generated, reset recurrence cutoff
	if(max(mutation_recurrence_total$Total) < recurrence_cutoff)
	{
		message <- paste0('The recurrence cutoff specified exceeds the recurrence seen in the data, resetting this value to equal max recurrence:', max(mutation_recurrence_total$Total))
		warning(message)
		recurrence_cutoff <- max(mutation_recurrence_total$Total)
	}
	
	# obtain a vector of gene names where the genes have mutation recurrence greater than or equal to the recurrence cutoff
	gene_above_recurrence <- sort(row.names(subset(mutation_recurrence_total, mutation_recurrence_total$Total >= recurrence_cutoff)))
	
	# add NA to the end of 'gene_above_recurrence' vector, allowing for all samples having NA as a gene name to be retained in the subset below
	gene_above_recurrence <- c(gene_above_recurrence, NA)
	
	# subset the original data frame based on the following: keep gene if it is in the gene vector in "mutation_recurrence_subset"
	subset_data_frame <- data_frame[(data_frame$gene %in% gene_above_recurrence), ]
	
	return(subset_data_frame)
}

gene_sort <- function(data_frame)
{
	
	#########################################################################################################################################
	####### Create a table of genes and their mutation status, sum the mutation totals for each gene and sort giving a vector of gene #######
	####### names sorted by mutation frequency, note that the function expects a data frame in long format with columns 'gene', 'trv_type' ##
	#########################################################################################################################################
	
	gene_mutation_table <- table(data_frame[,c('gene', 'trv_type')])
	gene_order <- names(sort(rowSums(gene_mutation_table)))
	return(gene_order)
}

convert_to_boolean <- function(x)
{
	
	###########################################################################################################
	################ Function to convert counts to a boolean based on if there or not #########################
	###########################################################################################################
	
	# if greater than 1 return 1
	i = which(x > 1)
	x[i] = 1
	return(x)
}

sample_sort <- function(data_frame)
{
	
	#########################################################################################################################################
	######## Function to perform a hiearchial sort on samples based on the presence of mutations in an ordered list of genes, expects #######
	######## a dataframe in long format with columns 'sample', 'gene', 'trv_type'                                                     #######
	#########################################################################################################################################
	
	require(reshape2)
	
	# recast the data going from long format to wide format, values in this data are counts of a mutation call
	wide_data <- dcast(data_frame, sample ~ gene, fun.aggregate = length, value.var="trv_type")
	
	# apply a boolean function to convert the data frame values to 1's and 0's
	wide_data[,-1] <- t(apply(wide_data[,-1], 1, convert_to_boolean))
	wide_boolean <- wide_data
	
	# move the first column (sample names) to row names
	wide_boolean_format <- data.frame(wide_boolean[,-1], row.names=wide_boolean[,1])
	
	# reverse the columns so that genes with highest mutation's are listed first (assuming gene_sort has been run on the data frame)
	test <- wide_boolean_format[,rev(1:ncol(wide_boolean_format))]

	# if there are any NA values present in a sample at the gene put that NA gene last so samples with NA are plotted last
	if(any(colnames(test) == 'NA.'))
	{
		# Find which column has the NA header
		NA_index <- which(colnames(test) == 'NA.')
		
		# Append NA column to end of data frame
		NA_gene <- test[,NA_index]
		test <- test[,-NA_index]
		test <- cbind(test, NA_gene)
	}
	
	# hiearchial sort on all column's (i.e. genes) such that samples are rearranged if there is a mutation in that gene
	sample_order <- rownames(test[do.call(order, as.list(-test)),])
	
	return(sample_order)
}

add_gene_to_NA <- function(data_frame)
{
	########################################################################################################
	######### Function to replace na values in a dataframe with a gene name, i.e. when the samples #########
	######### with NA values for gene/trv_type are read in add a gene to these NA values so that   #########
	######### they are plotted on the x axis instead of an "NA" gene appearing                     #########
	########################################################################################################
	
	# Get The gene with the most Mutations and add the NA samples to that gene (ensures that the NA's are added in as gene with most mutations will always be plotted)
	top_gene <- tail(levels(data_frame$gene), 1)
	data_frame$gene <- replace(data_frame$gene, is.na(data_frame$gene), top_gene)
	
	return(data_frame)
}

plot_heatmap <- function(data_frame, grid, label_x, gene_label_size, file_type)
{
	
	#############################################################################################################
	################ Plotting call to return a mutation heatmap, see comments/notes below for specifics #########
	#############################################################################################################
	
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ NOTES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
	# 1. scale_x_discret(drop=FALSE), added to ggplot2 call to ensure samples remove by 'mutation_recurrence_subset' are still plotted as empty tiles #
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
	
	####################################### define layers for ggplot call ########################################
	
	# Grid options
	vertical_grid <- geom_vline(xintercept = seq(.5, nlevels(data_frame$sample), by=1), linetype='solid', colour='grey80', size=.01)
	
	if(length(unique(data_frame$gene)) == 1)
	{
		horizontal_grid <- geom_hline()
	} else {
		horizontal_grid <- geom_hline(yintercept = seq(1.5, length(unique(data_frame$gene)), by=1), linetype='solid', colour='grey80', size=.01)
	}
	
	if(toupper(file_type) == toupper('TGI'))
	{
		
	# Declare a color palette
	palette <- c('#A80100', '#CF5A59', '#A80079', '#CF59AE', '#4f00A8', '#9159CF', '#000000', '#006666', '#00A8A8', '#79F2F2', '#59CF74', '#002AA8', '#5977CF', '#F37812', '#F2B079', '#888811', '#FDF31C', '#8C8C8C')
	
	# Create Legend labels
	breaks <- c("nonsense", "in_frame_del_stop_gain", "frame_shift_del", "frame_shift_ins", "in_frame_del", "in_frame_ins", "nonstop", "splice_site_del", "splice_site_ins", "splice_site", "missense", "5_prime_flanking_region", "3_prime_flanking_region", "3_prime_untranslated_region", "5_prime_untranslated_region", "rna", "intronic", "silent")
	labels <- c("Nonsense", "In Frame Deletion, Stop Gain", "Frame Shift Deletion", "Frame Shift Insertion", "In Frame Deletion", "In Frame Insertion", "Nonstop", "Splice Site Deletion", "Splice Site Insertion", "Splice Site", "Missense", "5' Flank", "3' Flank", "3' UTR", "5' UTR", "RNA", "Intronic", "Silent")
	
	} else if(toupper(file_type) == toupper('MAF'))
	{
	
	# Declare a color palette
	palette <- c('#A80100', '#CF5A59', '#A80079', '#CF59AE', '#4f00A8', '#9159CF', '#000000', '#59CF74', '#00A8A8', '#79F2F2', '#006666', '#002AA8', '#5977CF', '#F37812', '#F2B079', '#888811', '#FDF31C')
	
	'#006666'
	'#59CF74'	
	
	# Create Legend Labels
	breaks <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Splice_Site", "Missense_Mutation", "5\'Flank", "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron", "IGR", "Silent", "Targeted_Region")
	labels <- c("Nonsense", "Frame Shift Insertion", "Frame Shift Deletion", "In Frame Insertion", "In Frame Deletion", "Nonstop", "Splice Site", "Missense", "5' Flank", "3' Flank", "5' UTR", "3' UTR", "RNA", "Intron", "Intergenic Region", "Silent", "Targeted Region")
	}
	
	# Create Legend
	legend <- scale_fill_manual(name="Mutation Type", values=palette, breaks=breaks, labels=labels, drop=FALSE)
	
	# X Label
	x_label <- xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
	
	# Plot Title
	title <- ggtitle(title)

	# Theme, Boolean, if specified to plot x labels, define theme such that labels are plotted
	if(label_x == TRUE)
	{
		theme <-  theme(axis.ticks=element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill='white', colour='white'), axis.text.x=element_text(angle=50, hjust=1), axis.text.y=element_text(size=gene_label_size, colour='black', face='italic'), axis.title.y=element_blank(), axis.title.x=element_text(size=10), legend.title=element_text(size=14), plot.title=element_blank())
	} else {
		theme <-  theme(axis.ticks=element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill='white', colour='white'), axis.text.x=element_blank(), axis.text.y=element_text(size=gene_label_size, colour='black', face='italic'), axis.title.y=element_blank(), axis.title.x=element_text(size=20), legend.title=element_text(size=14), plot.title=element_blank(), panel.border=element_rect(colour='grey80', fill=NA, size=.1), legend.position=("right"))
	}
	
	# ggplot call
	if(grid == TRUE)
	{
		p1 <- ggplot(data_frame, aes(sample, gene)) + geom_tile(aes(fill=trv_type), position="identity") + theme + legend + ggtitle(title) + x_label + vertical_grid + horizontal_grid + scale_x_discrete(drop=FALSE) + title
	} else {
		p1 <- ggplot(data_frame, aes(sample, gene)) + geom_tile(aes(fill=trv_type), position="identity") + theme + legend + ggtitle(title) + x_label + scale_x_discrete(drop=FALSE) + title
	}
	
	return(p1)
}

plot_bar <- function(data_frame)
{
	#####################################################################################################################
	####### Function to create the left margin bar plot, function plots the percentage of samples with a mutation #######
	#####################################################################################################################
	
	# Require the scales package for use of ..count..
	require(scales)
	
	# Convert all silent mutations to Synonymous, and all else to non-synonymous
	data_frame$trv_type <- as.character(data_frame$trv_type)
	data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <- 'Non Synonymous'
	data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <- 'Synonymous'
	data_frame$trv_type <- factor(data_frame$trv_type, levels=c('Synonymous', 'Non Synonymous'))
	
	# Define the number of samples for the Percentage calculation (Note, to pass a variable outside of aes into aes it needs to be defined again)
	total_number_sample <- nlevels(data_frame$sample)
	
	# Define Theme and various other layers to be passed to ggplot
	theme <- theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position=('none'))
	y_limits <- ylim(100, 0)
	y_label <- ylab('% Samples With Mutation')	
	legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c('Synonymous', 'Non Synonymous'), drop=FALSE)
	
	# Plotting call
	p1 <- ggplot(na.omit(data_frame), aes(x=gene, total_number_sample=total_number_sample, y=(..count..)/total_number_sample * 100, fill=trv_type), environment = environment()) + geom_bar(position='stack', alpha=.75) + coord_flip() + theme + y_label + scale_y_reverse() + legend	
	
	return(p1)
}

plot_mutation_recurrence <- function(data_frame, coverage_space)
{
	#############################################################################################################################
	#################################### Function to plot the top margin barplot ################################################
	#############################################################################################################################
	
	# Add in mutations per MB calculation
	data_frame$mutation_per_MB <- data_frame$mutation_total/coverage_space * 1000000
	
	#print(data_frame)
	
	# Alter GGplot2 Theme 
	theme <- theme(panel.border =  element_blank(), axis.line =  element_line(), panel.background=element_rect(fill='white'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_text(size=14))
	
	# Add Legend
	legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c("Synonymous", "Non Synonymous"), drop=FALSE)
	
	# add y label
	y_label <- ylab('Mutations per MB')
	
	# ggplot2 call
	p1 <- ggplot(data_frame, aes(x=sample, y=mutation_per_MB, fill=trv_type)) + geom_bar(stat='identity', alpha=.75) + theme + y_label + legend
	
	return(p1)
}

align_y <- function(p2, p1, p3, title)
{
	
	#############################################################################################################
	############## Function to take three ggplots and align the plotting space on the y_axis and x_axis #########
	#############################################################################################################
	
	require(gridExtra)
	require(gtable)
	
	# define the ggplot's as grobs and create a blank plot
	gA <- ggplotGrob(p2)
	gB <- ggplotGrob(p1)
	gC <- ggplotGrob(p3)
	blankPanel<-grid.rect(gp=gpar(col="white"))
	
	# Adjust the grob widths so p1 and p3 plots line up
	maxwidth = grid::unit.pmax(gB$widths[2:5,], gC$widths[2:5,])
	gC$widths[2:5] <- as.list(maxwidth)
	gB$widths[2:5] <- as.list(maxwidth)
	
	# Adjust the grob heights so p1, and p2 plots line up
	maxheight = grid::unit.pmax(gA$heights[2:5,], gB$heights[2:5,])
	gA$heights[2:5] <- as.list(maxheight)
	gB$heights[2:5] <- as.list(maxheight)
	
	# plot the grobs with grid.arrange
	if(is.character(title))
	{
		p1 <- grid.arrange(blankPanel, gC, gA, gB, ncol=2, nrow=2, widths=c(1,4), heights=c(1,4), main=textGrob(title, gp=gpar(fontsize=20)))
	} else {
		p1 <- grid.arrange(blankPanel, gC, gA, gB, ncol=2, nrow=2, widths=c(1,4), heights=c(1,4))
	}
	

	return(p1)
}

mutation_heatmap <- function(data_frame, recurrence_cutoff = 0, grid = TRUE, label_x = FALSE, title ='NULL', gene_label_size=8, coverage_space=30000000, file_type='TGI')
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
