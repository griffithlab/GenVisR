#' Construct Lolliplot from data
#'
#' Construct a Lolliplot from object of class data frame giving observed
#' mutations and an ensembl transcript id
#' @name lolliplot
#' @param x object of class data frame containing column names "transcript_name"
#', "gene", and "amino_acid_change" with rows denoting mutations for top track
#' @param y object of class data frame containing columns "transcript_name", and
#' "amino_acid_change" with rows denoting mutations for bottom track (optional)
#' @param z Object of class data frame containing columns "description",
#' "start", "stop" for manual denotation of regions of interest
#' @param fillCol character string giving the name of the column to shade
#' variants on, required
#' @param labelCol character string specifying column containing text
#' information to be plotted, (optional) defaults to NULL
#' @param plot_text_angle integer specifying angle of text to be plotted if
#' labelCol is specified
#' @param plot_text_size integer specifying the size of the text plotted if
#' labelCol is specified
#' @param point_size integer specifying the size of points plotted
#' @param gene_colour color specifying the background fill of the plotted gene
#' @param obsA.rep.fact repulsive factor for plotted mutations observed track
#' @param obsA.rep.dist.lmt repulsive distance limit for plotted mutations
#' observed track
#' @param obsA.attr.fact attraction factor for plotted mutations observed track
#' @param obsA.adj.max maximum position change for each iteration observed track
#' @param obsA.adj.lmt position adjustment limit which simulation stops observed
#' track
#' @param obsA.iter.max maximum iterations beyond which to stop the simulation
#' observed track
#' @param obsB.rep.fact repulsive factor for plotted mutations cosmic track
#' @param obsB.rep.dist.lmt repulsive distance limit for plotted mutations
#' cosmic track
#' @param obsB.attr.fact attraction factor for plotted mutations cosmic track
#' @param obsB.adj.max maximum position change for each iteration cosmic track
#' @param obsB.adj.lmt position adjustment limit which simulation stops cosmic
#' track
#' @param obsB.iter.max maximum iterations beyond which to stop the simulation
#' cosmic track
#' @param plot_sidechain boolean value whether to plot the amino acid sidechains
#' instead of protein domains
#' @param species Character string specifying the species corresponding to the
#' given ensembl transcript id
#' @param max_VerticalLolli Integer specifying the upper limit of lolli points to allow to be stacked vertically
#' @param layers additional ggplot2 layers to plot
#' @param paletteA Character vector specifying colours for gene features
#' @param paletteB Character vector specifying colours for lolli features
#' @examples
#' # Create input data
#' data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
#' data <- as.data.frame(cbind(data, 'ENST00000269305'))
#' colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')
#'
#' # Call lolliplot
#' lolliplot(data)
#' @return object of class ggplot2
#' @export

lolliplot <- function(x, y=NULL, z=NULL, fillCol=NULL, labelCol=NULL,
                      plot_text_angle=45, plot_text_size=5, point_size=4,
                      gene_colour='#999999', obsA.rep.fact=5000,
                      obsA.rep.dist.lmt=500, obsA.attr.fact=.1, obsA.adj.max=.1,
                      obsA.adj.lmt=.5, obsA.iter.max=50000, obsB.rep.fact=5000,
                      obsB.rep.dist.lmt=500, obsB.attr.fact=.1, obsB.adj.max=.1,
                      obsB.adj.lmt=.5, obsB.iter.max=50000,
                      plot_sidechain=FALSE, species="hsapiens",
                      max_VerticalLolli=NULL, layers=NULL, paletteA=NULL, paletteB=NULL)
{
    # Perform quality check
    input <- lolliplot_qual(x, y, z)
    x <- input[[1]]
    y <- input[[2]]
    z <- input[[3]]

    # extract transcript id and subset data y on that id if it exists
    transcriptID <- as.character(x$transcript_name[1])
    if(!is.null(y))
    {
        y <- y[y$transcript_name == transcriptID,]
    }

    # extract HUGO gene name
    gene <- as.character(x$gene[1])

    # Obtain length of protein
    result <- lolliplot_transcriptID2codingSeq(transcriptID, species=species)
    codingSeq <- result$coding
    cdsLen <- result$cds_length
        
    
    # Get the sequence length in AA, perform quality checks along the way
    residueSeq <- lolliplot_DNAconv(codingSeq, to="residue")
    residueSeq <- residueSeq[-which(residueSeq %in% c("OPAL",
                                                      "OCHRE",
                                                      "AMBER"))]
    proteinLength <- length(residueSeq)    
    
    # obtain amino acid sequence and format if it is requested to plot the
    # sidechain
    if(plot_sidechain==TRUE)
    {
        AAsequence <- lolliplot_DNAconv(codingSeq, to="sidechain")
    } else {
        AAsequence <- NULL
    }
    
    # if z is specified plot that instead of fetching the domain information
    if(!is.null(z))
    {
        geneData <- lolliplot_constructGene(gene, z, proteinLength)
    } else {
        # extract protien domain data
        protien_domain <- lolliplot_fetchDomain(transcriptID, species=species)
        
        # construct gene from data collected
        geneData <- lolliplot_constructGene(gene, protien_domain, proteinLength)
    }
 
    # construct data frame of observed mutations for top track
    observed_mutation <- lolliplot_mutationObs(x, 'top', fillCol, labelCol,
                                               obsA.rep.fact, obsA.rep.dist.lmt,
                                               obsA.attr.fact, obsA.adj.max,
                                               obsA.adj.lmt, obsA.iter.max)
    observed_mutation <- lolliplot_reduceLolli(observed_mutation,
                                               max_VerticalLolli)

    # construct data frame of observed mutations for bottom track
    if(!is.null(y))
    {
        observed_mutation2 <- lolliplot_mutationObs(y, 'bottom', fillCol,
                                                    labelCol, obsB.rep.fact,
                                                    obsB.rep.dist.lmt,
                                                    obsB.attr.fact,
                                                    obsB.adj.max, obsB.adj.lmt,
                                                    obsB.iter.max)
        observed_mutation2 <- lolliplot_reduceLolli(observed_mutation2,
                                                   max_VerticalLolli)
    } else {
        observed_mutation2 <- NULL
    }


    
    # construct the lolliplot
    plot <- lolliplot_buildMain(geneData, length, observed_mutation,
                                observed_mutation2,fillCol, labelCol,
                                plot_text_angle, plot_text_size,point_size,
                                gene_colour, AAsequence,
                                plot_sidechain=plot_sidechain, layers=layers,
                                paletteA=paletteA, paletteB=paletteB)

    return(plot)
}
