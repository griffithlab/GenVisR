#' build chromosome
#'
#' given a data frame with cytogenetic band locations plot chromosome in ggplot
#' @name ideoView_buildMain
#' @param data_frame a data frame with columns chrom, chromStart, chromEnd,
#' name, gieStain, height_min, height_max, alternate, bandcenter, text_y, arm
#' @param chromosome character string specifying UCSC chromosome to plot one of
#' chr... or all
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @param layers additional ggplot2 layers to plot
#' @return ggplot object
#' @import ggplot2
#' @noRd

ideoView_buildMain <- function(data_frame, chromosome,
                               chr_txt_angle=chr_txt_angle,
                               chr_txt_size=chr_txt_size, layers=NULL)
{
    # define theme layer for ggplot
    theme <- theme(axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   legend.position='right')
    
    # modify ggplot legend defaults
    legend <- scale_fill_manual("Giemsa stain",
                                values=c('gneg'='grey100', 'stalk'='brown3',
                                         'acen'='brown4', 'gpos'='grey0',
                                         'gvar'='grey0', 'gpos1'='#FFFFFF',
                                         'gpos2'='#FCFCFC', 'gpos3'='#F9F9F9',
                                         'gpos4'='#F7F7F7', 'gpos5'='#F4F4F4',
                                         'gpos6'='#F2F2F2', 'gpos7'='#EFEFEF',
                                         'gpos8'='#ECECEC', 'gpos9'='#EAEAEA',
                                         'gpos10'='#E7E7E7', 'gpos11'='#E5E5E5',
                                         'gpos12'='#E2E2E2', 'gpos13'='#E0E0E0',
                                         'gpos14'='#DDDDDD', 'gpos15'='#DADADA',
                                         'gpos16'='#D8D8D8', 'gpos17'='#D5D5D5',
                                         'gpos18'='#D3D3D3', 'gpos19'='#D0D0D0',
                                         'gpos20'='#CECECE', 'gpos21'='#CBCBCB',
                                         'gpos22'='#C8C8C8', 'gpos23'='#C6C6C6',
                                         'gpos24'='#C3C3C3', 'gpos25'='#C1C1C1',
                                         'gpos26'='#BEBEBE', 'gpos27'='#BCBCBC',
                                         'gpos28'='#B9B9B9', 'gpos29'='#B6B6B6',
                                         'gpos30'='#B4B4B4', 'gpos31'='#B1B1B1',
                                         'gpos32'='#AFAFAF', 'gpos33'='#ACACAC',
                                         'gpos34'='#AAAAAA', 'gpos35'='#A7A7A7',
                                         'gpos36'='#A4A4A4', 'gpos37'='#A2A2A2',
                                         'gpos38'='#9F9F9F', 'gpos39'='#9D9D9D',
                                         'gpos40'='#9A9A9A', 'gpos41'='#979797',
                                         'gpos42'='#959595', 'gpos43'='#929292',
                                         'gpos44'='#909090', 'gpos45'='#8D8D8D',
                                         'gpos46'='#8B8B8B', 'gpos47'='#888888',
                                         'gpos48'='#858585', 'gpos49'='#838383',
                                         'gpos50'='#808080', 'gpos51'='#7E7E7E',
                                         'gpos52'='#7B7B7B', 'gpos53'='#797979',
                                         'gpos54'='#767676', 'gpos55'='#737373',
                                         'gpos56'='#717171', 'gpos57'='#6E6E6E',
                                         'gpos58'='#6C6C6C', 'gpos59'='#696969',
                                         'gpos60'='#676767', 'gpos61'='#646464',
                                         'gpos62'='#616161', 'gpos63'='#5F5F5F',
                                         'gpos64'='#5C5C5C', 'gpos65'='#5A5A5A',
                                         'gpos66'='#575757', 'gpos67'='#545454',
                                         'gpos68'='#525252', 'gpos69'='#4F4F4F',
                                         'gpos70'='#4D4D4D', 'gpos71'='#4A4A4A',
                                         'gpos72'='#484848', 'gpos73'='#454545',
                                         'gpos74'='#424242', 'gpos75'='#404040',
                                         'gpos76'='#3D3D3D', 'gpos77'='#3B3B3B',
                                         'gpos78'='#383838', 'gpos79'='#363636',
                                         'gpos80'='#333333', 'gpos81'='#303030',
                                         'gpos82'='#2E2E2E', 'gpos83'='#2B2B2B',
                                         'gpos84'='#292929', 'gpos85'='#262626',
                                         'gpos86'='#242424', 'gpos87'='#212121',
                                         'gpos88'='#1E1E1E', 'gpos89'='#1C1C1C',
                                         'gpos90'='#191919', 'gpos91'='#171717',
                                         'gpos92'='#141414', 'gpos93'='#121212',
                                         'gpos94'='#0F0F0F', 'gpos95'='#0C0C0C',
                                         'gpos96'='#0A0A0A', 'gpos97'='#070707',
                                         'gpos98'='#050505', 'gpos99'='#020202',
                                         'gpos100'='#000000'))

    # define additional layers for ggplot including p arm/q arm text labels,
    # segments connecting label to chromosome
    P_arm_text <- geom_text(data=subset(data_frame,
                                        data_frame$alternate == 'top'),
                            mapping=aes_string(x='band_center',
                                               y='text_y',
                                               label='name'),
                            angle=chr_txt_angle,
                            hjust=0,
                            size=chr_txt_size)
    
    Q_arm_text <- geom_text(data=subset(data_frame,
                                        data_frame$alternate == 'bottom'),
                            mapping=aes_string(x='band_center',
                                               y='text_y',
                                               label='name'),
                            angle=chr_txt_angle,
                            hjust=1,
                            size=chr_txt_size)
    
    text_line_seg_p <- geom_segment(data=subset(data_frame,
                                                data_frame$alternate == 'top'),
                                    mapping=aes_string(x='band_center',
                                                       y='text_y',
                                                       xend='band_center',
                                                       yend='height_max'))
    text_line_seg_q <- geom_segment(data=subset(data_frame,
                                                data_frame$alternate == 'bottom'),
                                    mapping=aes_string(x='band_center',
                                                       y='text_y',
                                                       xend='band_center',
                                                       yend='height_min'))
    # Create Y label
    ylabel <- ylab(chromosome)
    
    # Allow an additional layer from the user
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }

    # define the chromosome main body plot
    chr <- ggplot(data_frame,
                  aes_string(xmin='chromStart', xmax='chromEnd',
                             ymin='height_min', ymax='height_max')) +
        geom_rect(aes_string(fill='gieStain'), colour='black') + ylim(-1.2, 1.2)

    #plot the resulting layers
    chr <- chr + P_arm_text + Q_arm_text + text_line_seg_p + text_line_seg_q +
        theme_bw() + theme + ylabel + legend + layers

    return(chr)
}
