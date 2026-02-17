#combined coverage plots for rna-seq and ribo-seq 

# Command line options ---------------------------------------------------
library("optparse")
option_list <- list(
  make_option( c("-g", "--genes"), type = "character", default = NULL, metavar = "character",
    help = "File with all gene_ids of interest. Each id should be on a new line"),
  make_option( c("-a", "--annotation"), type = "character", default = NULL, metavar = "character",
    help = "gtf annotation file"),
  make_option( c("--rna"), type = "character", default = NULL, metavar = "character",
    help = "rnaseq coverage file: bigWig (recommended) or bedgraph format"),
  make_option( c("--ribo"), type = "character", default = NULL, metavar = "character",
    help = "riboseq coverage file: bigWig (recommended) or bedgraph format"),
  make_option( c("--psite_F"), type = "character", default = NULL, metavar = "character",
    help = "Psite file: bigwig (recommended) or wig format")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
#stop the script if not all options are used
if( length(opt) < 6 ){
  print("ERROR: all arguments must be selected. Refer to -h ")
  stop()
}

#If running interactively, run this and enter file paths manually 
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  opt <- list()
  opt$annotation <- "data/test_data/test.gtf"
  opt$rna <- "data/test_data/test_rna.bedgraph"
  opt$ribo <- "data/test_data/test_ribo.bedgraph"
  opt$psite <- "data/test_data/test_psite.wig"
  opt$genes <- "data/test_data/test_genes.txt"
}

cat( "My variables:", "\n", 
 "annotation = ", opt$annotation, "\n",
 "rna = ", opt$rna, "\n",
 "ribo = ", opt$ribo, "\n",
 "psite = ", opt$psite, "\n",
 "genes = ", opt$genes, "\n" 
)

# Dependencies -----------------------------------------------------------
library("tools") 
library("rtracklayer") #importing gtf and bedgraphs
library("tidyverse")
library("ggtranscript") #geom_gene function for visualizing transcript annotations
library("cowplot") #stick together multiple plots


# Main -------------------------------------------------------------------

#detect file format of a genome coverage file. This function is used with Rtracklayer::import()
identify_file_format <- function(file_name){
  file_extension <- file_name %>% file_ext() %>% str_to_lower()
  if( file_extension %in% c("bedgraph", "bg") ){ 
    file_format <- "bedGraph"
  }else if( file_extension %in% c("bigwig", "bw") ){
    file_format <- "BigWig"
  }else if( file_extension %in% c("wig") ){
    file_format <- "wig"
  }else( 
    paste0("ERROR: invalid filetype for ", file_name)
  )  
  return(file_format)
}

#get gene info
my_gtf <- import(opt$annotation) 
my_gene_list <- readLines(opt$genes)

create_plots <- function(my_gene_id){
  #extract gene information
  gene_gtf <- my_gtf[ mcols(my_gtf)$gene_id %in% my_gene_id ] 
  gene_name <- unique( gene_gtf$gene_name )
  if( length(gene_gtf) <1){ 
    print("gene_name not found in annotation file")
    stop()
  }
  gene_info <- gene_gtf[ gene_gtf$type == "gene" ] %>% as.data.frame()
  transcripts_gtf <- gene_gtf[gene_gtf$type != "gene"] %>% as.data.frame()
  x_axis_limits <- c( gene_info$start, gene_info$end) #all plots must have same coordinate cartesians to align

  #import coverage info that overlaps with the gene gtf
  coverage_rna <- import(opt$rna, format = identify_format(opt$rna), which = gene_gtf) 
  coverage_ribo  <- import(opt$ribo, format = identify_format(opt$ribo), which = gene_gtf)
  coverage_psite <- import(opt$psite, format = identify_format(opt$psite), which = gene_gtf)

  #If the region of interest has 0 coverage, fill it in with default values
  check_coverage <- function(coverage_file){    
    if( length(coverage_file) == 0){      
      coverage_file <- gene_info 
      coverage_file$score <- 0
    }
    return( as.data.frame(coverage_file) ) 
  }
  coverage_rna <- check_coverage(coverage_rna)
  coverage_ribo <- check_coverage(coverage_ribo)
  coverage_psite <- check_coverage(coverage_psite)

  #add a data type, used for ggplot
  coverage_ribo$type <- "ribo"
  coverage_psite$type <- "ribo"
  coverage_rna$type <- "rna"

  #plot transcript features
  features_plot <- ggplot( transcripts_gtf, aes(xstart = start, xend = end, y = transcript_id) )  + theme_bw()+
    geom_intron( aes(strand = strand), arrow.min.intron.length = 100000) +
    geom_range(data = filter(transcripts_gtf, type == "exon"), fill = "white", height = 0.25) +
    labs(y = "") + 
    coord_cartesian(xlim = x_axis_limits) +
    theme(legend.position = "none",
          plot.margin = unit( c(0,0,0,0), "mm") )
  
  #plot rna coverage 
  rna_plot <- ggplot(coverage_rna) + theme_bw() +
    geom_rect( aes( xmin = start, xmax = end, ymin = 0, ymax = score), 
      fill = "#756bb1",
      color = NA) +
    labs(y = "coverage", x = "") +
    facet_wrap(~type, strip.position = "right") +
    coord_cartesian(xlim = x_axis_limits) +
    theme(legend.position = "none",
          plot.margin = unit( c(0,0,0,0), "mm"), 
          axis.text.x = element_blank()
        )
  
  #plot riboseq coverage
  ribo_plot <- ggplot(coverage_ribo, aes(x = start, y = score) ) + theme_bw() +
    geom_rect( aes( xmin = start, xmax = end, ymin = 0, ymax = score), 
      fill = "#31a354",
      color = NA) +
    labs(y = "coverage", x = "") + 
    geom_col( data = coverage_psite, aes(x = start, y = score), color = "red" ) +
    facet_wrap(~type, strip.position = "right") +
    coord_cartesian(xlim = x_axis_limits) +  
    theme(legend.position = "none",
          plot.margin = unit( c(0,0,0,0), "mm"),
          axis.text.x = element_blank()
        )

  #combine plots and save
  final_plot <- plot_grid( rna_plot, ribo_plot, features_plot, nrow = 3, rel_heights = c(2, 2, 1), align = "v", axis = "lr")  
  pdf( file = paste0( my_gene_id, "_", gene_name, "-coverage_plot.pdf"), width = 5, height = 5 )
  print(final_plot)
  dev.off()

  png( file = paste0( my_gene_id, "_", gene_name, "-coverage_plot.png"), res = 300, width = 1200, height = 1200)
  print(final_plot)
  dev.off()

}
map(my_gene_list, create_plots)
