#!/usr/bin/env Rscript
#combined coverage plots for rna-seq and ribo-seq 

#Command line options
library("optparse")
option_list <- list(
  make_option( c("-g", "--gene_id"), type = "character", default = NULL, metavar = "character",
    help = "gene id"),
  make_option( c("-a", "--annotation"), type = "character", default = NULL, metavar = "character",
    help = "gtf annotation file"),
  make_option( c("--rna"), type = "character", default = NULL, metavar = "character",
    help = "rnaseq coverage file in bedgraph format"),
  make_option( c("--ribo"), type = "character", default = NULL, metavar = "character",
    help = "riboseq coverage file in bedgraph format"),
  make_option( c("--psite"), type = "character", default = NULL, metavar = "character",
    help = "Psite file in wiggle format")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
#stop the script if not all options are used
if( length(opt) < 6 ){
  print("ERROR: some arguments are blank. Refer to -h ")
  stop()
}

#If running interactively, run this and enter file paths manually 
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  opt <- list()
  opt$annotation <- "data/Arabidopsis_thaliana.TAIR10.45.gtf"
  opt$rna <- "data/RNA-A1-1_Nu_F.bedgraph"
  opt$ribo <- "data/RFP-A1-1B.bedgraph"
  opt$psite <- "data/RFP-A1-1B_psite_bedgraph_fw.wig"
  opt$gene_id <- "AT2G46830"
}

#dependencies
library("rtracklayer") #importing gtf and bedgraphs
library("tidyverse")
library("ggtranscript") #geom_gene function for visualizing transcript annotations
library("cowplot") #stick together multiple plots


cat( "My variables:", "\n", 
 "annotation = ", opt$annotation, "\n",
 "rna = ", opt$rna, "\n",
 "ribo = ", opt$ribo, "\n",
 "psite = ", opt$psite, "\n",
 "gene_id = ", opt$gene_id, "\n" 
)

#get gene info
my_gtf <- import(opt$annotation) 
gene_gtf <- my_gtf[ mcols(my_gtf)$gene_id %in% opt$gene_id ] 
gene_name <- unique( gene_gtf$gene_name )
if( length(gene_gtf) <1){ 
  print("gene_name not found in annotation file")
  stop()
}
gene_info <- gene_gtf[ gene_gtf$type == "gene" ] %>% as.data.frame()
transcripts_gtf <- gene_gtf[gene_gtf$type != "gene"] %>% as.data.frame()
x_axis_limits <- c( gene_info$start, gene_info$end) #all plots must have same coordinate cartesians to align

#import coverage info that overlaps with the gene gtf
coverage_rna <- import(opt$rna, format = "bedGraph", which = gene_gtf) 
coverage_ribo  <- import(opt$ribo, format = "bedGraph", which = gene_gtf)
coverage_psite <- import(opt$psite, format = "wig", which = gene_gtf)
#add a data type, used for ggplot
coverage_ribo$type <- "ribo"
coverage_psite$type <- "ribo"
coverage_rna$type <- "rna"

#convert to dataframe
coverage_ribo <- as.data.frame(coverage_ribo)
coverage_rna <- as.data.frame(coverage_rna)  
coverage_psite <- as.data.frame(coverage_psite)

#features plot
features_plot <- ggplot( transcripts_gtf, aes(xstart = start, xend = end, y = transcript_id) )  + theme_bw()+
  geom_intron( aes(strand = strand), arrow.min.intron.length = 100000) +
  geom_range(data = filter(transcripts_gtf, type == "exon"), fill = "white", height = 0.25) +
  labs(y = "") + 
  coord_cartesian(xlim = x_axis_limits) +
  theme(legend.position = "none",
        plot.margin = unit( c(0,0,0,0), "mm") )

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

#combine plots
final_plot <- plot_grid( rna_plot, ribo_plot, features_plot, nrow = 3, rel_heights = c(2, 2, 1), align = "v", axis = "lr")
pdf( file = paste0( opt$gene_id, "_", gene_name, "-coverage_plot.pdf"), width = 5, height = 5 )
print(final_plot)
dev.off()
