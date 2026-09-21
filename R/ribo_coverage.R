#combined coverage plots for rna-seq and ribo-seq 

# Command line options ---------------------------------------------------
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = NULL, metavar = "character",
    help = "File with all gene_ids of interest. Each id should be on a new line"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL, metavar = "character",
    help = "gtf annotation file"),
  make_option("--rna", type = "character", default = NULL, metavar = "character",
    help = "rnaseq coverage file: bigWig (recommended) or bedgraph format"),
  make_option("--ribo", type = "character", default = NULL, metavar = "character",
    help = "riboseq coverage file: bigWig (recommended) or bedgraph format"),
  make_option("--psite", type = "character", default = NULL, metavar = "character",
    help = "Psite file: bigwig (recommended) or wig format")
)

required <- c("genes", "annotation", "rna", "ribo", "psite")

if (interactive()) {
  # When running interactively, set file paths manually
  opt <- list(
    annotation = "data/test_data/test.gtf",
    rna        = "data/test_data/test_rna.bedgraph",
    ribo       = "data/test_data/test_ribo.bedgraph",
    psite      = "data/test_data/test_psite.wig",
    genes      = "data/test_data/test_genes.txt"
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  missing_opts <- required[vapply(opt[required], is.null, logical(1))]
  if (length(missing_opts) > 0) {
    stop("All arguments must be supplied. Missing: ",
         paste0("--", missing_opts, collapse = ", "), ". Refer to -h", call. = FALSE)
  }
}

not_found <- required[!file.exists(unlist(opt[required]))]
if (length(not_found) > 0) {
  stop("File not found for: ", paste(not_found, collapse = ", "), call. = FALSE)
}

cat("My variables:\n", paste0(" ", names(opt[required]), " = ", unlist(opt[required]), "\n"), sep = "")

# Main Dependencies ------------------------------------------------------
suppressPackageStartupMessages({
  library("tools")
  library("rtracklayer")    # importing gtf, bedgraph, bigwig, wig
  library("GenomicRanges")
  library("ggplot2")
  library("dplyr")
  library("purrr")
  library("ggtranscript")   # geom_range / geom_intron for transcript annotations
  library("cowplot")        # stick together multiple plots
})

# Helper functions -------------------------------------------------------

# Detect file format of a coverage file, for use with rtracklayer::import()
identify_format <- function(file_name) {
  switch(tolower(file_ext(file_name)),
    bedgraph = , bg = "bedGraph",
    bigwig   = , bw = "BigWig",
    wig      = "wig",
    stop("Invalid file type for ", file_name, call. = FALSE)
  )
}

# Returns a function(region) -> GRanges of coverage overlapping `region`.
# BigWig files are indexed, so they are queried per region.
# bedGraph/wig files cannot be queried efficiently, so they are read ONCE
# here and subset in memory (instead of re-reading the whole file per gene).
make_coverage_reader <- function(path) {
  fmt <- identify_format(path)
  if (fmt == "BigWig") {
    function(region) import(path, format = fmt, which = region)
  } else {
    cov <- import(path, format = fmt)
    function(region) subsetByOverlaps(cov, region)
  }
}

# Convert coverage to a data.frame; if the region has no coverage,
# fill it with a single zero-score row spanning the gene
coverage_df <- function(gr, type, region) {
  if (length(gr) == 0) {
    return(data.frame(start = start(region), end = end(region), score = 0, type = type))
  }
  df <- as.data.frame(gr)
  df$type <- type
  df
}

# One coverage panel (shared by RNA and Ribo plots)
coverage_panel <- function(df, fill, x_limits) {
  ggplot(df) + theme_bw() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = score),
              fill = fill, color = NA) +
    labs(y = "coverage", x = "") +
    facet_wrap(~type, strip.position = "right") +
    coord_cartesian(xlim = x_limits) +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.text.x = element_blank())
}

# Load data once ---------------------------------------------------------
# Only keep the feature types that are actually plotted
my_gtf <- import(opt$annotation, feature.type = c("gene", "transcript", "exon", "CDS"))
gtf_by_gene <- split(my_gtf, my_gtf$gene_id)   # split once instead of filtering per gene

my_gene_list <- readLines(opt$genes) %>% trimws() %>% unique()
my_gene_list <- my_gene_list[nzchar(my_gene_list)]

read_rna   <- make_coverage_reader(opt$rna)
read_ribo  <- make_coverage_reader(opt$ribo)
read_psite <- make_coverage_reader(opt$psite)

# Main function ----------------------------------------------------------
create_plots <- function(my_gene_id) {
  if (!my_gene_id %in% names(gtf_by_gene)) {
    stop("gene_id not found in annotation file", call. = FALSE)
  }
  gene_gtf <- gtf_by_gene[[my_gene_id]]

  # Region spanning the whole gene (used for coverage queries and x-axis)
  region <- range(gene_gtf, ignore.strand = TRUE)
  x_axis_limits <- c(start(region), end(region))

  coverage_rna   <- coverage_df(read_rna(region),   "rna",  region)
  coverage_ribo  <- coverage_df(read_ribo(region),  "ribo", region)
  coverage_psite <- coverage_df(read_psite(region), "ribo", region)  # same facet as ribo

  # Gene label for file names (fall back to the id if there is no gene_name)
  gene_name <- unique(na.omit(as.character(mcols(gene_gtf)$gene_name)))
  file_stub <- paste(unique(c(my_gene_id, gene_name[seq_len(min(1, length(gene_name)))])),
                     collapse = "_")

  # Transcript features
  transcripts_gtf <- as.data.frame(gene_gtf[gene_gtf$type != "gene"])

  features_plot <- ggplot(transcripts_gtf, aes(xstart = start, xend = end, y = transcript_id)) +
    theme_bw() +
    geom_intron(data = filter(transcripts_gtf, type == "transcript"),
                aes(strand = strand), arrow.min.intron.length = 100000) +
    geom_range(data = filter(transcripts_gtf, type == "exon"), fill = "white", height = 0.15) +
    geom_range(data = filter(transcripts_gtf, type == "CDS"),  fill = "grey",  height = 0.25) +
    labs(y = "") +
    coord_cartesian(xlim = x_axis_limits) +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "mm"))

  # Coverage panels
  rna_plot  <- coverage_panel(coverage_rna, "#756bb1", x_axis_limits)
  ribo_plot <- coverage_panel(coverage_ribo, "#31a354", x_axis_limits) +
    geom_col(data = coverage_psite, aes(x = start, y = score),
             fill = "red", color = NA, width = 1)

  # Combine and save
  final_plot <- plot_grid(rna_plot, ribo_plot, features_plot,
                          nrow = 3, rel_heights = c(2, 2, 1), align = "v", axis = "lr")

  ggsave(paste0(file_stub, "-coverage_plot.pdf"), final_plot, width = 5, height = 5)
  ggsave(paste0(file_stub, "-coverage_plot.png"), final_plot, width = 4, height = 4, dpi = 300)

  invisible(NULL)
}

walk(my_gene_list, function(g) {
  tryCatch(create_plots(g),
           error = function(e) message("Skipping ", g, ": ", conditionMessage(e)))
})