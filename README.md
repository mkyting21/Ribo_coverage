# Ribo_coverage
Visualize Ribo-seq and RNA-seq coverage over gene annotations.
Outputs a PDF figure over Gene of interest.

## Example usage:
```
Rscript Ribo_coverage.R \
--gene_id AT1G01010 \
--annotation test.gtf \
--rna test_rna.bedgraph \
--ribo test_ribo.bedgraph \
--psite test_psite.wig
```
