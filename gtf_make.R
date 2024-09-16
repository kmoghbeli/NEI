## Utility file to make GTF file from a Genbank record
library(tidyverse)
library(genbankr)
library(rentrez)

## Create GTF file to use with this tutorial: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
# KOS Genome: https://www.ncbi.nlm.nih.gov/nuccore/JQ673480
# RE Genome: https://www.ncbi.nlm.nih.gov/nuccore/KF498959.1

gba <- GBAccession("KF498959.1")
gb <- readGenBank(gba, partial=TRUE)

# echo -e 'GFP\tunknown\texon\t1\t922\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf
# GFP     unknown exon    1       922     .       +       .       gene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";

genes <- genes(gb)

## GTF format
# https://useast.ensembl.org/info/website/upload/gff.html
# http://mblab.wustl.edu/GTF22.html

lapply(genes, 
       function(x, contig_name) {
          paste(contig_name, 
                "uknown", 
                "exon", 
                x@ranges@start, 
                x@ranges@start + x@ranges@width - 1, 
                sep = "\t") 
       }, 
       "RE")


