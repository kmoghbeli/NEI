## Utility file to make GTF file from a Genbank record
library(tidyverse)

data_dir <- "../data_objects/viral_genomes/"

virus <- "re"

gtf <- readr::read_delim(paste0(data_dir, virus, ".fixed2.gtf"), 
                         col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"), 
                         comment = "#")

# Now add the "gene_id", "transcript_id", and "gene_name" attributes for all "exon" features (per the note here:) 
# https://kb.10xgenomics.com/hc/en-us/articles/115003327112-How-can-we-add-genes-to-a-reference-package-for-Cell-Ranger


#gsub(".*gene=([^;]+);.*", '\\1', "afsadsafafafdsa;gene=XYZ.43943XXX.d;sfajskffasklfajkfj=")

gtf.fixed <- gtf %>% 
  filter("exon" == feature) %>% 
  mutate(gene = gsub(".*gene=([^;]+);.*", '\\1', attributes), 
         gene = paste0("HSV-", gene), 
         attributes = paste0("gene_id \"", gene, "\"; transcript_id \"", gene, "\"; gene_name \"", gene, "\"")) %>% 
  select(!gene)
                                     
gtf.fixed %>% readr::write_delim(paste0(data_dir, virus, ".fixed.gtf"), delim = "\t", col_names = FALSE, quote = "none", escape = "none")
  

#------------------------------------------------------------------------------------------------------------------------------------------------
## Create GTF file to use with this tutorial: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
# KOS Genome: https://www.ncbi.nlm.nih.gov/nuccore/JQ673480
# RE Genome: https://www.ncbi.nlm.nih.gov/nuccore/ON960060.1
# (OLD INCORRECT) RE Genome: https://www.ncbi.nlm.nih.gov/nuccore/KF498959.1

library(genbankr)
library(rentrez)

gba <- GBAccession("ON960060.1")
gb <- readGenBank(gba, partial=TRUE)

# echo -e 'GFP\tunknown\texon\t1\t922\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf
# GFP     unknown exon    1       922     .       +       .       gene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";

genes <- genes(gb)

## GTF format
# https://useast.ensembl.org/info/website/upload/gff.html
# http://mblab.wustl.edu/GTF22.html

gtf_addition <-
  paste(
    lapply(seq_along(genes),
           function(i, contig_name) {
             paste(contig_name,
                   "KMgtfmake",
                   "exon",
                   genes[i]@ranges@start,
                   genes[i]@ranges@start + genes[i]@ranges@width - 1,
                   ".", 
                   as.character(genes[i]@strand@values), 
                   ".", 
                   paste0("gene_id \"", genes[i]$gene, "; transcript_id \"", genes[i]$gene, "; gene_name \"", genes[i]$gene, "; gene_biotype \"protein_coding\""),
                   sep = "\t")
           },
           "RE"),
    sep = "",
    collapse = "\n")

write(gtf_addition, file = paste0(data_dir, "re.fixed2.gtf"))


