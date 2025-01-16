library <- function(...) suppressPackageStartupMessages(base::library(...))
library(PopGenome)
VCF_split_into_scaffolds(unlist(snakemake@input), unlist(snakemake@output))
