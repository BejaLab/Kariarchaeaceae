library <- function(...) suppressPackageStartupMessages(base::library(...))

library(dplyr)
library(tidyr)
library(tools)

proteinortho_file <- unlist(snakemake@input)
with(snakemake@params, {
    genomes <<- genomes
    fraction <<- fraction
})
output_file <- unlist(snakemake@output)

counts <- read.table(proteinortho_file, header = T, comment.char = "", sep = "\t", check.names = F, na.strings = "*") %>%
    rename_with(~file_path_sans_ext(.)) %>%
    select(any_of(genomes)) %>%
    mutate(num_genomes = rowSums(!is.na(.))) %>%
    filter(num_genomes >= length(genomes) * fraction) %>%
    {colSums(!is.na(.))}
pcts <- counts %>%
    {100 * .[names(.) != "num_genomes"] / .["num_genomes"]} %>%
    {data.frame(genome = names(.), pct = unname(.))}
write.table(pcts, file = output_file, sep = "\t", col.names = F, row.names = F, quote = F)
