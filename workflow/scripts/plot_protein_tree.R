library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(readxl)
library(phangorn)
library(ggnewscale)
library(ggplot2)

with(snakemake@input, {
    newick_file <<- newick # analysis/rhodopsins/RAxML_bipartitions.txt
    metadata_file <<- metadata # metadata/rhodopsins.xlsx
})
# jplace_file <- "analysis/proteins_cat_phylogeny/pplacer.jplace"
# newick_file <- "analysis/proteins_cat_phylogeny/pplacer.newick"
# metadata_file <- "metadata/genomes.xlsx"
# part_ref_file <- "analysis/proteins_cat_phylogeny/treeshrink.part"
# part_query_file <- "analysis/proteins_cat_phylogeny/pplacer_input.part"
# aln_file <- "analysis/proteins_cat_phylogeny/pplacer_input.fasta"

add_mrca <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    mrca <- mutate(tree, my_column = !!as.name(colname)) %>%
        mutate(my_column = ifelse(my_column == "", NA, my_column)) %>%
        group_by(my_column) %>%
        mutate(is.tip = label %in% treedata@phylo$tip.label) %>%
        mutate(no_data = all(is.na(my_column))) %>%
        mutate(mrca = get_mrca(treedata@phylo, node[is.tip])) %>%
        mutate(mrca = ifelse(no_data | is.na(mrca), node, mrca)) %>%
        group_by(mrca) %>%
        mutate(enough_tips = sum(is.tip) > 0) %>%
        mutate(ifelse(node == mrca & enough_tips, first(na.omit(my_column)), NA)) %>%
        pull
    tree[[paste0(colname, "_mrca")]] <- mrca
    return(tree)
}
get_mrca <- function(phylo, tips) {
    getMRCA(phylo, tips) %>%
        replace(is.null(.), NA)
}
to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}

metadata <- read_excel(metadata_file)
tree <- read.tree(newick_file) %>%
    midpoint %>%
    as_tibble %>%
    mutate(bootstrap = ifelse(node %in% parent, label, NA)) %>%
    mutate(bootstrap = as.numeric(bootstrap)) %>%
    left_join(metadata, by = "label") %>%
    #mutate(name = ifelse(label == name, label, sprintf("%s (%s)", label, name))) %>%
    add_mrca(clade)
wrap_float <- function(x) {
    recode(format(round(x, 2), nsmall = 2), "  NA" = "")
}

shapes <- list(
    water = "circle",
    sediment = "square"
)
p <- ggtree(to_treedata(tree), layout = "rectangular") +
    scale_color_manual(values = c("black", "red")) +
    new_scale_color() +
    geom_text2(aes(subset = bootstrap > 1, x = branch, label = bootstrap), color = "black", size = 2, vjust = -0.5) +

    #geom_tippoint(aes(subset = !is.na(rhodopsin), color = rhodopsin, x = x + 0.025), shape = "diamond") +

    #new_scale_color() + new_scale("shape") +
    geom_tippoint(aes(subset = !is.na(activity), color = activity, x = x + 0.050), size = 1) +
    # scale_shape_manual(values = shapes) +

    #new_scale_fill() +
    #geom_fruit(data = part_ref, geom = geom_bar, mapping = aes(y = genome, x = num, fill = type), orientation = "y", stat = "identity") +
    #new_scale_fill() +
    #geom_fruit(data = part_query, geom = geom_tile, mapping = aes(y = genome, x = gene, fill = gene), axis.params = c(axis = "x", text = "gene", text.size = 2, text.angle = 90, hjust = 1)) +
    
    geom_tiplab(aes(label = alias), size = 2, offset = 0.075) +
    geom_cladelab(mapping = aes(subset = !is.na(clade_mrca), node = node, label = clade_mrca), align = T, offset = 0.6) +
    geom_treescale(width = 0.2) +
    theme(legend.position = "bottom")

ggsave(plot_file, p, width = 3, height = 5)
