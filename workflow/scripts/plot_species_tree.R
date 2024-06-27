library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(readxl)
library(phangorn)
library(ggnewscale)
library(ggplot2)
library(ggtreeExtra)

with(snakemake@input, {
    jplace_file <<- jplace
    newick_file <<- newick
    metadata_file <<- metadata
    genes_file <<- genes
})
plot_file <- unlist(snakemake@output)
# jplace_file <- "analysis/pplacer/pplacer.jplace"
# newick_file <- "analysis/pplacer/pplacer.newick"
# metadata_file <- "metadata/genomes.xlsx"
# genes_file <- "analysis/pplacer/pplacer_input.txt"

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
genes <- read.table(genes_file, col.names = c("genome", "gene", "set"))
ref_genes <- filter(genes, set == "reference") %>%
    mutate(total = n_distinct(gene)) %>%
    group_by(genome) %>%
    summarize(present = n(), missing = first(total) - present, .groups = "drop") %>%
    gather(type, num, -genome)
query_genes <- filter(genes, set == "query")

metadata <- read_excel(metadata_file) %>%
    separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")
jplace <- read.jplace(jplace_file)

lwr <- arrange(jplace@placements, -likelihood) %>%
    select(label = name, like_weight_ratio) %>%
    distinct(label, .keep_all = T)

tree <- read.tree(newick_file) %>%
    midpoint %>%
    as_tibble %>%
    mutate(bootstrap = ifelse(node %in% parent, label, NA)) %>%
    mutate(bootstrap = as.numeric(bootstrap)) %>%
    left_join(lwr, by = "label") %>%
    left_join(metadata, by = c(label = "genome")) %>%
    mutate(name = ifelse(label == name, label, sprintf("%s (%s)", label, name))) %>%
    add_mrca(genus) %>%
    add_mrca(species)
wrap_float <- function(x) {
    recode(format(round(x, 2), nsmall = 2), "  NA" = "")
}

shapes <- list(
    water = "circle",
    sediment = "square"
)
p <- ggtree(to_treedata(tree), aes(color = !is.na(like_weight_ratio)), layout = "rectangular") +
    scale_color_manual(values = c("black", "red")) +
    new_scale_color() +
    geom_text2(aes(subset = bootstrap > 90, x = branch, label = bootstrap), color = "black", size = 2, vjust = -0.5) +
    geom_text2(aes(x = branch, label = wrap_float(like_weight_ratio)), color = "red", size = 2, vjust = 1.5) +
    geom_tippoint(aes(subset = !is.na(rhodopsin), color = rhodopsin, x = x + 0.025), shape = "diamond") +

    new_scale_color() + new_scale("shape") +
    geom_tippoint(aes(shape = habitat1, color = habitat2, x = x + 0.050)) +
    # scale_shape_manual(values = shapes) +

    new_scale_fill() +
    geom_fruit(data = ref_genes, geom = geom_bar, mapping = aes(y = genome, x = num, fill = type), orientation = "y", stat = "identity") +
    new_scale_fill() +
    geom_fruit(data = query_genes, geom = geom_tile, mapping = aes(y = genome, x = gene, fill = gene), axis.params = c(axis = "x", text = "gene", text.size = 2, text.angle = 90, hjust = 1)) +
    
    geom_tiplab(aes(label = name), size = 2, offset = 0.075) +
    geom_cladelab(mapping = aes(subset = !is.na(genus_mrca), node = node, label = genus_mrca), align = T, offset = 0.6) +
    geom_treescale(width = 0.2) +
    theme(legend.position = "bottom") +
    xlim(0, 2)

ggsave(plot_file, p, width = 5, height = 5)
