library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(readxl)
library(phangorn)
library(ggnewscale)
library(ggplot2)
library(tanggle)
library(stringr)

with(snakemake@input, {
    tree_file <<- newick # analysis/rhodopsins/RAxML_bipartitions.txt
    network_file <<- nexus # analysis/rhodopsins/splitstree.nex
    metadata_file <<- metadata # metadata/rhodopsins.xlsx
    residues_file <<- residues # analysis/rhodopsins/residues.txt
})
with(snakemake@output, {
    plot_file <<- plot
    jtree_file <<- jtree
})

min_support <- 50
color_range <- c("blue", "red")
color_na <- "gray"

network <- read.nexus.networx(network_file)
tree <- read.tree(tree_file)

residues <- read.table(residues_file, header = T)
node_data <- read_excel(metadata_file) %>%
    left_join(residues, by = "label") %>%
    mutate(node = match(label, network$tip.label)) %>%
    mutate(fenestration_residue = str_sub(fenestration, 1, 1))

network$translate$label <- filter(node_data, !is.na(node)) %>%
    arrange(node) %>%
    pull(alias)

supports <- with(tree, c(rep(NA, length(tip.label)), node.label))
edge_data <- data.frame(parent = network$edge[,1], node = network$edge[,2], tree.node = createLabel(network, tree, tree$edge[,2], "edge")) %>%
    mutate(support = as.numeric(supports[tree.node])) %>%
    mutate(support = ifelse(support < min_support, NA, support))

add_edge_data <- function(p, more_data) {
    p$data <- left_join(p$data, more_data, by = c("node", "parent"))
    return(p)
}

p <- ggsplitnet(network) %>% add_edge_data(edge.labels) +
    aes(color = support) +
    scale_colour_gradient(low = "blue", high = "red") +
    geom_tiplab2()
ggsave(plot_file, p)


treedata <- as_tibble(tree) %>%
    mutate(support = ifelse(node %in% parent, label, NA)) %>%
    mutate(support = as.numeric(support))
p <- ggtree(as.treedata(tree), layout = "ape") +
    geom_text2(aes(subset = support > 50, x = branch, label = bootstrap), color = "black", size = 2, vjust = -0.5) +

    #geom_tippoint(aes(subset = !is.na(rhodopsin), color = rhodopsin, x = x + 0.025), shape = "diamond") +

    #new_scale_color() + new_scale("shape") +
    geom_tippoint(aes(subset = !is.na(activity), color = activity, x = x + 0.050), size = 1) +
    # scale_shape_manual(values = shapes) +

    #new_scale_fill() +
    #geom_fruit(data = part_ref, geom = geom_bar, mapping = aes(y = genome, x = num, fill = type), orientation = "y", stat = "identity") +
    #new_scale_fill() +
    #geom_fruit(data = part_query, geom = geom_tile, mapping = aes(y = genome, x = gene, fill = gene), axis.params = c(axis = "x", text = "gene", text.size = 2, text.angle = 90, hjust = 1)) +
    
    geom_tiplab(aes(label = alias), size = 2, offset = 0.075) +
    geom_treescale(width = 0.2) +
    theme(legend.position = "bottom")

ggsave(plot_file, p, width = 3, height = 5)
