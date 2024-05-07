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
    # tree_file <<- newick # analysis/rhodopsins/RAxML_bipartitions.txt
    network_file <<- network # analysis/rhodopsins/splitstree.nex
    metadata_file <<- metadata # metadata/rhodopsins.xlsx
    residues_file <<- residues # analysis/rhodopsins/residues.txt
})
plot_file <- unlist(snakemake@output)

min_support <- 50

network <- read.nexus.networx(network_file)

synonyms <- list(
    `chloride pump` = "anion pump",
    `[proton] pump` = "proton pump"
)
residues <- read.table(residues_file, header = T)
node_data <- read_excel(metadata_file) %>%
    left_join(residues, by = "label") %>%
    mutate(activity = recode(activity, !!!synonyms)) %>%
    mutate(node = 1 + match(label, network$tip.label)) %>% # I do not quite understand why we need to add 1 here
    mutate(fen1 = str_sub(fenestration, 1, 1), fen2 = str_sub(fenestration, 2, 2)) %>%
    mutate(alias = ifelse(is.na(alias), label, alias)) %>%
    rename(lab = label)

aliases <- with(node_data, setNames(alias, lab))
network$translate$label <- aliases[network$tip.label]  

# Support values from the tree are not yet used
#tree <- read.tree(tree_file)
#supports <- with(tree, c(rep(NA, length(tip.label)), node.label))
#edge_data <- data.frame(parent = network$edge[,1], node = network$edge[,2], tree.node = createLabel(network, tree, tree$edge[,2], "edge")) %>%
#    mutate(support = as.numeric(supports[tree.node])) %>%
#    mutate(support = ifelse(support < min_support, NA, support))

x_ext <- Vectorize(function(x0, y0, x1, y1, d) {
    dx <- x1 - x0
    dy <- y1 - y0
    return(x1 + dx * d / sqrt(dx^2 + dy^2))
}, vectorize.args = c("x0", "y0", "x1", "y1"))
y_ext <- Vectorize(function(x0, y0, x1, y1, d) {
    dx <- x1 - x0
    dy <- y1 - y0
    return(y1 + dy * d / sqrt(dx^2 + dy^2))
}, vectorize.args = c("x0", "y0", "x1", "y1"))

add_edge_data <- function(p, more_data) {
    p$data <- left_join(p$data, more_data, by = c("node", "parent"))
    return(p)
}
add_node_data <- function(p, more_data) {
    p$data <- left_join(p$data, more_data, by = "node")
    return(p)
}

activity_colors <- list(
    `proton pump` = "#00ba38ff",
    `sodium pump` = "#619cffff",
    `anion pump` = "#f8766dff",
    `sensory` = "#ff00ffff"
)
clade_colors <- list(
    NQ = "#00c266ff",
    P1 = "#00c2c2ff",
    XR = "#ff5cccff",
    ACB = "#f57066ff",
    PR = "#00a3ffff",
    ESR = "#d4aa00ff",
    P4 = "#66ff00ff",
    MACR = "#a7ac93ff",
    HeimdallR = "#cc00ffff",
    `Proteo-SR` = "#ff2ad4ff",
    TwR = "#c83737ff"
)
antenna_colors <- list(
    `=O` = "blue",
    `-OH` = "green"
)
fen_colors <- list(
    G = "#ffa500ff",
    W = "#0000ffff",
    F = "#0000ffff",
    Y = "#0000ffff",
    L = "#00ffffff"
)
motif_colors <- list(
    DTE = "#00ba38ff",
    DTK = "#00ba38ff",
    DTT = "#00ba38ff",
    NDQ = "#619cffff",
    NTQ = "#f8766dff"
)

p <- ggsplitnet(network, size = 0.1) %>%
    # add_edge_data(edge_data) %>%
    add_node_data(node_data) +
    geom_point2(aes(subset = !is.na(activity), x = x, y = y, color = activity)) +
    scale_color_manual(values = activity_colors) +
    new_scale_color() +
    geom_tiplab2(aes(label = motif, color = motif, x = x_ext(xend, yend, x, y, 0.01), y = y_ext(xend, yend, x, y, 0.01)), size = 2) +
    # scale_color_manual(values = "red", labels = "DTE") +
    new_scale_color() +
    geom_tiplab2(aes(label = fen1,  color = fen1,  x = x_ext(xend, yend, x, y, 0.055), y = y_ext(xend, yend, x, y, 0.055)), size = 2.5) +
    scale_color_manual(values = fen_colors) +
    new_scale_color() +
    geom_tiplab2(aes(label = alias, x = x_ext(xend, yend, x, y, 0.08), y = y_ext(xend, yend, x, y, 0.08)), color = "gray30") +
    new_scale_color() +
    geom_tiplab2(aes(subset = !is.na(antenna_type), x = x_ext(xend, yend, x, y, 0.2), y = y_ext(xend, yend, x, y, 0.2), color = antenna_type, label = antenna_type), size = 3) +
    scale_color_manual(values = antenna_colors) +
    new_scale_color() +
    geom_treescale(width = 0.1) +
    coord_fixed()
ggsave(plot_file, p, width = 10, height = 8)
