library <- function(...) suppressPackageStartupMessages(base::library(...))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(ggforce)
library(tools)
library(grid)
library(ggstar)
library(ggnewscale)
library(sf)

with(snakemake@input, {
    genomes_file <<- genomes
    shape_dir <<- shape
    logan_file <<- logan
})
with(snakemake@output, {
    plot_file <<- plot
    data_file <<- data
})
with(snakemake@params, {
    min_coverage <<- min_coverage
    taxon <<- taxon
})

logan_data <- read.csv(logan_file, header = T) %>%
    mutate(data_type = ifelse(as.logical(is_transcriptome), "Metatranscriptome", "Metagenome"))

logan_data_present <- filter(logan_data, cov_fraction >= min_coverage)
logan_data_absent  <- filter(logan_data, cov_fraction <  min_coverage)

genomes <- read_excel(genomes_file) %>%
    filter(is.na(redundant), grepl(taxon, classification)) %>%
    mutate(rhodopsin = gsub("[*]", "", rhodopsin))

my_ggplot_world <- function(shape_dir) {
    world.land <- st_read(shape_dir)
    ggplot() +
        geom_sf(data = world.land, color = NA, fill = "lightgray") +
        coord_sf(expand = F) +
        scale_x_continuous(breaks = seq(-180, 180, 60)) +
        scale_y_continuous(breaks = seq(-90, 90, 30))
}

width <- 16
height <- 8

colors <- list(
    Metagenome = "#0c3d8aff",
    Metatranscriptome = "#00c0ffff"
)

write.csv(logan_data_present, data_file)
p <- my_ggplot_world(shape_dir) +
    # geom_point(aes(x = lon_parsed, y = lat_parsed, color = data_type), size = 0.01, shape = 20, data = logan_data_absent) +
    geom_point(aes(x = lon_parsed, y = lat_parsed, size = cov_fraction, color = data_type), data = logan_data_present) +
    # scale_size(nice.breaks = T) + # breaks = c(min_coverage, 0.01, 0.1, 1)) +
    scale_size(limits = c(min_coverage, 1), breaks = c(min_coverage, min_coverage * 4, min_coverage * 16, 0.25, 0.5, 1)) + 
    scale_color_manual(values = colors) +
    new_scale("color") +
    geom_text_repel(aes(x = lon, y = lat, label = name, color = rhodopsin), data = genomes, max.overlaps = Inf, nudge_y = 5, force = 10, min.segment.length = 0, point.size = 0, max.iter = Inf) +
    scale_colour_brewer(palette = "Set1", na.value = "black") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot_file, p, width = width, height = height)
