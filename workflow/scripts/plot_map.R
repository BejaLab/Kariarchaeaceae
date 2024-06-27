library <- function(...) suppressPackageStartupMessages(base::library(...))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(ggforce)
library(tools)
library(rgdal)
library(grid)
library(ggstar)
library(ggnewscale)

with(snakemake@input, {
    genomes_file <<- genomes # "metadata/genomes.xlsx"
    shape_file <<- shape # "analysis/maps/land/ne_110m_land.shp"
    jgi_img_profiles_file <<- jgi_img_profiles # "analysis/profiles/JGI_IMG_unrestricted_profiles.csv"
    jgi_img_matches_file  <<- jgi_img_matches # "analysis/profiles/JGI_IMG_unrestricted_matches.csv"
    jgi_img_samples_file  <<- jgi_img_samples # "analysis/profiles/JGI_IMG_unrestricted_samples.csv"
    om_rgc_profiles_file  <<- om_rgc_profiles # "analysis/profiles/OM-RGC_v2_orfs_profiles.csv"
    om_rgc_matches_file   <<- om_rgc_matches # "analysis/profiles/OM-RGC_v2_orfs_matches.csv"
    om_rgc_samples_file   <<- om_rgc_samples # "analysis/profiles/OM-RGC_v2_orfs_samples.csv"
})
output_file <- unlist(snakemake@output)
target_clade <- unlist(snakemake@params['clade'])

jgi_img_matches <- read.csv(jgi_img_matches_file)
jgi_img_samples <- read.csv(jgi_img_samples_file, na.strings = "_")
jgi_img_profiles <- read.csv(jgi_img_profiles_file)

exclude <- "\\b(deep|subsurface|sediment|mesocosm|aphotic|oxygen minimum|O2min|bottom|250 m)\\b"
max_depth <- 200
jgi_img <- left_join(jgi_img_profiles, jgi_img_samples, by = "Genome.ID") %>%
    left_join(jgi_img_matches, by = "Scaffold.Name") %>%
    mutate(Description = paste(Isolation, Habitat, Ecosystem, Ecosystem.Category, Ecosystem.Subtype, Ecosystem.Type, Specific.Ecosystem)) %>%
    filter(Ecosystem.Category == "Aquatic", Depth.In.Meters <= max_depth | is.na(Depth.In.Meters), !grepl(exclude, Description, i = T)) %>%
    mutate(Sample = as.character(Genome.ID))

layers <- c("SRF", "DCM")
om_rgc_profiles <- read.csv(om_rgc_profiles_file)
om_rgc_matches <- read.csv(om_rgc_matches_file)
om_rgc_samples <- read.csv(om_rgc_samples_file)

om_rgc <- left_join(om_rgc_profiles, om_rgc_samples, by = "sample") %>%
    left_join(om_rgc_matches, by = c(OMRGC_ID = "om_rgc")) %>%
    filter(Layer %in% layers) %>%
    mutate(Sample = sample)

all_profiles <- bind_rows(list(om_rgc = om_rgc, jgi_img = jgi_img), .id = "Source") %>%
    group_by(Source, Sequencing.Strategy, Latitude, Longitude, Layer, clade) %>%
    summarize(value = sum(value), .groups = "drop_last") %>%
    summarize(target_pct = ifelse(any(clade == target_clade), 100 * value[clade == target_clade] / sum(value), 0), .groups = "drop") # ratio is the abundance of the target_clade relative to all the other clades

genomes <- read_excel(genomes_file) %>%
    filter(!is.na(rhodopsin), is.na(redundant)) %>%
    mutate(rhodopsin = gsub("[*]", "", rhodopsin))

get_repel_coords <- function(.data, map_g, width, height, ...) {
    grid.newpage()
    Total <- 1
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_point(aes(x, y), data = .data) +
        geom_text_repel(aes(x, y), data = .data, max.overlaps = Inf, ...)
    panel_params <- ggplot_build(g)$layout$panel_params[[1]]
    xrg <- panel_params$x.range
    yrg <- panel_params$y.range

    textrepeltree <- ggplotGrob(g) %>%
        grid.force(draw = F) %>%
        getGrob("textrepeltree", grep = T)
    children <- childNames(textrepeltree) %>%
        grep("textrepelgrob", ., value = T)

    get_xy <- function(n) {
        grob <- getGrob(textrepeltree, n)
        data.frame(
            x.repel = xrg[1] + diff(xrg) * convertX(grob$x, "native", valueOnly = T),
            y.repel = yrg[1] + diff(yrg) * convertY(grob$y, "native", valueOnly = T)
        )
    }
    lapply(children, get_xy) %>%
        bind_rows %>%
        cbind(.data) %>%
        mutate(theta = atan2(y - y.repel, x - x.repel), x.segm = x.repel + Total * cos(theta), y.segm = y.repel + Total * sin(theta))
}
my_ggplot_world <- function(shape_file, xlim = c(-185, 185), ylim = c(-85, 90)) {
    map_base <- file_path_sans_ext(shape_file)
    world.land <- readOGR(dirname(shape_file), basename(file_path_sans_ext(shape_file)))
    continents.simple <- fortify(world.land)
    all.continents <- data_frame(id = unique(as.character(continents.simple$id)))
    ggplot(all.continents) +
        geom_map(data = continents.simple, map = continents.simple, aes(map_id = id), color = "lightgray", fill = "lightgray") +
        coord_quickmap(xlim = xlim, ylim = ylim, expand = F)
}

xlim <- c(-185, 185)
ylim <- c(-85, 90)

g <- my_ggplot_world(shape_file, xlim, ylim)

width <- 16
height <- 8

profile_data <- rename(all_profiles, x = Longitude, y = Latitude) %>%
    filter(!is.na(x), !is.na(y), x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])
profile_data_non_zero <- filter(profile_data, target_pct > 0) %>%
    get_repel_coords(g, width, height, label = "o", size = 0.5, force_pull = 20)
profile_data_zero <- group_by(profile_data, Source, Sequencing.Strategy, x, y) %>%
    filter(all(target_pct == 0)) %>%
    summarize(.groups = "drop") %>%
    get_repel_coords(g, width, height, label = "o", size = 0.5, force_pull = 20)

genome_data <- filter(genomes, lon >= xlim[1], lon <= xlim[2], lat >= ylim[1], lat <= ylim[2])

colors <- list(
    Metagenome = "#0c3d8aff",
    Metatranscriptome = "#d4aa00ff"
)
shapes_non_zero <- c(
    om_rgc = 0,
    jgi_img = 5
)
shapes_zero <- c(
    om_rgc = 3,
    jgi_img = 4
)
p <- g +
    geom_star(aes(x = lon, y = lat, color = rhodopsin, size = 2 - is.na(rhodopsin)), starshape = 1, data = genome_data) +
    scale_size(range = c(3, 5)) +
    new_scale("shape") + new_scale("color") + new_scale("fill") + new_scale("size") +
    geom_text_repel(aes(x = lon, y = lat, label = name, color = is.na(rhodopsin)), data = genome_data, max.overlaps = Inf, min.segment.length = 0, point.size = 5, max.iter = Inf, force = 20, force_pull = 50) +
    scale_color_manual(values = list("TRUE" = "darkgray", "FALSE" = colors$Metagenome)) +
    new_scale("color") + new_scale("size") +

    geom_point(aes(x = x.repel, y = y.repel, color = Sequencing.Strategy, shape = Source), size = 1, stroke = 0.1, data = profile_data_zero) +
    scale_shape_manual(values = shapes_zero) +
    scale_color_manual(values = colors) +
    new_scale("shape") +

    geom_point(aes(x = x.repel, y = y.repel, color = Sequencing.Strategy, size = target_pct, shape = Source), data = profile_data_non_zero) +
    scale_shape_manual(values = shapes_non_zero) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(2, 10), guide = guide_legend(override.aes = list(shape = 0))) +

    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(output_file, p, width = width, height = height)
