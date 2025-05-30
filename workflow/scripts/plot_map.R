library <- function(...) suppressPackageStartupMessages(base::library(...))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(ggforce)
library(tools)
library(sf)
library(grid)
library(ggstar)
library(ggnewscale)

with(snakemake@input, {
    shape_dir <<- shape # "analysis/maps/land/ne_110m_land.shp"
    jgi_img_profiles_files <<- jgi_img_profiles # "analysis/profiles/JGI_IMG_unrestricted_profiles.csv"
    jgi_img_matches_files  <<- jgi_img_matches # "analysis/profiles/JGI_IMG_unrestricted_matches.csv"
    jgi_img_samples_files  <<- jgi_img_samples # "analysis/profiles/JGI_IMG_unrestricted_samples.csv"
    om_rgc_profiles_files  <<- om_rgc_profiles # "analysis/profiles/OM-RGC_v2_orfs_profiles.csv"
    om_rgc_matches_files   <<- om_rgc_matches # "analysis/profiles/OM-RGC_v2_orfs_matches.csv"
    om_rgc_samples_files   <<- om_rgc_samples # "analysis/profiles/OM-RGC_v2_orfs_samples.csv"
})
with(snakemake@output, {
    plot_file <<- plot
    om_rgc_file <<- om_rgc
    jgi_img_file <<- jgi_img
    all_profiles_file <<- all_profiles
})
target_clade <- unlist(snakemake@params['clade'])
total_clade <- unlist(snakemake@params['total_clade']) # can be NULL in which case the sum of all clades is taken
max_depth <- unlist(snakemake@params['max_depth'])
size_breaks <- unname(unlist(snakemake@params['size_breaks']))

jgi_img_matches <- bind_rows(lapply(jgi_img_matches_files, read.csv)) %>%
    distinct(Scaffold.Name, .keep_all = T)
jgi_img_samples <- bind_rows(lapply(jgi_img_samples_files, read.csv, colClasses = c(Depth.In.Meters = "character"), na.strings = "_")) %>%
    distinct(Genome.ID, .keep_all = T)
jgi_img_profiles <- bind_rows(lapply(jgi_img_profiles_files, read.csv))

exclude <- "\\b(deep|subsurface|sediment|mesocosm|aphotic|oxygen minimum|O2min|bottom|250 m|enrichment|brine|lake|anoxic|anoxygenic|soil|brackish|crenothrix|cyanobacterial|enviromental|freshwater|glacier|spring|hypersaline|ice|lentic|wetlands|pond|alkaline|saline|watersheds)\\b"
rescue <- "\\b(Freshwater to marine saline gradient|Bioluminescent Bay)\\b"
jgi_img <- left_join(jgi_img_profiles, jgi_img_samples, by = "Genome.ID") %>%
    left_join(jgi_img_matches, by = "Scaffold.Name") %>%
    mutate(Description = paste(Isolation, Habitat, Ecosystem, Ecosystem.Category, Ecosystem.Subtype, Ecosystem.Type, Specific.Ecosystem)) %>%
    mutate(Depth.Num = as.numeric(Depth.In.Meters)) %>%
    filter(Ecosystem.Category == "Aquatic", Depth.Num <= max_depth | is.na(Depth.Num), grepl(rescue, Description, i = T) | !grepl(exclude, Description, i = T)) %>%
    mutate(Sample = as.character(Genome.ID))
write.csv(jgi_img, row.names = F, file = jgi_img_file)

layers <- c("SRF", "DCM")
om_rgc_profiles <- bind_rows(lapply(om_rgc_profiles_files, read.csv))
om_rgc_matches <- bind_rows(lapply(om_rgc_matches_files, read.csv)) %>%
    distinct(om_rgc, .keep_all = T)
om_rgc_samples <- bind_rows(lapply(om_rgc_samples_files, read.csv)) %>%
    distinct(sample, .keep_all = T)

om_rgc <- left_join(om_rgc_profiles, om_rgc_samples, by = "sample") %>%
    left_join(om_rgc_matches, by = c(OMRGC_ID = "om_rgc")) %>%
    filter(Layer %in% layers) %>%
    mutate(Sample = sample, Depth.In.Meters = Depth..nominal)
write.csv(om_rgc, row.names = F, file = om_rgc_file)

pct_value <- function(value, clade, target_clade, total_clade) {
    sum_value <- ifelse(is.null(total_clade), sum(value), sum(value[clade == total_clade]))
    ifelse(any(clade == target_clade), 100 * value[clade == target_clade] / sum(value), 0)
}

all_profiles <- bind_rows(list(om_rgc = om_rgc, jgi_img = jgi_img), .id = "Source") %>%
    group_by(Sample, Source, Sequencing.Strategy, Latitude, Longitude, Layer, Depth.In.Meters, clade) %>%
    summarize(value = sum(value), .groups = "drop_last") %>%
    summarize(target_pct = pct_value(value, clade, target_clade, total_clade), .groups = "drop") # ratio is the abundance of the target_clade relative to all the other clades
write.csv(all_profiles, row.names = F, file = all_profiles_file)

get_repel_coords <- function(.data, map_g, width, height, ...) {
    grid.newpage()
    Total <- 1
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_point(aes(x, y), data = .data) +
        geom_text_repel(aes(x, y), data = .data, max.overlaps = Inf, ...)
    panel_params <- ggplot_build(g)$layout$panel_params[[1]]
    xrg <- panel_params$x_range
    yrg <- panel_params$y_range

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
my_ggplot_world <- function(shape_dir) {
    world.land <- st_read(shape_dir)
    ggplot() +
        geom_sf(data = world.land, color = NA, fill = "lightgray") +
        coord_sf(expand = F) +
        scale_x_continuous(breaks = seq(-180, 180, 60)) +
        scale_y_continuous(breaks = seq(-90, 90, 30))
}

width <- 12
height <- 9

g <- my_ggplot_world(shape_dir)

profile_data <- rename(all_profiles, x = Longitude, y = Latitude) %>%
    filter(!is.na(x), !is.na(y)) # , x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])
profile_data_non_zero <- filter(profile_data, target_pct > 0) %>%
    get_repel_coords(g, width, height, label = "o", size = 0.5, force_pull = 20) %>%
    arrange(Sequencing.Strategy)
profile_data_zero <- group_by(profile_data, Source, Sequencing.Strategy, x, y) %>%
    filter(all(target_pct == 0)) %>%
    summarize(.groups = "drop") %>%
    get_repel_coords(g, width, height, label = "o", size = 0.5, force_pull = 20)

colors <- list(
    Metagenome = "#0c3d8aff",
    Metatranscriptome = "#00c0ffff"
)
colors_fill <- list(
    Metagenome = "#6d88c4ff",
    Metatranscriptome = NA
)
shapes_non_zero <- c(
    om_rgc = 22,
    jgi_img = 23
)
shapes_zero <- c(
    om_rgc = 3,
    jgi_img = 4
)

size_lim_1 <- min(profile_data_non_zero$target_pct)
size_lim_2 <- max(profile_data_non_zero$target_pct, size_breaks)

p <- g +

    geom_point(aes(x = x, y = y, color = Sequencing.Strategy, shape = Source), size = 1.5, stroke = 0.3, data = profile_data_zero) +
    scale_shape_manual(values = shapes_zero) +
    scale_color_manual(values = colors) +
    new_scale("shape") +

    geom_point(aes(x = x, y = y, color = Sequencing.Strategy, fill = Sequencing.Strategy, size = target_pct, shape = Source), stroke = 0.7, data = profile_data_non_zero) +
    scale_shape_manual(values = shapes_non_zero) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors_fill, na.value = NA) +
    scale_size(range = c(2.5, 14), breaks = rep(size_breaks, each = 2), limits = c(size_lim_1, size_lim_2), guide = guide_legend(ncol = 2, override.aes = list(shape = shapes_non_zero))) +

    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")
ggsave(plot_file, p, width = width, height = height)
