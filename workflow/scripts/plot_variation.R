library <- function(...) suppressPackageStartupMessages(base::library(...))
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(gggenes)
library(scales)
library(ggrepel)

with(snakemake@input, {
    bed_file <<- bed
    gff_file <<- gff
    cov_file <<- cov
})
with(snakemake@params, {
    min_scaf_len <<- min_scaf_len
    min_cds_len <<- min_cds_len
    target <<- target
    win_size <<- win_size
    discard_dep_percentile <<- discard_dep_percentile
})
with(snakemake@output, {
    genome_file <<- genome
    scaffold_file <<- scaffold
    genome_csv_file <<- genome_csv
    scaffold_csv_file <<- scaffold_csv
})

get_attr <- function(.data, into) {
    extract(.data, attribute, into = into, regex = paste0("\\b", into, "=([^;]+)"), convert = T, remove = F)
}
products <- list(
    bacteriorhodopsin = "rhodopsin",
    `hypothetical protein` = ""
)
interval_join <- function(left, right) {
    mutate(data.table(left), start_pos = start, end_pos = end)[right, on = .(seqname, start_pos <= pos, end_pos >= pos)]
}
generate_sliding_windows <- function(interval_size, window_size) {
    lapply(1:(interval_size - window_size + 1), function(i) c(start = i, end = i + window_size - 1))
}
quantile_lines <- c(
    `5%` = "dotted",
    `25%` = "dashed",
    `50%` = "solid",
    `75%` = "dashed",
    `95%` = "dotted"
)
strand_colors <- c(
    `+` = "#0096ffff",
    `-` = "#808000ff"
)

col_names <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gff <- read.table(gff_file, sep = "\t", col.names = col_names, quote = "") %>%
    mutate(seqname = URLdecode(seqname)) %>%
    mutate(seqname = sub("_length=\\d+", "", seqname))
cds <- filter(gff, feature == "CDS") %>%
    get_attr("product") %>% get_attr("locus_tag") %>%
    extract("locus_tag", "gene_num", regex = "_(\\d+)$") %>%
    select(seqname, gene_num, strand, start, end, frame, product) %>%
    mutate(label = str_wrap(recode(product, !!!products), 20)) %>%
    mutate(cds_len = end - start + 1)
scaffolds <- filter(gff, feature == "region") %>%
    select(seqname, scaf_len = end)
coverage <- read.csv(cov_file) %>%
    separate(Locus, into = c("seqname", "pos"), sep = ":", convert = T)

col_names <- c("seqname", "start0", "end")
bed <- read.table(bed_file, sep = "\t", col.names = col_names) %>%
    mutate(seqname = sub("_length=\\d+", "", seqname)) %>%
    uncount(end - start0, .id = "offset") %>%
    mutate(pos = start0 + offset)

cds_data_bed <- interval_join(cds, bed) %>%
    filter(cds_len >= min_cds_len) %>%
    group_by(seqname, start, end, strand, gene_num, cds_len, product) %>%
    summarize(total = n(), .groups = "drop") %>%
    mutate(measure = "Number of variants per bp", value = total / cds_len)
cds_data_dep <- interval_join(cds, coverage) %>%
    filter(cds_len >= min_cds_len) %>%
    group_by(seqname, start, end, strand, gene_num, cds_len, product) %>%
    summarize(total = sum(Total_Depth), .groups = "drop") %>%
    mutate(measure = "Depth per bp", value = total / cds_len) %>%
    mutate(depth_too_high = value > quantile(value, discard_dep_percentile))
cds_data <- bind_rows(cds_data_bed, cds_data_dep) %>%
    left_join(scaffolds, by = "seqname") %>%
    mutate(is_target = gene_num %in% target) %>%
    mutate(pos = (start + end) / 2) %>%
    group_by(measure) %>%
    mutate(is_outlier = value > quantile(value, 0.75) + 1.5 * IQR(value))
cds_quants <- summarize(cds_data, value = quantile(value, c(0.25, 0.5, 0.75)), quantile = names(value), .groups = "drop")

genome_int <- 100000
scale_fun <- function(x) sprintf("%d", x)

cds_data_filtered <- filter(cds_data, scaf_len >= min_scaf_len)
p <- ggplot(cds_data_filtered, aes(x = pos, y = value)) +
    geom_line(linewidth = 0.1) +
    geom_point(aes(color = strand, size = is_target, shape = is_target)) +
    geom_point(aes(x = pos, y = value + 0.01), shape = 8, size = 2, data = filter(cds_data_filtered, is_outlier)) +
    geom_hline(aes(yintercept = value, linetype = quantile), color = "#008000ff", data = cds_quants) +
    scale_linetype_manual(values = quantile_lines) +
    scale_color_manual(values = strand_colors) +
    facet_grid(measure ~ seqname, scale = "free", space = "free_x", switch = "y") +
    scale_x_continuous(labels = scale_fun, breaks = function(z) if (z[2] > genome_int) c(1, seq(genome_int, range(z)[2], by = genome_int)) else c(1,genome_int)) +
    xlab("Position, bp") +
    theme_bw() +
    theme(axis.title.y = element_blank(), strip.placement = "outside", strip.background = element_blank())
ggsave(genome_file, height = 7, width = 15)
pivot_wider(cds_data,
    id_cols = c("seqname", "start", "end", "strand", "gene_num", "product"), names_from = "measure", values_from = c("total", "value", "is_outlier")
) %>% write.csv(file = genome_csv_file)

target_cds <- filter(cds, gene_num == target)

windows <- rowwise(scaffolds) %>%
    mutate(sliding_windows = list(generate_sliding_windows(scaf_len, win_size))) %>%
    unnest(cols = sliding_windows) %>%
    unnest_wider(sliding_windows)

win_data_bed <- interval_join(windows, bed) %>%
    group_by(seqname, start, end) %>%
    summarize(num_vars = n(), .groups = "drop") %>%
    mutate(measure = "Number of variants per bp", value = num_vars / (end - start + 1))
#win_data_dep <- interval_join(windows, coverage) %>%
#    group_by(seqname, start, end) %>%
#    summarize(total_depth = sum(Total_Depth), .groups = "drop") %>%
#    mutate(measure = "Depth per bp", value = total_depth / (end - start + 1))
win_data <- bind_rows(win_data_bed) %>%
    left_join(scaffolds, by = "seqname") %>%
    mutate(pos = (start + end) / 2)
win_quants <- group_by(win_data, measure) %>%
    summarize(value = quantile(value, c(0.05, 0.25, 0.5, 0.75, 0.95)), quantile = names(value), .groups = "drop")
scaf_cds <- filter(cds, seqname == target_cds$seqname) %>%
    mutate(pos = (start + end)/2)
scaffold_win_data <- filter(win_data, seqname == target_cds$seqname)

scaf_int <- 1000
gene_track_y <- max(win_data$value)
p <- ggplot(scaffold_win_data, aes(x = pos, y = value)) +
    geom_gene_arrow(aes(xmin = start, xmax = end, forward = strand == "+"), y = gene_track_y, arrowhead_height = unit(25, "mm"), arrow_body_height = unit(20, "mm"), data = scaf_cds) +
    geom_line(linewidth = 0.1) +
    geom_text(aes(x = pos, label = label), y = gene_track_y, data = scaf_cds, angle = 90, lineheight = .6) +
    geom_hline(aes(yintercept = value, linetype = quantile), color = "#008000ff", data = win_quants) +
    scale_linetype_manual(values = quantile_lines) +
    xlab("Position, bp") +
    ylab("Number of variants per bp") +
    ylim(c(0, gene_track_y)) +
    scale_x_continuous(labels = scale_fun, breaks = function(z) c(1, seq(scaf_int, range(z)[2], by = scaf_int))) +
    theme_bw()
ggsave(scaffold_file, height = 4, width = 15)

pivot_wider(scaffold_win_data, id_cols = c("seqname", "pos", "start", "end"), names_from = "measure", values_from = "value") %>%
    write.csv(file = scaffold_csv_file)
