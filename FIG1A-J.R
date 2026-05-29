#Figure 1
#Figure 1A-J

pkgs <- c(
  "dplyr", "ggplot2", "patchwork", "sf", "rnaturalearth",
  "rnaturalearthdata", "colorspace", "scales", "vegan",
  "ggnewscale", "viridis", "ggridges"
)
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE))

get_script_dir <- function() {
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(cmd_file) > 0 && file.exists(cmd_file[1])) {
    return(dirname(normalizePath(cmd_file[1], winslash = "/", mustWork = TRUE)))
  }
  frame_files <- vapply(
    sys.frames(),
    function(x) if (!is.null(x$ofile)) x$ofile else NA_character_,
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0 && file.exists(tail(frame_files, 1))) {
    return(dirname(normalizePath(tail(frame_files, 1), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
output_dir <- script_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

find_file <- function(filename) {
  candidates <- c(
    file.path(script_dir, filename),
    file.path(dirname(script_dir), filename),
    filename,
    file.path("D:/R代码", filename)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("File not found: ", filename)
  }
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

save_panel <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
}

iso_list <- c(
  "JO","SD","AZ","EE","MZ","GT","MX","CM","TZ","ML","RO","NG","KZ","EC","TH","RS","VE",
  "RU","CF","MG","MW","PE","SG","BR","PL","TR","LT","HR","PT","UG","AE","IN","IT","ZW",
  "AU","IE","DK","ES","BD","CO","KE","HU","NL","MY","NO","GR","AM","CH","CA","SE","AT",
  "IS","FR","GB","NZ","FI","DE","IL","BE","CN","KR","US","GH","JP"
)

base_colors <- c(
  "#cc79a7", "#009e74", "#56b3e9", "#0071b2",
  "#000000", "#e69d00", "#f0e442", "#d55e00"
)
expanded_colors <- qualitative_hcl(
  n = 64,
  palette = "Dark 3",
  h = c(0, 360),
  c = 60,
  l = 65
)
names(expanded_colors) <- iso_list

stage_colors <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

meta_file <- find_file(if (file.exists(file.path(script_dir, "sample_metadata.tsv"))) "sample_metadata.tsv" else "Supplementary_Table_44_sample_metadata.tsv")
meta <- read.delim(meta_file, sep = "\t", check.names = FALSE) %>%
  dplyr::mutate(sample = paste(project, srr, sep = "_"))

country_stage <- read.csv(find_file("country_stage.csv"), stringsAsFactors = FALSE) %>%
  dplyr::transmute(iso = location_name, stage = as.character(stage))

regions_raw <- read.csv(find_file("regions.csv"), header = FALSE, stringsAsFactors = FALSE)
regions <- regions_raw[, seq_len(min(3, ncol(regions_raw)))]
colnames(regions) <- c("iso_a2", "country", "region")
regions <- regions %>%
  dplyr::filter(iso_a2 != "iso") %>%
  dplyr::mutate(iso_a2 = as.character(iso_a2))

filtered <- readRDS(find_file("filtered.rds"))
nmds <- readRDS(find_file("nmds.rds"))
points <- readRDS(find_file("pcoa_points.rds"))

points.annotated <- points %>%
  dplyr::left_join(meta, by = "sample") %>%
  dplyr::left_join(country_stage, by = c("iso" = "iso")) %>%
  dplyr::mutate(stage = factor(stage, levels = c("1", "2", "3")))

eigen_total <- sum(nmds$eigen) / 100
xlims <- c(-10.5, 12.5)
ylims <- c(-8, 11)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world$region_wb[world$admin %in% c("Australia", "New Zealand")] <- "Australia and New Zealand"

incidence_map_df <- data.frame(
  location_name = c(
    "Australia and New Zealand",
    "Europe & Central Asia",
    "North America",
    "South Asia",
    "Middle East & North Africa",
    "Latin America & Caribbean",
    "East Asia & Pacific",
    "Sub-Saharan Africa",
    "Oceania"
  ),
  incidence = c(
    23.8575598, 16.3063269, 16.3063269, 5.2794844,
    3.3016793, 1.9733635, 1.3079423, 1.0791487, 0.6343523
  )
)

prevalence_map_df <- data.frame(
  location_name = c(
    "Australia and New Zealand",
    "Europe & Central Asia",
    "North America",
    "South Asia",
    "East Asia & Pacific",
    "Middle East & North Africa",
    "Latin America & Caribbean",
    "Sub-Saharan Africa",
    "Oceania"
  ),
  prevalence = c(
    269.268860, 198.090941, 198.090941, 43.507169,
    9.795618, 35.192616, 15.253657, 8.713955, 4.510866
  )
)

incidence_bar_df <- data.frame(
  region = c(
    "Australia & New Zealand",
    "Europe & Northern America",
    "Central & Southern Asia",
    "Northern Africa & Western Asia",
    "Latin America & Caribbean",
    "Eastern & South-Eastern Asia",
    "Sub-Saharan Africa",
    "Oceania"
  ),
  incidence = c(
    23.8575598, 16.3063269, 5.2794844, 3.3016793,
    1.9733635, 1.3079423, 1.0791487, 0.6343523
  )
)

prevalence_bar_df <- data.frame(
  region = c(
    "Australia & New Zealand",
    "Europe & Northern America",
    "Central & Southern Asia",
    "Northern Africa & Western Asia",
    "Latin America & Caribbean",
    "Eastern & South-Eastern Asia",
    "Sub-Saharan Africa",
    "Oceania"
  ),
  prevalence = c(
    269.268860, 198.090941, 43.507169, 35.192616,
    15.253657, 9.795618, 8.713955, 4.510866
  )
)

build_world_panel <- function(plot_df, value_col, title, legend_title, breaks) {
  ggplot(plot_df) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.15) +
    scale_fill_gradientn(
      colours = c("#deebf7", "#9ecae1", "#6baed6", "#3182bd", "#08306b"),
      breaks = breaks,
      limits = range(breaks),
      name = legend_title,
      na.value = "grey92"
    ) +
    labs(title = title, subtitle = "Data Source: GBD Estimates") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 7)
    )
}

build_region_bar <- function(df, xcol, title, fill_col, xlab) {
  df$region <- factor(df$region, levels = rev(df$region))
  ggplot(df, aes(x = .data[[xcol]], y = region)) +
    geom_col(fill = fill_col) +
    theme_minimal(base_size = 10) +
    labs(x = xlab, y = "Region", title = title) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 6)
    )
}

toplot <- rnaturalearth::ne_countries(scale = "medium", type = "map_units", returnclass = "sf") %>%
  dplyr::left_join(regions, by = c("iso_a2"))
toplot$iso_a2 <- toplot$iso_a2_eh
toplot <- toplot %>%
  dplyr::filter(iso_a2 %in% regions$iso_a2)

toplot_stage <- toplot %>%
  dplyr::left_join(country_stage, by = c("iso_a2" = "iso"))

color_df <- data.frame(
  ISO = names(expanded_colors),
  Color = unname(expanded_colors)
) %>%
  dplyr::mutate(
    Row = rep(8:1, each = 8),
    Col = rep(1:8, 8)
  )

stage_meta <- meta %>%
  dplyr::left_join(country_stage, by = c("iso" = "iso")) %>%
  dplyr::filter(!is.na(stage)) %>%
  dplyr::mutate(stage = factor(stage, levels = c("1", "2", "3")))

working <- as.data.frame(filtered)
working$sample <- rownames(working)
working_stage <- working %>%
  dplyr::inner_join(stage_meta %>% dplyr::select(sample, stage), by = "sample")
rownames(working_stage) <- working_stage$sample
abundance_only <- working_stage %>% dplyr::select(-sample, -stage)

diversity_df <- data.frame(
  sample = rownames(working_stage),
  stage = working_stage$stage,
  shannon = vegan::diversity(abundance_only, index = "shannon"),
  simpson = vegan::diversity(abundance_only, index = "simpson"),
  species_count = vegan::specnumber(abundance_only)
)

country_meta <- meta %>%
  dplyr::filter(!is.na(iso), iso %in% names(expanded_colors)) %>%
  dplyr::select(sample, iso, region)

country_counts_df <- country_meta %>%
  dplyr::count(iso, name = "Freq") %>%
  dplyr::arrange(dplyr::desc(Freq))

working_country <- working %>%
  dplyr::inner_join(country_meta, by = "sample") %>%
  dplyr::filter(!region %in% c("unknown", "", "Oceania"))
rownames(working_country) <- working_country$sample

country_abundance <- working_country %>% dplyr::select(-sample, -iso, -region)
country_diversity_df <- data.frame(
  sample = rownames(working_country),
  iso = working_country$iso,
  shannon = vegan::diversity(country_abundance, index = "shannon"),
  simpson = vegan::diversity(country_abundance, index = "simpson"),
  species_count = vegan::specnumber(country_abundance)
) %>%
  dplyr::filter(!iso %in% c("HK", "BW", "TW"))

stage_heat <- function(stage_id, xvar = "mds1", yvar = "mds2", xlab = "", ylab = "") {
  toplot_local <- points.annotated %>%
    dplyr::filter(stage == stage_id) %>%
    dplyr::select(dplyr::all_of(c(xvar, yvar)))
  background <- points.annotated %>%
    dplyr::select(dplyr::all_of(c(xvar, yvar)))
  colnames(toplot_local) <- c("v1", "v2")
  colnames(background) <- c("v1", "v2")
  ggplot(toplot_local, aes(x = v1, y = v2)) +
    geom_bin_2d(data = background, bins = 30) +
    scale_colour_gradient(low = "#d9d9d9", high = "#d9d9d9", aesthetics = "fill") +
    ggnewscale::new_scale_fill() +
    geom_bin_2d(data = toplot_local, bins = 30) +
    scale_fill_viridis_c(option = "D") +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme_bw() +
    labs(x = xlab, y = ylab, title = paste("Stage", stage_id)) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold")
    )
}

mds_density_panel <- function(mds_col, show_y = FALSE) {
  toplot_local <- points.annotated %>%
    dplyr::select(stage, dplyr::all_of(mds_col)) %>%
    dplyr::rename(selected = dplyr::all_of(mds_col)) %>%
    dplyr::filter(!is.na(stage))
  toplot_local$stage <- factor(toplot_local$stage, levels = c("3", "2", "1"))
  p <- ggplot(toplot_local, aes(x = selected, y = stage, fill = stage)) +
    ggridges::stat_density_ridges(
      quantile_lines = TRUE,
      quantiles = 2,
      scale = 1.4,
      alpha = 0.95,
      color = "black",
      linewidth = 0.25
    ) +
    scale_fill_manual(values = stage_colors) +
    theme_bw() +
    labs(x = toupper(mds_col), y = if (show_y) "Stage" else "")
  if (show_y) {
    p + theme(legend.position = "none")
  } else {
    p + theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
}

#Figure 1A
plot_a_df <- dplyr::left_join(world, prevalence_map_df, by = c("region_wb" = "location_name"))
p_a <- build_world_panel(
  plot_df = plot_a_df,
  value_col = "prevalence",
  title = "Global Distribution of IBD Prevalence (by Region)",
  legend_title = "IBD Prevalence Rate",
  breaks = c(0, 50, 100, 150, 200, 250)
)
p_a
#Figure 1B
plot_b_df <- dplyr::left_join(world, incidence_map_df, by = c("region_wb" = "location_name"))
p_b <- build_world_panel(
  plot_df = plot_b_df,
  value_col = "incidence",
  title = "Global Distribution of IBD Incidence Rate (by Region)",
  legend_title = "IBD Incidence Rate",
  breaks = c(0, 5, 10, 15, 20)
)
p_b
#Figure 1C
p_c <- build_region_bar(
  incidence_bar_df,
  xcol = "incidence",
  title = "Global IBD Incidence by Region",
  fill_col = "#B2DF8A",
  xlab = "IBD Incidence"
)
p_c
#Figure 1D
p_d <- build_region_bar(
  prevalence_bar_df,
  xcol = "prevalence",
  title = "Global IBD Prevalence by Region",
  fill_col = "#A6CEE3",
  xlab = "IBD Prevalence"
)
p_d
#Figure 1E
p_e_map <- ggplot(data = toplot) +
  geom_sf(aes(fill = iso_a2), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = c(expanded_colors, Other = "gray90")) +
  coord_sf(crs = "+proj=eqearth +wktext") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
p_e_legend <- ggplot(color_df, aes(x = Col, y = Row, fill = Color)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = ISO), color = "black", size = 2.1, fontface = "bold") +
  scale_fill_identity() +
  coord_fixed() +
  theme_void()
p_e <- p_e_map + p_e_legend + patchwork::plot_layout(widths = c(2.6, 1.2))
p_e
#Figure 1F
p_f <- ggplot(data = toplot_stage) +
  geom_sf(aes(fill = factor(stage)), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = stage_colors, name = NULL) +
  coord_sf(crs = "+proj=eqearth +wktext") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm")
  )
p_f
#Figure 1G
stage_counts <- stage_meta %>%
  dplyr::count(stage, name = "samples") %>%
  dplyr::mutate(stage = factor(stage, levels = c("1", "2", "3")))
p_g <- ggplot(stage_counts, aes(x = samples, y = stage, fill = stage)) +
  geom_col(alpha = 0.55) +
  scale_fill_manual(values = stage_colors) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  labs(x = "Samples", y = "Stage") +
  theme(legend.position = "none")
p_g
#Figure 1H
build_stage_violin <- function(df, value_col, xlab, show_y = FALSE) {
  p <- ggplot(df, aes(x = .data[[value_col]], y = stage, fill = stage)) +
    geom_violin(trim = FALSE, linewidth = 0.25) +
    stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
    scale_fill_manual(values = stage_colors) +
    theme_bw() +
    labs(x = xlab, y = if (show_y) "Stage" else "")
  if (show_y) {
    p + theme(legend.position = "none")
  } else {
    p + theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
}
p_h_shannon <- build_stage_violin(diversity_df, "shannon", "Shannon diversity", TRUE)
p_h_simpson <- build_stage_violin(diversity_df, "simpson", "Simpson diversity", FALSE)
p_h_species <- build_stage_violin(diversity_df, "species_count", "Species count", FALSE)
p_h <- p_h_shannon + p_h_simpson + p_h_species + patchwork::plot_layout(ncol = 3)
p_h
#Figure 1I
p_i_scatter <- ggplot(points.annotated, aes(x = mds1, y = mds2, color = stage)) +
  geom_point(size = 0.08, alpha = 0.32) +
  scale_color_manual(values = stage_colors, name = "Stage") +
  scale_x_continuous(limits = xlims, expand = c(0, 0)) +
  scale_y_continuous(limits = ylims, expand = c(0, 0)) +
  theme_bw() +
  labs(
    x = "",
    y = paste0("MDS2 (", round(nmds$eigen[2] / eigen_total, 1), "%)"),
    title = paste0("MDS1 (", round(nmds$eigen[1] / eigen_total, 1), "%)")
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
p_i_stage1 <- stage_heat("1")
p_i_stage2 <- stage_heat("2")
p_i_stage3 <- stage_heat(
  "3",
  xlab = paste0("MDS1 (", round(nmds$eigen[1] / eigen_total, 1), "%)"),
  ylab = "MDS2"
)
p_i <- p_i_scatter + p_i_stage1 + p_i_stage2 + p_i_stage3 +
  patchwork::plot_layout(ncol = 4, widths = c(1.05, 1, 1, 1))

#Figure 1J
p_j1 <- mds_density_panel("mds1", TRUE)
p_j2 <- mds_density_panel("mds2", FALSE)
p_j3 <- mds_density_panel("mds3", FALSE)
p_j4 <- mds_density_panel("mds4", FALSE)
p_j <- p_j1 + p_j2 + p_j3 + p_j4 + patchwork::plot_layout(ncol = 4)
p_j

#Figure 22
samplecount_plot_iso_df <- country_counts_df
samplecount_plot_iso_df$iso <- factor(
  samplecount_plot_iso_df$iso,
  levels = rev(samplecount_plot_iso_df$iso)
)
p_s22 <- ggplot(samplecount_plot_iso_df, aes(y = iso, x = Freq, fill = iso)) +
  geom_col(alpha = 0.5) +
  geom_text(
    aes(x = max(Freq) * 0.01, y = iso, label = iso),
    hjust = 0,
    size = 3.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = c(0, 20000, 50000, 100000),
    limits = c(0, max(samplecount_plot_iso_df$Freq) * 1.05),
    expand = c(0, 5000),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  scale_fill_manual(values = expanded_colors, name = "Country") +
  theme_bw() +
  labs(x = "Samples", y = "Country") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    rect = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
p_s22
#Figure 23
iso_order_samplecount <- rev(country_counts_df$iso)
country_diversity_df$iso <- factor(country_diversity_df$iso, levels = iso_order_samplecount)
p_s23 <- ggplot(country_diversity_df, aes(x = shannon, y = iso, fill = iso)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
  scale_fill_manual(values = expanded_colors, name = "Country") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Shannon diversity")
p_s23
#Figure 24
p_s24 <- ggplot(country_diversity_df, aes(x = simpson, y = iso, fill = iso)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
  scale_fill_manual(values = expanded_colors, name = "Country") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Simpson diversity")
p_s24
#Figure 25
p_s25 <- ggplot(country_diversity_df, aes(x = species_count, y = iso, fill = iso)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
  scale_fill_manual(values = expanded_colors, name = "Country") +
  coord_cartesian(xlim = c(0, 180)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Species count")
p_s25
#Figure 1C-F composite
p_country_diversity_legacy <- p_s23 / p_s24 / p_s25
figure1_design <- "
aaabbb..
aaabbb..
cccddd..
eeeeefff
gghhhhhh
iiiiiiii
iiiiiiii
jjjjjjjj
"

p_figure1 <- patchwork::wrap_plots(
  a = p_a,
  b = p_b,
  c = p_c,
  d = p_d,
  e = p_e,
  f = p_f,
  g = p_g,
  h = p_h,
  i = p_i,
  j = p_j,
  design = figure1_design
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  ) &
  theme(
    plot.margin = margin(4, 4, 4, 4)
  )

save_panel(p_a, "Figure1A_global_prevalence_map.pdf", 6.2, 3.2)
save_panel(p_b, "Figure1B_global_incidence_map.pdf", 6.2, 3.2)
save_panel(p_c, "Figure1C_global_incidence_bar.pdf", 4.2, 3.2)
save_panel(p_d, "Figure1D_global_prevalence_bar.pdf", 4.2, 3.2)
save_panel(p_e, "Figure1E_country_classification_map.pdf", 9.5, 3.8)
save_panel(p_f, "Figure1F_stage_map.pdf", 6.2, 3.6)
save_panel(p_g, "Figure1G_stage_samplecount.pdf", 3.8, 2.8)
save_panel(p_h, "Figure1H_stage_diversity_violin.pdf", 10.5, 3.1)
save_panel(p_i, "Figure1I_stage_pcoa.pdf", 13.5, 3.6)
save_panel(p_j, "Figure1J_stage_mds_density.pdf", 13.5, 3.8)
save_panel(p_figure1, "Figure1_A-J_composite.pdf", 14, 17.5)
save_panel(p_s22, "Supplementary_Figure22_country_samplecount.pdf", 7, 8)
save_panel(p_s23, "Supplementary_Figure23_country_shannon_violin.pdf", 9.6, 10)
save_panel(p_s24, "Supplementary_Figure24_country_simpson_violin.pdf", 9.6, 10)
save_panel(p_s25, "Supplementary_Figure25_country_speciescount_violin.pdf", 9.6, 10)
save_panel(p_s22, "Figure1C_samplecount_iso.pdf", 7, 8)
save_panel(p_country_diversity_legacy, "Figure1D-F_diversity_iso.pdf", 9.6, 30.2)
