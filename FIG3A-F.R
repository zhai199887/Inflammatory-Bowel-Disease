
#Figure 3: Taxa are differentially abundant between IBD stage-associated regions
#IBD

#Output panels:
#Fig3A - Stage comparison bubble chart
#Fig3B - Stage p-value heatmap + mean heatmap + bar chart
#Fig3C - Ridge plots for top 10 genera across 3 stages
#Fig3D - LEfSe: Stage 1 vs Stage 2 bar chartLEfSe 1 vs 2
#Fig3E - LEfSe: Stage 2 vs Stage 3 bar chartLEfSe 2 vs 3
#Fig3F - LEfSe: Stage 1 vs Stage 3 bar chartLEfSe 1 vs 3

#Supplementary figures:
#SuppFig_country_bubble - Country-level bubble chart
#SuppFig_country_heatmap - Country-level p-value + mean heatmap
#SuppFig_country_ridges - Country-level ridge plots
#SuppFig_LEfSe_cladogram - LEfSe cladogram all 3 stages
#SuppFig_LEfSe_bar_all3 - LEfSe bar all 3 stages

#Supplementary tables:
#Supplementary_Table_LEfSe_Stage1vs2_results.csv
#Supplementary_Table_LEfSe_Stage2vs3_results.csv
#Supplementary_Table_LEfSe_Stage1vs3_results.csv



#Package loading


`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

pkgs <- c("dplyr","readr","stringr","tibble","tidyr","ggplot2",
          "patchwork","cowplot","ggridges","viridis","colorspace",
          "scales","grid","microeco","magrittr")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "))
invisible(lapply(pkgs, library, character.only = TRUE))


#Path setup


script_dir <- tryCatch({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
}, error = function(e) getwd())

output_dir <- script_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Data search directories (priority order)
data_dirs <- unique(c(
  script_dir,
  dirname(script_dir),
  "D:/R\u4ee3\u7801",
  "D:/IHME-GBD_2021_DATA-29fb49f4-1/\u6587\u7ae0/\u65b0fig/table 1-44/table 1-44"
))

find_file <- function(filename, extra_dirs = character()) {
  dirs <- unique(c(data_dirs, extra_dirs))
  hit  <- Filter(file.exists, file.path(dirs, filename))[1]
  if (is.na(hit)) stop("File not found: ", filename)
  hit
}

read_csv_safely <- function(filename, ...) {
  readr::read_csv(find_file(filename), show_col_types = FALSE, ...)
}

message("Output directory: ", output_dir)


#Color definitions


#CSV"stage1" / "stage2" / "stage3"
#Stage colors with "stage" prefix, matching bubble chart CSV labels
stage_colors <- c("stage1" = "#1B9E77", "stage2" = "#D95F02", "stage3" = "#7570B3")

#LEfSe"1" / "2" / "3"
#Stage colors short keys "1","2","3" for heatmaps, ridge plots, LEfSe bars
sc <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

iso_list <- c(
  "JO","SD","AZ","EE","MZ","GT","MX","CM","TZ","ML","RO","NG","KZ","EC","TH","RS","VE",
  "RU","CF","MG","MW","PE","SG","BR","PL","TR","LT","HR","PT","UG","AE","IN","IT","ZW",
  "AU","IE","DK","ES","BD","CO","KE","HU","NL","MY","NO","GR","AM","CH","CA","SE","AT",
  "IS","FR","GB","NZ","FI","DE","IL","BE","CN","KR","US","GH","JP"
)
expanded_colors <- setNames(
  colorspace::qualitative_hcl(length(iso_list), palette = "Dark 3",
                               h = c(0, 360), c = 60, l = 65),
  iso_list
)


#Data loading


message("Loading taxa names, genus matrix and metadata...")

taxa_names <- read_csv_safely("taxaNames.csv") %>%
  dplyr::distinct(full.taxon, taxon)

genus_raw <- readRDS(find_file("genus.rds"))
if (!is.data.frame(genus_raw)) genus_raw <- as.data.frame(genus_raw)

metadata_raw <- readRDS(find_file("metadata_for_diffAbundance.rds"))
metadata <- tibble::as_tibble(metadata_raw)

#Normalize column names
if (!"sample" %in% names(metadata) && all(c("project","srr") %in% names(metadata)))
  metadata <- dplyr::mutate(metadata, sample = paste(project, srr, sep = "_"))
for (col in c("iso","stage")) {
  variants <- paste0(col, c("", ".x", ".y"))
  found    <- variants[variants %in% names(metadata)][1]
  if (!is.na(found) && found != col)
    metadata <- dplyr::mutate(metadata, !!col := .data[[found]])
}
metadata <- metadata %>%
  dplyr::filter(!is.na(sample)) %>%
  dplyr::distinct(sample, .keep_all = TRUE)

#Relative abundance
make_rel <- function(m) { rs <- rowSums(m, na.rm=TRUE); rs[rs==0] <- 1; sweep(m,1,rs,"/") }

genus_tbl    <- tibble::as_tibble(genus_raw, rownames = "sample")
common_samps <- intersect(genus_tbl$sample, metadata$sample)
genus_tbl    <- dplyr::filter(genus_tbl, sample %in% common_samps)
metadata     <- dplyr::filter(metadata,   sample %in% common_samps) %>%
  dplyr::arrange(match(sample, genus_tbl$sample))
genus_tbl    <- dplyr::arrange(genus_tbl, match(sample, metadata$sample))
stopifnot(identical(genus_tbl$sample, metadata$sample))

genus_mat <- as.data.frame(dplyr::select(genus_tbl, -sample), check.names = FALSE)
row.names(genus_mat) <- metadata$sample
genus_mat[genus_mat == 0] <- 1
genus_rel <- make_rel(genus_mat)

#Rename taxon columns
ord          <- match(colnames(genus_rel), taxa_names$full.taxon)
keep         <- !is.na(ord)
genus_rel    <- genus_rel[, keep, drop = FALSE]
colnames(genus_rel) <- taxa_names$taxon[ord[keep]]

rename_taxa_df <- function(df, col = "taxon") {
  for (i in seq_len(nrow(taxa_names)))
    df[[col]][df[[col]] == taxa_names$full.taxon[i]] <- taxa_names$taxon[i]
  df
}

message("Data loaded: ", nrow(genus_rel), " samples, ", ncol(genus_rel), " genera")


#Shared plotting functions


#Bubble chart --------------------------------------------------
#fig4A-C.Rscale_radius
#Recreates the original single-plot layout with diagonal colored labels and scale_radius
build_overview_bubble <- function(data_file, fill_values, output_file,
                                  label_levels, bubble_size = 3,
                                  height_in = 5, width_in = 5) {
  taxa_counts <- read_csv_safely(data_file)

#fig4A-C.R
#Derive region order from data (matches original fig4A-C.R) so all pairs form lower triangle
  plot_region_order <- unique(c(as.character(taxa_counts$region1),
                                 as.character(taxa_counts$region2)))

#Diagonal label data rev
  labels <- tibble::tibble(
    region2     = rev(plot_region_order),
    regionLabel = rev(plot_region_order)
  )

#taxa_counts regionLabel
  if ("regionLabel" %in% names(taxa_counts))
    taxa_counts <- dplyr::select(taxa_counts, -regionLabel)

  num_taxa <- max(taxa_counts$freq, na.rm = TRUE)

#Bubble chart (without legend, to be inset)
  bubble <- ggplot() +
    geom_point(
      data = taxa_counts,
      aes(x = factor(region1, levels = rev(plot_region_order)),
          y = factor(region2, levels = rev(plot_region_order)),
          size = freq),
      fill = "grey", shape = 21, stroke = 0
    ) +
    geom_text(
      data = taxa_counts,
      aes(x = factor(region1, levels = rev(plot_region_order)),
          y = factor(region2, levels = rev(plot_region_order)),
          label = freq),
      size = bubble_size, vjust = 0
    ) +
    geom_label(
      data = labels,
      aes(x = factor(region2, levels = rev(plot_region_order)),
          y = factor(region2, levels = rev(plot_region_order)),
          label = regionLabel,
          fill  = factor(regionLabel, levels = rev(plot_region_order))),
      size = bubble_size
    ) +
    scale_fill_manual(values = fill_values, aesthetics = c("colour","fill"), drop = FALSE) +
    scale_radius(
      range  = c(2, 22),
      limits = c(0, max(70, num_taxa)),
      breaks = unique(sort(c(10, 20, 30, 50, num_taxa)))
    ) +
    labs(size = "Number of\nDifferentially\nAbundant Taxa", x = NULL, y = NULL) +
    guides(
      size  = guide_legend(override.aes = list(fill="grey", stroke=0.25),
                           label.position = "bottom", title.position = "top", order = 1),
      fill  = "none",
      color = "none"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      plot.margin     = grid::unit(c(0,0,0,0),"pt")
    )


#Extract legend and inset into bubble panel (matches original)
  size_legend <- cowplot::get_legend(bubble)
  bubble <- bubble + theme(legend.position = "none") +
    patchwork::inset_element(size_legend, 0.15, 0.6, 0.4, 0.95,
                              align_to = "panel", ignore_tag = TRUE)

#wrap_elements(bubble + label + plot_layout(heights=c(39,1)))
#Bottom color strip combined with bubble using plot_layout (matches original)
  lp <- ggplot(labels) +
    geom_tile(aes(x = factor(region2, levels = rev(plot_region_order)), y = 1,
                  fill = factor(regionLabel, levels = rev(plot_region_order)))) +
    scale_fill_manual(values = fill_values, aesthetics = c("colour","fill"), drop = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = grid::unit(c(0,0,0,0),"pt")) +
    scale_x_discrete(expand = c(0,0), drop = FALSE) +
    scale_y_discrete(expand = c(0,0))

  overview <- patchwork::wrap_elements(bubble + lp + patchwork::plot_layout(heights = c(39,1)))

  ggsave(output_file, plot = overview, device = "pdf",
         height = height_in, width = width_in, units = "in")
  message("Saved: ", basename(output_file))
  invisible(overview)
}

#P / P-value heatmap --------------------------------------------
#inset_element
#Bottom color strip overlaid via inset_element, matching original code style
build_pval_heatmap <- function(pvals, taxa_means, x_var, x_levels, fill_values) {
  pvals_plot <- dplyr::filter(pvals, !is.na(.data[[x_var]]), !is.na(taxon))
  pvals_plot[[x_var]] <- factor(as.character(pvals_plot[[x_var]]), levels = x_levels)
  hm <- ggplot(pvals_plot) +
    geom_tile(aes(x = .data[[x_var]], y = factor(taxon, levels = rev(taxa_means$taxon)),
                  fill = log10(p.adj))) +
    scale_fill_gradient2(low = "red", high = "white") +
    theme(axis.text.x = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), text = element_text(size=8),
          legend.position = "top", legend.title = element_text(size=9)) +
    scale_x_discrete(expand = c(0,0), drop = FALSE) +
    scale_y_discrete(expand = c(0,0)) +
    guides(fill = guide_colorbar(title="Adjusted P Value",
                                  title.position="top", title.hjust=0.5))
#inset_element
#Color strip: inset_element at bottom of heatmap panel (matches original)
  ldf <- tibble::tibble(x_ = factor(x_levels, levels = x_levels))
  names(ldf) <- x_var
  lp <- ggplot(ldf) +
    geom_tile(aes(x = .data[[x_var]], y = 1, fill = .data[[x_var]])) +
    scale_fill_manual(values = fill_values, aesthetics = c("colour","fill"), drop = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = grid::unit(c(0,0,0,0),"pt")) +
    scale_x_discrete(expand = c(0,0), drop = FALSE) +
    scale_y_discrete(expand = c(0,0))
  combined <- hm + patchwork::inset_element(lp, 0, -0.02, 1, 0,
                                             align_to = "panel", ignore_tag = TRUE)
  list(plot = hm, legend = cowplot::get_legend(hm), combined = combined)
}

#Mean abundance heatmap -------------------------------------
build_mean_heatmap <- function(mean_df, taxa_means, x_var, x_levels, fill_values) {
  mean_df[[x_var]] <- factor(as.character(mean_df[[x_var]]), levels = x_levels)
  hm <- ggplot(mean_df) +
    geom_tile(aes(x = .data[[x_var]], y = factor(taxon, levels = rev(taxa_means$taxon)),
                  fill = log10(mean))) +
    viridis::scale_fill_viridis(option = "D") +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position = "top",
          legend.text = element_text(size=8), legend.title = element_text(size=9),
          plot.margin = grid::unit(c(0,0,0,0),"pt")) +
    scale_x_discrete(expand = c(0,0), drop = FALSE) +
    scale_y_discrete(expand = c(0,0)) +
    guides(fill = guide_colorbar(title="Mean Relative Abundance",
                                  title.position="top", title.hjust=0.5))
  lp <- ggplot(dplyr::distinct(mean_df, .data[[x_var]])) +
    geom_tile(aes(x = .data[[x_var]], y = 1, fill = .data[[x_var]])) +
    scale_fill_manual(values = fill_values, aesthetics = c("colour","fill"), drop = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = grid::unit(c(0,0,0,0),"pt")) +
    scale_x_discrete(expand = c(0,0), drop = FALSE) +
    scale_y_discrete(expand = c(0,0))
  combined <- hm + patchwork::inset_element(lp, 0, -0.02, 1, 0,
                                             align_to = "panel", ignore_tag = TRUE)
  list(plot = hm, legend = cowplot::get_legend(hm), combined = combined)
}

#Mean bar chart --------------------------------------------
build_mean_bar <- function(taxa_means) {
  ggplot(taxa_means) +
    geom_col(aes(x = taxa_means, y = factor(taxon, levels = rev(taxon))), width = 0.6) +
    xlab("Mean Relative\nAbundance") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_text(size=7),
          axis.text.x = element_text(size=6),
          panel.spacing = grid::unit(c(0,0,0,0),"pt"),
          plot.margin   = grid::unit(c(0,0,0,0),"pt")) +
    scale_x_continuous(breaks = c(0, 0.1))
}

#Ridge plot ---------------------------------------------------
build_ridges <- function(genus_long, group_var, group_levels, fill_values) {
  genus_long[[group_var]] <- factor(as.character(genus_long[[group_var]]), levels = group_levels)
  ggplot(genus_long) +
    ggridges::geom_density_ridges(
      aes(x = log10(rel_abundance), y = .data[[group_var]], fill = .data[[group_var]]),
      quantile_lines = TRUE, quantiles = 2
    ) +
    scale_fill_manual(values = fill_values, aesthetics = c("colour","fill"), drop = FALSE) +
    labs(x = "Relative Abundance", y = "Density") +
    theme(legend.position = "none",
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(color = "black"),
          plot.title   = element_text(size = 8),
          plot.margin  = grid::unit(c(10,10,10,10),"pt")) +
    scale_x_continuous(labels = scales::math_format(), breaks = c(-5,-1)) +
    facet_wrap(~taxon, ncol = 5)
}

#LEfSe / Custom LEfSe bar with stage colors ---
#trans_diff
#Extracts trans_diff$res_diff and draws colored horizontal bar chart
plot_lefse_custom <- function(lefse_obj, group1, group2, colors,
                               top_n = 30, title = NULL) {
  res <- lefse_obj$res_diff
  if (is.null(res) || nrow(res) == 0) {
    message("No LEfSe results found"); return(ggplot() + theme_void())
  }

#LDA / Find LDA score column
  score_col <- intersect(c("Score","LDA_score","lda","LDA"), colnames(res))[1]
  if (is.na(score_col)) {
    message("No score column found in res_diff. Columns: ", paste(colnames(res), collapse=", "))
    return(ggplot() + theme_void())
  }

  res <- res[!is.na(res[[score_col]]), ]
  res <- res[order(abs(res[[score_col]]), decreasing = TRUE), ]
  res <- head(res, top_n)

#Find taxon name column
  taxon_col <- intersect(c("Taxon","Taxa","taxon","taxa","feature"), colnames(res))[1]

#group1=group2
#Positive = enriched in group1, negative = enriched in group2
  res$display_score <- ifelse(res$Group == group1, res[[score_col]], -res[[score_col]])
  res$fill_group    <- as.character(res$Group)

#Order taxa by score
  res[[taxon_col]] <- factor(res[[taxon_col]],
                              levels = res[[taxon_col]][order(res$display_score)])

  p <- ggplot(res) +
    geom_col(aes(x = display_score, y = .data[[taxon_col]], fill = fill_group), width = 0.8) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(
      values = colors,
      name   = "Group",
      labels = setNames(paste0("Stage ", names(colors)), names(colors))
    ) +
    labs(x = "LDA Score (log10)", y = NULL, title = title) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.y     = element_text(size = 7),
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 10, hjust = 0.5)
    )
  p
}


#Data preparation


#Stage-level data -------------------------------------------
#sc
#Uses short format "1","2","3" to match sc color keys
prepare_stage_data <- function() {
  pvals <- read_csv_safely("Supplementary_Table_13_DA_LMM_emmeans_results_stage.csv") %>%
    rename_taxa_df() %>%
    dplyr::mutate(
      region1 = gsub("^stage", "", as.character(region1)),
      region2 = gsub("^stage", "", as.character(region2))
    ) %>%
    dplyr::filter(!is.na(taxon), !is.na(region2), !is.na(p.adj)) %>%
    dplyr::filter(region2 == "3")

#stage
  stage_meta <- metadata %>%
    dplyr::transmute(sample, stage = as.character(stage)) %>%
    dplyr::filter(stage %in% c("1","2","3"))

  stage_genus <- genus_rel[stage_meta$sample, , drop = FALSE]
  stage_meta  <- dplyr::arrange(stage_meta, match(sample, rownames(stage_genus)))

  mean_df <- as.data.frame(stage_genus) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(stage_meta, by = "sample") %>%
    tidyr::pivot_longer(-c(sample, stage), names_to = "taxon", values_to = "mean") %>%
    dplyr::group_by(stage, taxon) %>%
    dplyr::summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")

  taxa_means <- tibble::tibble(
    taxon      = colnames(stage_genus),
    taxa_means = colMeans(stage_genus, na.rm = TRUE)
  ) %>% dplyr::arrange(dplyr::desc(taxa_means))

  top10 <- taxa_means$taxon[1:min(10, nrow(taxa_means))]
  genus_long <- as.data.frame(stage_genus) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(stage_meta, by = "sample") %>%
    tidyr::pivot_longer(-c(sample, stage), names_to = "taxon", values_to = "rel_abundance") %>%
    dplyr::filter(taxon %in% top10) %>%
    dplyr::mutate(taxon = factor(taxon, levels = top10))

  list(pvals = pvals, mean_df = mean_df, taxa_means = taxa_means, genus_long = genus_long)
}

#ISO-level data ---------------------------------------------
prepare_iso_data <- function() {
  pvals <- read_csv_safely("Supplementary_Table_14_DA_LMM_emmeans_results_iso.csv") %>%
    rename_taxa_df() %>%
    dplyr::filter(!is.na(taxon), !is.na(region2), !is.na(p.adj))

  iso_meta <- metadata %>%
    dplyr::transmute(sample, iso = as.character(iso)) %>%
    dplyr::filter(iso %in% iso_list)

  iso_genus <- genus_rel[iso_meta$sample, , drop = FALSE]
  iso_meta  <- dplyr::arrange(iso_meta, match(sample, rownames(iso_genus)))

  mean_df <- as.data.frame(iso_genus) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(iso_meta, by = "sample") %>%
    tidyr::pivot_longer(-c(sample, iso), names_to = "taxon", values_to = "mean") %>%
    dplyr::group_by(iso, taxon) %>%
    dplyr::summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")

  taxa_means <- tibble::tibble(
    taxon      = colnames(iso_genus),
    taxa_means = colMeans(iso_genus, na.rm = TRUE)
  ) %>% dplyr::arrange(dplyr::desc(taxa_means))

  genus_long <- as.data.frame(iso_genus) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(iso_meta, by = "sample") %>%
    tidyr::pivot_longer(-c(sample, iso), names_to = "taxon", values_to = "rel_abundance") %>%
    dplyr::filter(taxon %in% taxa_means$taxon[1:min(10, nrow(taxa_means))])

  list(pvals = pvals, mean_df = mean_df, taxa_means = taxa_means, genus_long = genus_long)
}


#Fig3A - / Stage comparison bubble chart

#IBDstage1 / 2 / 3

#Shows differentially abundant taxa counts between pairs of IBD stages.
#Bubble size / number = count of significantly different taxa; diagonal = colored stage labels.

message("\n=== Fig3A: Stage bubble chart ===")
tryCatch({
  p_A <- build_overview_bubble(
    data_file    = "diff_taxa_counts_for_my5AB.csv",
    fill_values  = stage_colors,
    output_file  = file.path(output_dir, "Figure3A_stage_bubble.pdf"),
    label_levels = names(stage_colors),#stage1","stage2","stage3"
    bubble_size  = 4,
    height_in    = 5,
    width_in     = 5
  )
}, error = function(e) message("ERROR Fig3A: ", e$message))
p_A

#Fig3B - / Stage heatmap (p-value + mean abundance)

#Plog10stage2stage3stage1
#log10stage1 / 2 / 3

#Left heatmap: adjusted p-value (log10), columns = stage2 & stage3 vs stage1
#Right heatmap: mean relative abundance (log10), 3 columns stage1 / 2 / 3
#Far right: mean bar chart

message("\n=== Fig3B: Stage heatmap ===")
stage_data <- NULL
tryCatch({
  stage_data <- prepare_stage_data()

#Pregion2=="3"xregion1"1""2"stage1stage2
#x_levels for p-value heatmap: region1 values "1","2" (stage1, stage2 vs stage3)
  pval_x_levels <- sort(unique(stage_data$pvals$region1[!is.na(stage_data$pvals$region1)]))

  stage_mean_bar <- build_mean_bar(stage_data$taxa_means)
  stage_pval     <- build_pval_heatmap(stage_data$pvals, stage_data$taxa_means,
                                        "region1", pval_x_levels, sc)
  stage_mean     <- build_mean_heatmap(stage_data$mean_df, stage_data$taxa_means,
                                        "stage", c("1","2","3"), sc)

#fig4A-C.Rwrap_elements(p + m + mean.bar)inset
#Layout matching original: wrap_elements(p + m + mean.bar), legends as inset above
  bars_B <- patchwork::wrap_elements(
    stage_pval$combined + stage_mean$combined + stage_mean_bar
  )
  p_B <- bars_B +
    patchwork::inset_element(stage_pval$legend, 0.15, 1.02, 0.59, 1.07,
                              align_to = "plot", ignore_tag = TRUE) +
    patchwork::inset_element(stage_mean$legend, 0.60, 1.02, 1.00, 1.07,
                              align_to = "plot", ignore_tag = TRUE)

  ggsave(file.path(output_dir, "Figure3B_stage_heatmap.pdf"),
         plot = p_B, width = 14, height = 8, dpi = 300)
  message("Saved: Figure3B_stage_heatmap.pdf")
}, error = function(e) message("ERROR Fig3B: ", e$message))
p_B

#Fig3C - / Ridge plots for top 10 genera across 3 stages

#log10
#y "1","2","3"
#Each facet: log10 relative abundance density for one genus across 3 stages.
#Y-axis labels "1","2","3" with stage colors.

message("\n=== Fig3C: Stage ridge plots ===")
tryCatch({
  if (is.null(stage_data)) stage_data <- prepare_stage_data()
#stage "1","2","3" / stage column already "1","2","3"
  p_C <- build_ridges(stage_data$genus_long, "stage", c("1","2","3"), sc)
  ggsave(file.path(output_dir, "Figure3C_stage_ridges.pdf"),
         plot = p_C, width = 14, height = 6, dpi = 300)
  message("Saved: Figure3C_stage_ridges.pdf")
}, error = function(e) message("ERROR Fig3C: ", e$message))
p_C

#LEfSe / LEfSe data preparation

#microeco otu / group / tax LEfSe
#IBD
#group "1","2","3"
#Uses microeco; group column values are "1","2","3".

message("\n=== LEfSe data loading ===")
lefse_ready <- FALSE
dataset     <- NULL
tryCatch({
  otu   <- read.csv(find_file("otu.csv"),   row.names = 1)
  grp   <- read.csv(find_file("group.csv"), row.names = 1)
  tax   <- read.csv(find_file("tax.csv"),   row.names = 1)

  dataset <- microeco::microtable$new(sample_table = grp, otu_table = otu, tax_table = tax)
  dataset <- microeco::tidy_taxonomy(dataset)

  dataset_1v2 <- dataset$clone(deep = TRUE)
  dataset_2v3 <- dataset$clone(deep = TRUE)
  dataset_1v3 <- dataset$clone(deep = TRUE)

  dataset_1v2$sample_table <- subset(dataset$sample_table, group %in% c("1","2"))
  dataset_1v2$otu_table    <- dataset$otu_table[, rownames(dataset_1v2$sample_table)]

  dataset_2v3$sample_table <- subset(dataset$sample_table, group %in% c("2","3"))
  dataset_2v3$otu_table    <- dataset$otu_table[, rownames(dataset_2v3$sample_table)]

  dataset_1v3$sample_table <- subset(dataset$sample_table, group %in% c("1","3"))
  dataset_1v3$otu_table    <- dataset$otu_table[, rownames(dataset_1v3$sample_table)]

  lefse_ready <- TRUE
  message("LEfSe datasets ready")
}, error = function(e) message("ERROR loading LEfSe data: ", e$message))


#Fig3D - LEfSe Stage 1 vs Stage 2
#D - LEfSe 1 vs 2

#Stage1Stage2LDA
#Stage1= Stage2
#Horizontal bar chart of LDA scores: Stage1 vs Stage2.
#Positive (right) = enriched in Stage1 (green); Negative (left) = enriched in Stage2 (orange).

message("\n=== Fig3D: LEfSe Stage1 vs Stage2 ===")
tryCatch({
  if (!lefse_ready) stop("LEfSe data not ready")
  lefse_1v2 <- microeco::trans_diff$new(
    dataset = dataset_1v2, method = "lefse", group = "group",
    p_adjust_method = "none", alpha = 0.01)

#stage1=, stage2=
#Custom colored bar chart (stage1=green, stage2=orange)
  p_D <- plot_lefse_custom(lefse_1v2, group1 = "1", group2 = "2",
                            colors = sc[c("1","2")],
                            top_n = 30, title = "Stage 1 vs Stage 2")

  ggsave(file.path(output_dir, "Figure3D_LEfSe_Stage1vs2.pdf"),
         plot = p_D, width = 8, height = 10, dpi = 300)
  write.csv(lefse_1v2$res_diff,
            file.path(output_dir, "Supplementary_Table_LEfSe_Stage1vs2_results.csv"),
            row.names = FALSE)
  message("Saved: Figure3D_LEfSe_Stage1vs2.pdf")
  message("Saved: Supplementary_Table_LEfSe_Stage1vs2_results.csv")
}, error = function(e) message("ERROR Fig3D: ", e$message))
p_D

#Fig3E - LEfSe Stage 2 vs Stage 3
#E - LEfSe 2 vs 3

#Stage2= Stage3
#Positive (right) = enriched in Stage2 (orange); Negative (left) = enriched in Stage3 (purple).

message("\n=== Fig3E: LEfSe Stage2 vs Stage3 ===")
tryCatch({
  if (!lefse_ready) stop("LEfSe data not ready")
  lefse_2v3 <- microeco::trans_diff$new(
    dataset = dataset_2v3, method = "lefse", group = "group",
    p_adjust_method = "none", alpha = 0.01)

  p_E <- plot_lefse_custom(lefse_2v3, group1 = "2", group2 = "3",
                            colors = sc[c("2","3")],
                            top_n = 30, title = "Stage 2 vs Stage 3")

  ggsave(file.path(output_dir, "Figure3E_LEfSe_Stage2vs3.pdf"),
         plot = p_E, width = 8, height = 10, dpi = 300)
  write.csv(lefse_2v3$res_diff,
            file.path(output_dir, "Supplementary_Table_LEfSe_Stage2vs3_results.csv"),
            row.names = FALSE)
  message("Saved: Figure3E_LEfSe_Stage2vs3.pdf")
  message("Saved: Supplementary_Table_LEfSe_Stage2vs3_results.csv")
}, error = function(e) message("ERROR Fig3E: ", e$message))
p_E

#Fig3F - LEfSe Stage 1 vs Stage 3
#F - LEfSe 1 vs 3

#Stage1= Stage3
#Positive (right) = enriched in Stage1 (green); Negative (left) = enriched in Stage3 (purple).

message("\n=== Fig3F: LEfSe Stage1 vs Stage3 ===")
tryCatch({
  if (!lefse_ready) stop("LEfSe data not ready")
  lefse_1v3 <- microeco::trans_diff$new(
    dataset = dataset_1v3, method = "lefse", group = "group",
    alpha = 0.01)

  p_F <- plot_lefse_custom(lefse_1v3, group1 = "1", group2 = "3",
                            colors = sc[c("1","3")],
                            top_n = 30, title = "Stage 1 vs Stage 3")

  ggsave(file.path(output_dir, "Figure3F_LEfSe_Stage1vs3.pdf"),
         plot = p_F, width = 8, height = 10, dpi = 300)
  write.csv(lefse_1v3$res_diff,
            file.path(output_dir, "Supplementary_Table_LEfSe_Stage1vs3_results.csv"),
            row.names = FALSE)
  message("Saved: Figure3F_LEfSe_Stage1vs3.pdf")
  message("Saved: Supplementary_Table_LEfSe_Stage1vs3_results.csv")
}, error = function(e) message("ERROR Fig3F: ", e$message))
p_F

#D-F / Combine D-F panels


message("\n=== Combining D-F panels ===")
tryCatch({
  if (exists("p_D") && exists("p_E") && exists("p_F")) {
    p_DEF <- cowplot::plot_grid(p_D, p_E, p_F, ncol = 3, labels = c("D","E","F"), align = "h")
    ggsave(file.path(output_dir, "Figure3D-F_LEfSe_combined.pdf"),
           plot = p_DEF, width = 24, height = 10, dpi = 300)
    message("Saved: Figure3D-F_LEfSe_combined.pdf")
  }
}, error = function(e) message("ERROR combining D-F: ", e$message))


#LEfSe / Supp: Three-stage LEfSe cladogram

#LEfSe
#Run LEfSe on all three stages and plot cladogram (supplementary figure).

message("\n=== SuppFig: LEfSe cladogram (all 3 stages) ===")
tryCatch({
  if (!lefse_ready) stop("LEfSe data not ready")
  lefse_all <- microeco::trans_diff$new(
    dataset = dataset, method = "lefse", group = "group",
    alpha = 0.05, lefse_subgroup = NULL)

  supp_cladogram <- lefse_all$plot_diff_cladogram(
    use_taxa_num = 150, use_feature_num = 50,
    clade_label_level = 5, group_order = c("1","2","3"))
  ggsave(file.path(output_dir, "SuppFig_LEfSe_cladogram_all3stages.pdf"),
         plot = supp_cladogram, width = 14, height = 14, dpi = 300)
  message("Saved: SuppFig_LEfSe_cladogram_all3stages.pdf")

  supp_bar_all <- lefse_all$plot_diff_bar(
    use_number = 1:30, width = 0.8, group_order = c("1","2","3"))
  ggsave(file.path(output_dir, "SuppFig_LEfSe_bar_all3stages.pdf"),
         plot = supp_bar_all, width = 8, height = 10, dpi = 300)
  message("Saved: SuppFig_LEfSe_bar_all3stages.pdf")
}, error = function(e) message("ERROR SuppFig LEfSe: ", e$message))
supp_cladogram

#Supp: Country-level bubble chart


#Fig 16).

message("\n=== SuppFig: Country-level bubble chart ===")
tryCatch({
  supp_iso_bubble <- build_overview_bubble(
    data_file    = "diff_taxa_counts_for_my5A.csv",
    fill_values  = expanded_colors,
    output_file  = file.path(output_dir, "Supplementary_Figure16_country_bubble.pdf"),
    label_levels = iso_list,
    bubble_size  = 2.2,
    height_in    = 20,
    width_in     = 20
  )
}, error = function(e) message("ERROR SuppFig country bubble: ", e$message))
supp_iso_bubble

#Supp: Country-level heatmaps

#P17
#Fig 17).

message("\n=== SuppFig: Country-level heatmaps ===")
tryCatch({
  iso_data     <- prepare_iso_data()
  iso_mean_bar <- build_mean_bar(iso_data$taxa_means)
  iso_pval     <- build_pval_heatmap(iso_data$pvals, iso_data$taxa_means,
                                      "region2", iso_list, expanded_colors)
  iso_mean     <- build_mean_heatmap(iso_data$mean_df, iso_data$taxa_means,
                                      "iso", iso_list, expanded_colors)
  supp_iso_hm  <- patchwork::wrap_elements(
    iso_pval$combined + iso_mean$combined + iso_mean_bar) +
    patchwork::inset_element(iso_pval$legend, 0.15, 1.02, 0.59, 1.10,
                              align_to = "plot", ignore_tag = TRUE) +
    patchwork::inset_element(iso_mean$legend, 0.60, 1.02, 1.00, 1.10,
                              align_to = "plot", ignore_tag = TRUE)
  ggsave(file.path(output_dir, "Supplementary_Figure17_country_heatmap.pdf"),
         plot = supp_iso_hm, width = 20, height = 10, dpi = 300)
  message("Saved: Supplementary_Figure17_country_heatmap.pdf")
}, error = function(e) message("ERROR SuppFig country heatmap: ", e$message))
supp_iso_hm

#Supp: Country-level ridge plots


#Fig 18).

message("\n=== SuppFig: Country-level ridge plots ===")
tryCatch({
  if (!exists("iso_data")) iso_data <- prepare_iso_data()
  supp_iso_ridges <- build_ridges(iso_data$genus_long, "iso", iso_list, expanded_colors)
  ggsave(file.path(output_dir, "Supplementary_Figure18_country_ridges.pdf"),
         plot = supp_iso_ridges, width = 12, height = 18, dpi = 300)
  message("Saved: Supplementary_Figure18_country_ridges.pdf")
}, error = function(e) message("ERROR SuppFig country ridges: ", e$message))
supp_iso_ridges

#Done


message("\n========== Figure 3 outputs (", output_dir, ") ==========")
cat_if_exists <- function(f) {
  path <- file.path(output_dir, f)
  cat(sprintf("  [%s] %s\n", if (file.exists(path)) "OK" else "MISSING", f))
}
main_files <- c(
  "Figure3A_stage_bubble.pdf",
  "Figure3B_stage_heatmap.pdf",
  "Figure3C_stage_ridges.pdf",
  "Figure3D_LEfSe_Stage1vs2.pdf",
  "Figure3E_LEfSe_Stage2vs3.pdf",
  "Figure3F_LEfSe_Stage1vs3.pdf",
  "Figure3D-F_LEfSe_combined.pdf"
)
supp_files <- c(
  "SuppFig_LEfSe_cladogram_all3stages.pdf",
  "SuppFig_LEfSe_bar_all3stages.pdf",
  "Supplementary_Figure16_country_bubble.pdf",
  "Supplementary_Figure17_country_heatmap.pdf",
  "Supplementary_Figure18_country_ridges.pdf",
  "Supplementary_Table_LEfSe_Stage1vs2_results.csv",
  "Supplementary_Table_LEfSe_Stage2vs3_results.csv",
  "Supplementary_Table_LEfSe_Stage1vs3_results.csv"
)
message("--- Main figures ---")
invisible(lapply(main_files,  cat_if_exists))
message("--- Supplementary ---")
invisible(lapply(supp_files, cat_if_exists))


#Figure 3 / Complete Figure 3 combined layout (A-F)

#A-F
#Combines all panels A-F into one publication-ready figure.
#Layout:
#Row 1: [A (bubble)] [B (heatmap, wider)]
#Row 2: [C (ridge plots, full width)]
#Row 3: [D (LEfSe 1v2)] [E (LEfSe 2v3)] [F (LEfSe 1v3)]

message("\n=== Building complete Figure 3 combined layout ===")
tryCatch({
  if (!exists("p_A") || !exists("p_B") || !exists("p_C") ||
      !exists("p_D") || !exists("p_E") || !exists("p_F")) {
    stop("One or more panels missing; run full script first.")
  }

#A + BAB / Row 1: A + B (A narrow, B wide)
  row1 <- cowplot::plot_grid(
    p_A, p_B,
    ncol       = 2,
    rel_widths = c(1, 2.2),
    labels     = c("A", "B"),
    label_size = 14
  )

#C / Row 2: C full width
  row2 <- cowplot::plot_grid(
    p_C,
    ncol       = 1,
    labels     = "C",
    label_size = 14
  )

#D + E + F / Row 3: D + E + F equal width
  row3 <- cowplot::plot_grid(
    p_D, p_E, p_F,
    ncol       = 3,
    labels     = c("D", "E", "F"),
    label_size = 14,
    align      = "h"
  )

#Stack rows vertically
  fig3_complete <- cowplot::plot_grid(
    row1, row2, row3,
    ncol        = 1,
    rel_heights = c(1.4, 0.9, 1.5)
  )

  ggsave(
    filename = file.path(output_dir, "Figure3_complete_A-F.pdf"),
    plot     = fig3_complete,
    width    = 22,
    height   = 28,
    dpi      = 300
  )
  message("Saved: Figure3_complete_A-F.pdf")
  cat_if_exists("Figure3_complete_A-F.pdf")
}, error = function(e) message("ERROR building complete Fig3: ", e$message))
