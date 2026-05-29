
#Figure 4 / 4 - Master Script (Panels A–J)
#Systematic differences in microbiome composition between IBD and healthy individuals
#IBD

#Panel summary / :
#A – Genus-level bubble chart: pairwise number of DA genera (NC / UC / CD)
#NC / UC / CD
#B – Species-level bubble chart: pairwise number of DA species (NC / UC / CD)
#NC / UC / CD
#C – Genus-level ridgeline plots (NC / UC / CD, top 10 genera)
#NC / UC / CD10
#D – Species-level ridgeline plots (NC / UC / CD, top 10 species)
#NC / UC / CD10
#E – Genus-level ridgeline plots (NC / IBD, top 10 genera)
#NC / IBD10
#F – Species-level ridgeline plots (NC / IBD, top 10 species)
#NC / IBD10
#G – Genus-level heatmap composite: adj p-value + mean abundance (NC / UC / CD)
#p + NC / UC / CD
#H – Genus-level heatmap composite: adj p-value + mean abundance (NC / IBD)
#p + NC / IBD
#I – Species-level heatmap composite: adj p-value + mean abundance (NC / UC / CD)
#p + NC / UC / CD
#J – Species-level heatmap composite: adj p-value + mean abundance (NC / IBD)
#p + NC / IBD

#Input files required / :
#Supplementary_Table_20__relative_abundance.tsv
#Supplementary_Table_21__metadata_common.tsv
#taxaNamesIBD.csv
#diff_taxa_counts_for_my5ABCD.csv (genus-level pairwise DA counts)
#diff_taxa_counts_for_my5ABCD_species.csv (species-level pairwise DA counts)
#regionDF-diseasesubtype.csv
#regionDF-disease.csv
#Supplementary_Table_22_DA_LMM_emmeans_results_diseasesubtype.csv
#Supplementary_Table_23_DA_LMM_emmeans_results_disease.csv
#Supplementary_Table_24_DA_LMM_emmeans_results_disease_genus.csv
#Supplementary_Table_25_DA_LMM_emmeans_results_diseasesubtype_genus.csv




#Load packages /

pkgs <- c("readr", "dplyr", "ggplot2", "cowplot", "patchwork", "grid",
          "tidyr", "ggridges", "viridis", "stringr", "tibble")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE))



#Global settings and helper functions /

output_dir <- "figure4"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Disease color palettes /
#group (NC / UC / CD) /
diseasecolors_3 <- c(
  "NC" = "#1B9E77",#green /
  "UC" = "#D95F02",#orange /
  "CD" = "#7570B3"#purple /
)
#group (NC / IBD metadata may store healthy / IBD or NC / IBD) /
diseasecolors_2 <- c(
  "healthy" = "#1B9E77",#NC = healthy / NChealthy
  "IBD"     = "#D95F02"
)

#Relative abundance helper /
rel <- function(x) x / sum(x, na.rm = TRUE)
make_rel <- function(table) {
  data.frame(t(apply(table, 1, rel)))
}

#File finder /
find_file <- function(filename) {
  candidates <- c(filename, file.path("D:/R代码", filename))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("File not found: ", filename)
  hit
}



#Shared data loading /


#Metadata /
metadata <- read.delim(
  find_file("Supplementary_Table_21__metadata_common.tsv"), sep = "\t"
)
metadata$sample <- as.character(metadata$sample)

#Taxon name mapping /
taxaNames_all <- read.csv(find_file("taxaNamesIBD.csv"), stringsAsFactors = FALSE)
taxaNames_all$full.taxon <- gsub("\\|", ".", taxaNames_all$full.taxon)



#Helper: build relative-abundance matrix and join metadata


build_rel_matrix <- function(abund_wide, meta, group_col) {
#abund_wide: samples x taxa data.frame (numeric, possibly character)
  abund_wide[] <- lapply(abund_wide, function(x) as.numeric(trimws(as.character(x))))
  mat <- as.matrix(abund_wide) + 1e-6
  rel_df <- make_rel(mat)
  rownames(rel_df) <- rownames(abund_wide)
#Align with metadata /
  common <- intersect(meta$sample, rownames(rel_df))
  rel_df <- rel_df[common, , drop = FALSE]
  meta_sub <- meta[meta$sample %in% common, ]
  meta_sub <- meta_sub[match(rownames(rel_df), meta_sub$sample), ]
#Attach group column /
  rel_df[[group_col]] <- meta_sub[[group_col]]
  rel_df
}

#Helper: rename columns using taxaNames lookup
#taxaNames
rename_taxa_cols <- function(rel_df, taxaNames, pvals_taxon_vec) {
  tn <- taxaNames[taxaNames$taxon %in% pvals_taxon_vec, ]
#Keep only columns whose taxon names appear in pvals / pvals
  keep <- colnames(rel_df) %in% tn$taxon
  rel_df[, keep, drop = FALSE]
}




#Fig A / A
#Genus-level bubble chart: pairwise DA genera counts (NC / UC / CD)
#NC / UC / CD


#p_a

#Compute genus DA counts from Supplementary_Table_25 (AdjP_BH < 0.05)
#Supplementary_Table_25 DAAdjP_BH < 0.05
pvals_g3_raw_forA <- read_csv(
  find_file("Supplementary_Table_25_DA_LMM_emmeans_results_diseasesubtype_genus.csv"),
  show_col_types = FALSE
) %>%
  separate(Comparison, into = c("region1", "region2"), sep = " - ")
taxaCounts_genus <- pvals_g3_raw_forA %>%
  filter(AdjP_BH < 0.05) %>%
  group_by(region1, region2) %>%
  summarise(freq = n_distinct(Taxon), .groups = "drop") %>%
  mutate(regions = paste0(region1, "_", region2), regionLabel = region2)
NUM_TAXA_G <- n_distinct(pvals_g3_raw_forA$Taxon)#total evaluated genera /
plot_ord_g <- rev(unique(c(taxaCounts_genus$region1, taxaCounts_genus$region2)))
labels_g   <- data.frame(region2 = rev(plot_ord_g), regionLabel = rev(plot_ord_g))

bubble_g <- ggplot() +
  geom_point(
    data = taxaCounts_genus,
    aes(x = factor(region1, levels = rev(plot_ord_g)),
        y = factor(region2, levels = rev(plot_ord_g)),
        size = freq, fill = "grey"),
    shape = 21, stroke = 0
  ) +
  geom_text(
    data = taxaCounts_genus,
    aes(x = factor(region1, levels = rev(plot_ord_g)),
        y = factor(region2, levels = rev(plot_ord_g)),
        label = freq),
    size = 5, vjust = 0
  ) +
  geom_label(
    data = labels_g,
    aes(x = factor(region2, levels = rev(plot_ord_g)),
        y = factor(region2, levels = rev(plot_ord_g)),
        label = regionLabel, fill = regionLabel),
    size = 5
  ) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(
    legend.position     = "top",
    axis.text           = element_blank(),
    axis.ticks          = element_blank(),
    legend.title        = element_text(size = 8),
    panel.spacing       = unit(c(0, 0, 0, 0), "pt"),
    plot.margin         = unit(c(0, 0, 0, 0), "pt")
  ) +
  scale_radius(range = c(2, 22), limits = c(0, 70),
               breaks = c(10, 20, 30, 50, NUM_TAXA_G)) +
  guides(
    size  = guide_legend(override.aes = list(fill = "grey", stroke = .25),
                         label.position = "bottom", title.position = "top", order = 1),
    fill  = FALSE, color = FALSE
  ) +
  labs(size = "Number of \nDifferentially \nAbundant Taxa", x = NULL, y = NULL)

size_legend_g <- cowplot::get_legend(bubble_g)
bubble_g <- bubble_g + theme(legend.position = "none") +
  inset_element(size_legend_g, 0.15, 0.60, 0.40, 0.95, align_to = "panel")

label_tile_g <- ggplot(labels_g) +
  geom_tile(aes(x = factor(region2, levels = plot_ord_g),
                y = 1,
                fill = factor(regionLabel, levels = rev(plot_ord_g)))) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.spacing = unit(c(0, 0, 0, 0), "pt"),
        plot.margin   = unit(c(0, 0, 0, 0), "pt")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))

p_a <- wrap_elements(bubble_g + label_tile_g + plot_layout(heights = c(39, 1)))
p_a




#Fig B / B
#Species-level bubble chart: pairwise DA species counts (NC / UC / CD)
#NC / UC / CD


#p_b

#Load pre-computed species DA counts (NC / UC / CD pairwise)
#DANC / UC / CD
taxaCounts_species <- read_csv(find_file("diff_taxa_counts_for_my5ABCD.csv"),
                                show_col_types = FALSE)
NUM_TAXA_S  <- 42#total evaluated species /
plot_ord_s  <- rev(unique(c(taxaCounts_species$region1, taxaCounts_species$region2)))
labels_s    <- data.frame(region2 = rev(plot_ord_s), regionLabel = rev(plot_ord_s))

bubble_s <- ggplot() +
  geom_point(
    data = taxaCounts_species,
    aes(x = factor(region1, levels = rev(plot_ord_s)),
        y = factor(region2, levels = rev(plot_ord_s)),
        size = freq, fill = "grey"),
    shape = 21, stroke = 0
  ) +
  geom_text(
    data = taxaCounts_species,
    aes(x = factor(region1, levels = rev(plot_ord_s)),
        y = factor(region2, levels = rev(plot_ord_s)),
        label = freq),
    size = 5, vjust = 0
  ) +
  geom_label(
    data = labels_s,
    aes(x = factor(region2, levels = rev(plot_ord_s)),
        y = factor(region2, levels = rev(plot_ord_s)),
        label = regionLabel, fill = regionLabel),
    size = 5
  ) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(
    legend.position = "top",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.title    = element_text(size = 8),
    panel.spacing   = unit(c(0, 0, 0, 0), "pt"),
    plot.margin     = unit(c(0, 0, 0, 0), "pt")
  ) +
  scale_radius(range = c(2, 22), limits = c(0, 70),
               breaks = c(10, 20, 30, 50, NUM_TAXA_S)) +
  guides(
    size  = guide_legend(override.aes = list(fill = "grey", stroke = .25),
                         label.position = "bottom", title.position = "top", order = 1),
    fill  = FALSE, color = FALSE
  ) +
  labs(size = "Number of \nDifferentially \nAbundant Taxa", x = NULL, y = NULL)

size_legend_s <- cowplot::get_legend(bubble_s)
bubble_s <- bubble_s + theme(legend.position = "none") +
  inset_element(size_legend_s, 0.15, 0.60, 0.40, 0.95, align_to = "panel")

label_tile_s <- ggplot(labels_s) +
  geom_tile(aes(x = factor(region2, levels = plot_ord_s),
                y = 1,
                fill = factor(regionLabel, levels = rev(plot_ord_s)))) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.spacing = unit(c(0, 0, 0, 0), "pt"),
        plot.margin   = unit(c(0, 0, 0, 0), "pt")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))

p_b <- wrap_elements(bubble_s + label_tile_s + plot_layout(heights = c(39, 1)))
p_b




#Genus-level data preprocessing (shared for C, G, E, H)
#CGEH



#Aggregate species table to genus level /
abund_raw_genus <- read.delim(
  find_file("Supplementary_Table_20__relative_abundance.tsv"), sep = "\t"
)
abund_raw_genus$genus <- sub(".*(g__[^|]+).*", "\\1", abund_raw_genus$taxa)
abund_raw_genus$taxa  <- abund_raw_genus$genus
abund_raw_genus$genus <- NULL

abund_genus_agg <- abund_raw_genus %>%
  group_by(taxa) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

#Transpose: rows = samples, cols = genera / ==
genus_t <- as.data.frame(t(abund_genus_agg))
colnames(genus_t) <- genus_t[1, ]
genus_t <- genus_t[-1, , drop = FALSE]
rownames(genus_t) <- as.character(rownames(genus_t))

#Make relative abundance matrix /
genus_t[] <- lapply(genus_t, function(x) as.numeric(trimws(as.character(x))))
genus_mat  <- as.matrix(genus_t) + 1e-6
genus_rel_base <- make_rel(genus_mat)#samples x genera, relative abundances

#Keep only samples present in metadata /
common_genus <- intersect(metadata$sample, rownames(genus_rel_base))
genus_rel_base <- genus_rel_base[common_genus, , drop = FALSE]
meta_genus <- metadata[match(common_genus, metadata$sample), ]




#Fig C / C
#Genus-level ridgeline plots (NC / UC / CD), top 10 genera
#NC / UC / CD10


#p_c

#Load LMM p-values for genus × disease_subtype / ×LMM p
pvals_g3 <- read_csv(
  find_file("Supplementary_Table_25_DA_LMM_emmeans_results_diseasesubtype_genus.csv")
) %>%
  separate(Comparison, into = c("region1", "region2"), sep = " - ")

#Genus pvals already use short names (g__xxx) no taxaNames lookup needed
#pvals g__xxx taxaNames
all_taxa_g3 <- unique(pvals_g3$Taxon)#all genera from LMM / LMM

#Build genus rel matrix only taxa present in abundance table

valid_g3   <- intersect(colnames(genus_rel_base), all_taxa_g3)
genus_rel3 <- as.data.frame(genus_rel_base[, valid_g3, drop = FALSE])
genus_rel3$disease_subtype <- meta_genus$disease_subtype[
  match(rownames(genus_rel3), meta_genus$sample)]
genus_rel3$disease_subtype <- factor(genus_rel3$disease_subtype,
                                     levels = c("NC", "UC", "CD"))

#ALL pvals taxa sorted by mean abundance used for heatmap y-axis
#pvalsy0
taxa_g3_cols <- setdiff(colnames(genus_rel3), "disease_subtype")
taxa_means_g3 <- data.frame(
  taxon      = all_taxa_g3,
  taxa.means = sapply(all_taxa_g3, function(t) {
    if (t %in% taxa_g3_cols) mean(genus_rel3[[t]], na.rm = TRUE) else 0
  })
) %>% arrange(desc(taxa.means))
#Top 10 for ridges only taxa with abundance data

top10_g3 <- taxa_means_g3$taxon[taxa_means_g3$taxon %in% taxa_g3_cols][1:10]

#Long format for ridgeline /
genus_long_g3 <- genus_rel3[, c(top10_g3, "disease_subtype")] %>%
  pivot_longer(-disease_subtype, names_to = "taxon", values_to = "rel_abundance")
genus_long_g3$taxon <- factor(genus_long_g3$taxon, levels = top10_g3)

p_c <- wrap_elements(
  ggplot(data = genus_long_g3) +
    geom_density_ridges(
      aes(x = log10(rel_abundance),
          y = factor(disease_subtype, levels = c("NC", "UC", "CD")),
          fill = disease_subtype),
      quantile_lines = TRUE, quantiles = 2
    ) +
    scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
    theme(
      legend.position = "none",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      text            = element_text(color = "black"),
      plot.margin     = unit(c(10, 10, 10, 10), "pt")
    ) +
    labs(x = "Relative Abundance", y = "Density") +
    scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
    facet_wrap(~taxon, ncol = 5,
               labeller = labeller(taxon = function(x) sub("^g__", "", x)))
)
p_c




#Fig G / G
#Genus-level heatmap composite: adj p-value + mean abundance (NC / UC / CD)
#p + NC / UC / CD


#p_g

#Bar chart of mean relative abundances /
taxa_means_g3$taxon <- factor(taxa_means_g3$taxon, levels = taxa_means_g3$taxon)

mean_bar_g3 <- ggplot(data = taxa_means_g3) +
  geom_col(aes(x = taxa.means, y = taxon), fill = "gray30") +
  xlab("Mean Relative \nAbundance") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "pt"),
    plot.margin   = unit(c(0, 0, 0, 0), "pt"),
    axis.title.x  = element_text(size = 7),
    axis.text.x   = element_text(size = 6)
  ) +
  scale_y_discrete(limits = rev(levels(taxa_means_g3$taxon)))

#Adj p-value heatmap (vs NC) / pNC
pvals_g3_nc <- pvals_g3 %>%
  filter(region1 == "NC" | region2 == "NC") %>%
  mutate(
    region2 = ifelse(region1 == "NC", region2, region1),
    region1 = "NC"
  )

pval_hm_g3 <- ggplot(data = pvals_g3_nc %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = factor(region2, levels = c("UC", "CD")),
    y    = factor(Taxon, levels = rev(as.character(taxa_means_g3$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 6),
    axis.title    = element_blank(),
    axis.ticks    = element_blank(),
    legend.position = "top",
    legend.title  = element_text(size = 9),
    text          = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title          = "Adjusted P Value (vs NC)",
    title.position = "top", title.hjust = 0.5
  ))

pval_legend_g3 <- cowplot::get_legend(pval_hm_g3)
pval_hm_g3 <- pval_hm_g3 + theme(legend.position = "none")

#Color strip below p-value heatmap / p
pval_strip_g3 <- ggplot(
  data.frame(region = factor(c("UC", "CD"), levels = c("UC", "CD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_3[c("UC", "CD")]) +
  scale_x_discrete(limits = c("UC", "CD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_pval_g3 <- pval_hm_g3 +
  inset_element(pval_strip_g3, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

#Mean relative abundance heatmap (NC / UC / CD) /
regionDF_3 <- read_csv(find_file("regionDF-diseasesubtype.csv"))
if (exists("counts_g3")) rm(counts_g3)
for (i in seq_len(nrow(regionDF_3))) {
  grp   <- regionDF_3$region[i]
  idx   <- !is.na(genus_rel3$disease_subtype) & genus_rel3$disease_subtype == grp
  cnts  <- genus_rel3[idx, taxa_g3_cols, drop = FALSE]
  mnDF  <- data.frame(
    taxon  = colnames(cnts),
    mean   = colMeans(cnts, na.rm = TRUE),
    region = grp
  )
  counts_g3 <- if (!exists("counts_g3")) mnDF else rbind(counts_g3, mnDF)
}
counts_g3$region <- factor(counts_g3$region, levels = c("NC", "UC", "CD"))

mean_hm_g3 <- ggplot(counts_g3) +
  geom_tile(aes(
    x    = factor(region, levels = c("NC", "UC", "CD")),
    y    = factor(taxon,  levels = rev(as.character(taxa_means_g3$taxon))),
    fill = log10(mean)
  )) +
  scale_fill_viridis() +
  scale_x_discrete(limits = c("NC", "UC", "CD"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing  = unit(c(0, 0, 0, 0), "pt"),
    plot.margin    = unit(c(0, 0, 0, 0), "pt"),
    legend.position = "top",
    legend.title   = element_text(size = 9),
    legend.text    = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(
    title = "Mean Relative Abundance",
    title.position = "top", title.hjust = 0.5
  ))

mean_legend_g3 <- cowplot::get_legend(mean_hm_g3)
mean_hm_g3 <- mean_hm_g3 + theme(legend.position = "none")

mean_strip_g3 <- ggplot(
  data.frame(region = factor(c("NC", "UC", "CD"), levels = c("NC", "UC", "CD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_3) +
  scale_x_discrete(limits = c("NC", "UC", "CD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_mean_g3 <- mean_hm_g3 +
  inset_element(mean_strip_g3, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

#Combine into p_g / p_g
bars_g3 <- wrap_elements(p_pval_g3 + p_mean_g3 + mean_bar_g3)
blank_g  <- ggplot() + theme_void()

p_g <- blank_g + bars_g3 +
  inset_element(pval_legend_g3,  0.05, 1.01, 0.55, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  inset_element(mean_legend_g3,  0.56, 1.01, 1.00, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  plot_layout(ncol = 1)
p_g




#Fig E / E
#Genus-level ridgeline plots (NC / IBD), top 10 genera
#NC / IBD10


#p_e

pvals_g2 <- read_csv(
  find_file("Supplementary_Table_24_DA_LMM_emmeans_results_disease_genus.csv")
) %>%
  separate(Comparison, into = c("region1", "region2"), sep = " - ")

#Genus pvals already use short names (g__xxx) no taxaNames lookup needed
#pvals g__xxx taxaNames
all_taxa_g2 <- unique(pvals_g2$Taxon)

valid_g2   <- intersect(colnames(genus_rel_base), all_taxa_g2)
genus_rel2 <- as.data.frame(genus_rel_base[, valid_g2, drop = FALSE])
genus_rel2$disease <- meta_genus$disease[
  match(rownames(genus_rel2), meta_genus$sample)]
genus_rel2$disease <- factor(genus_rel2$disease, levels = c("healthy", "IBD"))

taxa_g2_cols <- setdiff(colnames(genus_rel2), "disease")
taxa_means_g2 <- data.frame(
  taxon      = all_taxa_g2,
  taxa.means = sapply(all_taxa_g2, function(t) {
    if (t %in% taxa_g2_cols) mean(genus_rel2[[t]], na.rm = TRUE) else 0
  })
) %>% arrange(desc(taxa.means))
top10_g2 <- taxa_means_g2$taxon[taxa_means_g2$taxon %in% taxa_g2_cols][1:10]

genus_long_g2 <- genus_rel2[, c(top10_g2, "disease")] %>%
  pivot_longer(-disease, names_to = "taxon", values_to = "rel_abundance")
genus_long_g2$taxon   <- factor(genus_long_g2$taxon,   levels = top10_g2)
genus_long_g2$disease <- factor(genus_long_g2$disease, levels = c("healthy", "IBD"))

p_e <- wrap_elements(
  ggplot(data = genus_long_g2) +
    geom_density_ridges(
      aes(x    = log10(rel_abundance),
          y    = factor(disease, levels = c("healthy", "IBD")),
          fill = disease),
      quantile_lines = TRUE, quantiles = 2
    ) +
    scale_fill_manual(values = diseasecolors_2, aesthetics = c("colour", "fill")) +
    theme(
      legend.position = "none",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      text            = element_text(color = "black"),
      plot.margin     = unit(c(10, 10, 10, 10), "pt")
    ) +
    labs(x = "Relative Abundance", y = "Density") +
    scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
    facet_wrap(~taxon, ncol = 5,
               labeller = labeller(taxon = function(x) sub("^g__", "", x)))
)
p_e




#Fig H / H
#Genus-level heatmap composite: adj p-value + mean abundance (NC / IBD)
#p + NC / IBD


#p_h

taxa_means_g2$taxon <- factor(taxa_means_g2$taxon, levels = taxa_means_g2$taxon)

mean_bar_g2 <- ggplot(data = taxa_means_g2) +
  geom_col(aes(x = taxa.means, y = taxon), fill = "gray30") +
  xlab("Mean Relative \nAbundance") +
  theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "pt"),
    plot.margin   = unit(c(0, 0, 0, 0), "pt"),
    axis.title.x  = element_text(size = 7),
    axis.text.x   = element_text(size = 6)
  ) +
  scale_y_discrete(limits = rev(levels(taxa_means_g2$taxon)))

pval_hm_g2 <- ggplot(data = pvals_g2 %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = factor(region2, levels = c("healthy", "IBD")),
    y    = factor(Taxon, levels = rev(as.character(taxa_means_g2$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 6),
    axis.title    = element_blank(), axis.ticks = element_blank(),
    legend.position = "top", legend.title = element_text(size = 9),
    text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title = "Adjusted P Value", title.position = "top", title.hjust = 0.5
  ))

pval_legend_g2 <- cowplot::get_legend(pval_hm_g2)
pval_hm_g2 <- pval_hm_g2 + theme(legend.position = "none")

pval_strip_g2 <- ggplot(
  data.frame(region = factor(c("healthy", "IBD"), levels = c("healthy", "IBD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_2) +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_pval_g2 <- pval_hm_g2 +
  inset_element(pval_strip_g2, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

regionDF_2 <- read_csv(find_file("regionDF-disease.csv"))
if (exists("counts_g2")) rm(counts_g2)
for (i in seq_len(nrow(regionDF_2))) {
  grp  <- regionDF_2$region[i]
  idx  <- !is.na(genus_rel2$disease) & genus_rel2$disease == grp
  cnts <- genus_rel2[idx, taxa_g2_cols, drop = FALSE]
  mnDF <- data.frame(
    taxon  = colnames(cnts),
    mean   = colMeans(cnts, na.rm = TRUE),
    region = grp
  )
  counts_g2 <- if (!exists("counts_g2")) mnDF else rbind(counts_g2, mnDF)
}
counts_g2$region <- factor(counts_g2$region, levels = c("healthy", "IBD"))

mean_hm_g2 <- ggplot(counts_g2) +
  geom_tile(aes(
    x    = factor(region, levels = c("healthy", "IBD")),
    y    = factor(taxon,  levels = rev(as.character(taxa_means_g2$taxon))),
    fill = log10(mean)
  )) +
  scale_fill_viridis() +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing   = unit(c(0, 0, 0, 0), "pt"),
    plot.margin     = unit(c(0, 0, 0, 0), "pt"),
    legend.position = "top",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(
    title = "Mean Relative Abundance", title.position = "top", title.hjust = 0.5
  ))

mean_legend_g2 <- cowplot::get_legend(mean_hm_g2)
mean_hm_g2 <- mean_hm_g2 + theme(legend.position = "none")

mean_strip_g2 <- ggplot(
  data.frame(region = factor(c("healthy", "IBD"), levels = c("healthy", "IBD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_2) +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_mean_g2 <- mean_hm_g2 +
  inset_element(mean_strip_g2, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

bars_g2 <- wrap_elements(p_pval_g2 + p_mean_g2 + mean_bar_g2)
blank_h  <- ggplot() + theme_void()

p_h <- blank_h + bars_g2 +
  inset_element(pval_legend_g2,  0.05, 1.01, 0.55, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  inset_element(mean_legend_g2,  0.56, 1.01, 1.00, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  plot_layout(ncol = 1)
p_h




#Species-level data preprocessing (shared for D, I, F, J)
#DIFJ



abund_raw_species <- read.delim(
  find_file("Supplementary_Table_20__relative_abundance.tsv"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)
species_t <- as.data.frame(t(abund_raw_species))
rownames(species_t) <- as.character(rownames(species_t))

#Make relative abundance matrix /
species_t[] <- lapply(species_t, function(x) as.numeric(trimws(as.character(x))))
species_mat       <- as.matrix(species_t) + 1e-6
species_rel_base  <- make_rel(species_mat)

#Align with metadata /
common_species    <- intersect(metadata$sample, rownames(species_rel_base))
species_rel_base  <- species_rel_base[common_species, , drop = FALSE]
meta_species      <- metadata[match(common_species, metadata$sample), ]

#Convert column names from | to . to match taxaNames_all$full.taxon format
#taxaNames_all$full.taxon
colnames(species_rel_base) <- gsub("\\|", ".", colnames(species_rel_base))




#Fig D / D
#Species-level ridgeline plots (NC / UC / CD), top 10 species
#NC / UC / CD10


#p_d

pvals_sp3 <- read_csv(
  find_file("Supplementary_Table_22_DA_LMM_emmeans_results_diseasesubtype.csv")
) %>%
  separate(Comparison, into = c("region1", "region2"), sep = " - ")

taxaNames_sp3 <- taxaNames_all[taxaNames_all$full.taxon %in% pvals_sp3$Taxon, ]
#Extract species short name (no prefix) from full taxonomy

taxaNames_sp3$species_name <- sub(".*\\.s__", "", taxaNames_sp3$full.taxon)
for (i in seq_len(nrow(taxaNames_sp3))) {
  pvals_sp3$Taxon[pvals_sp3$Taxon == taxaNames_sp3$full.taxon[i]] <- taxaNames_sp3$species_name[i]
}

valid_sp3    <- colnames(species_rel_base) %in% taxaNames_sp3$full.taxon
species_rel3 <- species_rel_base[, valid_sp3, drop = FALSE]
ord_sp3      <- match(colnames(species_rel3), taxaNames_sp3$full.taxon)
colnames(species_rel3) <- taxaNames_sp3$species_name[ord_sp3]

species_rel3 <- as.data.frame(species_rel3)
species_rel3$disease_subtype <- meta_species$disease_subtype[
  match(rownames(species_rel3), meta_species$sample)]
species_rel3$disease_subtype <- factor(species_rel3$disease_subtype,
                                       levels = c("NC", "UC", "CD"))

taxa_sp3_cols <- setdiff(colnames(species_rel3), "disease_subtype")
#All species from pvals (for heatmap y-axis) / pvalsy
all_taxa_sp3  <- unique(pvals_sp3$Taxon)
taxa_means_sp3 <- data.frame(
  taxon      = all_taxa_sp3,
  taxa.means = sapply(all_taxa_sp3, function(t) {
    if (t %in% taxa_sp3_cols) mean(species_rel3[[t]], na.rm = TRUE) else 0
  })
) %>% arrange(desc(taxa.means))
top10_sp3 <- taxa_means_sp3$taxon[taxa_means_sp3$taxon %in% taxa_sp3_cols][1:10]

species_long_sp3 <- species_rel3[, c(top10_sp3, "disease_subtype")] %>%
  pivot_longer(-disease_subtype, names_to = "taxon", values_to = "rel_abundance")
species_long_sp3$taxon <- factor(species_long_sp3$taxon, levels = top10_sp3)

p_d <- wrap_elements(
  ggplot(data = species_long_sp3) +
    geom_density_ridges(
      aes(x    = log10(rel_abundance),
          y    = factor(disease_subtype, levels = c("NC", "UC", "CD")),
          fill = disease_subtype),
      quantile_lines = TRUE, quantiles = 2
    ) +
    scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
    theme(
      legend.position = "none",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      text            = element_text(color = "black"),
      plot.margin     = unit(c(10, 10, 10, 10), "pt")
    ) +
    labs(x = "Relative Abundance", y = "Density") +
    scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
    facet_wrap(~taxon, ncol = 5)
)
p_d




#Fig I / I
#Species-level heatmap composite: adj p-value + mean abundance (NC / UC / CD)
#p + NC / UC / CD


#p_i

taxa_means_sp3$taxon <- factor(taxa_means_sp3$taxon, levels = taxa_means_sp3$taxon)

mean_bar_sp3 <- ggplot(data = taxa_means_sp3) +
  geom_col(aes(x = taxa.means, y = taxon), fill = "gray30") +
  xlab("Mean Relative \nAbundance") +
  theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "pt"),
    plot.margin   = unit(c(0, 0, 0, 0), "pt"),
    axis.title.x  = element_text(size = 7),
    axis.text.x   = element_text(size = 6)
  ) +
  scale_x_continuous(breaks = c(0, 0.1)) +
  scale_y_discrete(limits = rev(levels(taxa_means_sp3$taxon)))

pvals_sp3_nc <- pvals_sp3 %>%
  filter(region1 == "NC" | region2 == "NC") %>%
  mutate(
    region2 = ifelse(region1 == "NC", region2, region1),
    region1 = "NC"
  )

pval_hm_sp3 <- ggplot(data = pvals_sp3_nc %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = factor(region2, levels = c("UC", "CD")),
    y    = factor(Taxon, levels = rev(as.character(taxa_means_sp3$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 6),
    axis.title    = element_blank(), axis.ticks = element_blank(),
    legend.position = "top", legend.title = element_text(size = 9),
    text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels = function(x) paste0("s__", x)) +
  guides(fill = guide_colorbar(
    title = "Adjusted P Value (vs NC)", title.position = "top", title.hjust = 0.5
  ))

pval_legend_sp3 <- cowplot::get_legend(pval_hm_sp3)
pval_hm_sp3 <- pval_hm_sp3 + theme(legend.position = "none")

pval_strip_sp3 <- ggplot(
  data.frame(region = factor(c("UC", "CD"), levels = c("UC", "CD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_3[c("UC", "CD")]) +
  scale_x_discrete(limits = c("UC", "CD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_pval_sp3 <- pval_hm_sp3 +
  inset_element(pval_strip_sp3, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

if (exists("counts_sp3")) rm(counts_sp3)
for (i in seq_len(nrow(regionDF_3))) {
  grp  <- regionDF_3$region[i]
  idx  <- !is.na(species_rel3$disease_subtype) & species_rel3$disease_subtype == grp
  cnts <- species_rel3[idx, taxa_sp3_cols, drop = FALSE]
  mnDF <- data.frame(
    taxon  = colnames(cnts),
    mean   = colMeans(cnts, na.rm = TRUE),
    region = grp
  )
  counts_sp3 <- if (!exists("counts_sp3")) mnDF else rbind(counts_sp3, mnDF)
}
counts_sp3$region <- factor(counts_sp3$region, levels = c("NC", "UC", "CD"))

mean_hm_sp3 <- ggplot(counts_sp3) +
  geom_tile(aes(
    x    = factor(region, levels = c("NC", "UC", "CD")),
    y    = factor(taxon,  levels = rev(as.character(taxa_means_sp3$taxon))),
    fill = log10(mean)
  )) +
  scale_fill_viridis() +
  scale_x_discrete(limits = c("NC", "UC", "CD"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing   = unit(c(0, 0, 0, 0), "pt"),
    plot.margin     = unit(c(0, 0, 0, 0), "pt"),
    legend.position = "top",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(
    title = "Mean Relative Abundance", title.position = "top", title.hjust = 0.5
  ))

mean_legend_sp3 <- cowplot::get_legend(mean_hm_sp3)
mean_hm_sp3 <- mean_hm_sp3 + theme(legend.position = "none")

mean_strip_sp3 <- ggplot(
  data.frame(region = factor(c("NC", "UC", "CD"), levels = c("NC", "UC", "CD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_3) +
  scale_x_discrete(limits = c("NC", "UC", "CD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_mean_sp3 <- mean_hm_sp3 +
  inset_element(mean_strip_sp3, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

bars_sp3 <- wrap_elements(p_pval_sp3 + p_mean_sp3 + mean_bar_sp3)
blank_i  <- ggplot() + theme_void()

p_i <- blank_i + bars_sp3 +
  inset_element(pval_legend_sp3, 0.05, 1.01, 0.55, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  inset_element(mean_legend_sp3, 0.56, 1.01, 1.00, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  plot_layout(ncol = 1)
p_i




#Fig F / F
#Species-level ridgeline plots (NC / IBD), top 10 species
#NC / IBD10


#p_f

pvals_sp2 <- read_csv(
  find_file("Supplementary_Table_23_DA_LMM_emmeans_results_disease.csv")
) %>%
  separate(Comparison, into = c("region1", "region2"), sep = " - ")

taxaNames_sp2 <- taxaNames_all[taxaNames_all$full.taxon %in% pvals_sp2$Taxon, ]
#Extract species short name (no prefix) from full taxonomy

taxaNames_sp2$species_name <- sub(".*\\.s__", "", taxaNames_sp2$full.taxon)
for (i in seq_len(nrow(taxaNames_sp2))) {
  pvals_sp2$Taxon[pvals_sp2$Taxon == taxaNames_sp2$full.taxon[i]] <- taxaNames_sp2$species_name[i]
}

valid_sp2    <- colnames(species_rel_base) %in% taxaNames_sp2$full.taxon
species_rel2 <- species_rel_base[, valid_sp2, drop = FALSE]
ord_sp2      <- match(colnames(species_rel2), taxaNames_sp2$full.taxon)
colnames(species_rel2) <- taxaNames_sp2$species_name[ord_sp2]

species_rel2 <- as.data.frame(species_rel2)
species_rel2$disease <- meta_species$disease[
  match(rownames(species_rel2), meta_species$sample)]
species_rel2$disease <- factor(species_rel2$disease, levels = c("healthy", "IBD"))

taxa_sp2_cols <- setdiff(colnames(species_rel2), "disease")
#All species from pvals (for heatmap y-axis) / pvalsy
all_taxa_sp2  <- unique(pvals_sp2$Taxon)
taxa_means_sp2 <- data.frame(
  taxon      = all_taxa_sp2,
  taxa.means = sapply(all_taxa_sp2, function(t) {
    if (t %in% taxa_sp2_cols) mean(species_rel2[[t]], na.rm = TRUE) else 0
  })
) %>% arrange(desc(taxa.means))
top10_sp2 <- taxa_means_sp2$taxon[taxa_means_sp2$taxon %in% taxa_sp2_cols][1:10]

species_long_sp2 <- species_rel2[, c(top10_sp2, "disease")] %>%
  pivot_longer(-disease, names_to = "taxon", values_to = "rel_abundance")
species_long_sp2$taxon   <- factor(species_long_sp2$taxon,   levels = top10_sp2)
species_long_sp2$disease <- factor(species_long_sp2$disease, levels = c("healthy", "IBD"))

p_f <- wrap_elements(
  ggplot(data = species_long_sp2) +
    geom_density_ridges(
      aes(x    = log10(rel_abundance),
          y    = factor(disease, levels = c("healthy", "IBD")),
          fill = disease),
      quantile_lines = TRUE, quantiles = 2
    ) +
    scale_fill_manual(values = diseasecolors_2, aesthetics = c("colour", "fill")) +
    theme(
      legend.position = "none",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      text            = element_text(color = "black"),
      plot.margin     = unit(c(10, 10, 10, 10), "pt")
    ) +
    labs(x = "Relative Abundance", y = "Density") +
    scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
    facet_wrap(~taxon, ncol = 5)
)
p_f




#Fig J / J
#Species-level heatmap composite: adj p-value + mean abundance (NC / IBD)
#p + NC / IBD


#p_j

taxa_means_sp2$taxon <- factor(taxa_means_sp2$taxon, levels = taxa_means_sp2$taxon)

mean_bar_sp2 <- ggplot(data = taxa_means_sp2) +
  geom_col(aes(x = taxa.means, y = taxon), fill = "gray30") +
  xlab("Mean Relative \nAbundance") +
  theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "pt"),
    plot.margin   = unit(c(0, 0, 0, 0), "pt"),
    axis.title.x  = element_text(size = 7),
    axis.text.x   = element_text(size = 6)
  ) +
  scale_x_continuous(breaks = c(0, 0.1)) +
  scale_y_discrete(limits = rev(levels(taxa_means_sp2$taxon)))

pval_hm_sp2 <- ggplot(data = pvals_sp2 %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = factor(region2, levels = c("healthy", "IBD")),
    y    = factor(Taxon, levels = rev(as.character(taxa_means_sp2$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 6),
    axis.title    = element_blank(), axis.ticks = element_blank(),
    legend.position = "top", legend.title = element_text(size = 9),
    text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels = function(x) paste0("s__", x)) +
  guides(fill = guide_colorbar(
    title = "Adjusted P Value", title.position = "top", title.hjust = 0.5
  ))

pval_legend_sp2 <- cowplot::get_legend(pval_hm_sp2)
pval_hm_sp2 <- pval_hm_sp2 + theme(legend.position = "none")

pval_strip_sp2 <- ggplot(
  data.frame(region = factor(c("healthy", "IBD"), levels = c("healthy", "IBD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_2) +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_pval_sp2 <- pval_hm_sp2 +
  inset_element(pval_strip_sp2, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

if (exists("counts_sp2")) rm(counts_sp2)
for (i in seq_len(nrow(regionDF_2))) {
  grp  <- regionDF_2$region[i]
  idx  <- !is.na(species_rel2$disease) & species_rel2$disease == grp
  cnts <- species_rel2[idx, taxa_sp2_cols, drop = FALSE]
  mnDF <- data.frame(
    taxon  = colnames(cnts),
    mean   = colMeans(cnts, na.rm = TRUE),
    region = grp
  )
  counts_sp2 <- if (!exists("counts_sp2")) mnDF else rbind(counts_sp2, mnDF)
}
counts_sp2$region <- factor(counts_sp2$region, levels = c("healthy", "IBD"))

mean_hm_sp2 <- ggplot(counts_sp2) +
  geom_tile(aes(
    x    = factor(region, levels = c("healthy", "IBD")),
    y    = factor(taxon,  levels = rev(as.character(taxa_means_sp2$taxon))),
    fill = log10(mean)
  )) +
  scale_fill_viridis() +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing   = unit(c(0, 0, 0, 0), "pt"),
    plot.margin     = unit(c(0, 0, 0, 0), "pt"),
    legend.position = "top",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(
    title = "Mean Relative Abundance", title.position = "top", title.hjust = 0.5
  ))

mean_legend_sp2 <- cowplot::get_legend(mean_hm_sp2)
mean_hm_sp2 <- mean_hm_sp2 + theme(legend.position = "none")

mean_strip_sp2 <- ggplot(
  data.frame(region = factor(c("healthy", "IBD"), levels = c("healthy", "IBD")))
) +
  geom_tile(aes(x = region, y = 1, fill = region)) +
  scale_fill_manual(values = diseasecolors_2) +
  scale_x_discrete(limits = c("healthy", "IBD"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))

p_mean_sp2 <- mean_hm_sp2 +
  inset_element(mean_strip_sp2, 0, -0.02, 1, 0,
                align_to = "panel", ignore_tag = TRUE)

bars_sp2 <- wrap_elements(p_pval_sp2 + p_mean_sp2 + mean_bar_sp2)
blank_j  <- ggplot() + theme_void()

p_j <- blank_j + bars_sp2 +
  inset_element(pval_legend_sp2, 0.05, 1.01, 0.55, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  inset_element(mean_legend_sp2, 0.56, 1.01, 1.00, 1.06,
                align_to = "plot", ignore_tag = TRUE) +
  plot_layout(ncol = 1)
p_j




#Save individual panels A–J / A–J


ggsave(file.path(output_dir, "Figure4A_genus_bubble.pdf"),
       plot = p_a, width = 5, height = 5, device = "pdf")
ggsave(file.path(output_dir, "Figure4B_species_bubble.pdf"),
       plot = p_b, width = 5, height = 5, device = "pdf")
ggsave(file.path(output_dir, "Figure4C_genus_ridges_3groups.pdf"),
       plot = p_c, width = 10, height = 5, device = "pdf")
ggsave(file.path(output_dir, "Figure4D_species_ridges_3groups.pdf"),
       plot = p_d, width = 10, height = 5, device = "pdf")
ggsave(file.path(output_dir, "Figure4E_genus_ridges_2groups.pdf"),
       plot = p_e, width = 10, height = 4, device = "pdf")
ggsave(file.path(output_dir, "Figure4F_species_ridges_2groups.pdf"),
       plot = p_f, width = 10, height = 4, device = "pdf")
ggsave(file.path(output_dir, "Figure4G_genus_heatmap_3groups.pdf"),
       plot = p_g, width = 8, height = 12, device = "pdf")
ggsave(file.path(output_dir, "Figure4H_genus_heatmap_2groups.pdf"),
       plot = p_h, width = 6, height = 12, device = "pdf")
ggsave(file.path(output_dir, "Figure4I_species_heatmap_3groups.pdf"),
       plot = p_i, width = 8, height = 12, device = "pdf")
ggsave(file.path(output_dir, "Figure4J_species_heatmap_2groups.pdf"),
       plot = p_j, width = 6, height = 12, device = "pdf")




#Figure 4 (A–J) / 4A–J
#Layout / :
#Row 1: A B (bubble charts)
#Row 2: C D (ridges 3-group)
#Row 3: E F (ridges 2-group)
#Row 4: G H I J (heatmap composites)


fig4_layout <- "
AABB
CCDD
EEFF
GHIJ
"

fig4_combined <- p_a + p_b +
  p_c + p_d +
  p_e + p_f +
  p_g + p_h + p_i + p_j +
  plot_layout(design = fig4_layout) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(output_dir, "Figure4_combined_A-J.pdf"),
  plot   = fig4_combined,
  width  = 24, height = 36, dpi = 300, device = "pdf"
)

message("Figure 4 panels A–J saved to: ", output_dir)




#Figures /
#full pairwise heatmaps for genus and species, extended ridges)



#Supp: Full genus pairwise heatmap (NC / UC / CD all pairs)
#NC / UC / CD
supp_g3_full <- ggplot(data = pvals_g3 %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = region2,
    y    = factor(Taxon, levels = rev(as.character(taxa_means_g3$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.title    = element_blank(), axis.ticks = element_blank(),
    legend.position = "top", text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title = "Adjusted P Value (all pairwise)",
    title.position = "top", title.hjust = 0.5
  )) +
  labs(title = "Supplementary: Genus-level pairwise DA (NC/UC/CD)")

ggsave(
  file.path(output_dir, "Supp_Figure_genus_pairwise_heatmap_3groups.pdf"),
  plot = supp_g3_full, width = 8, height = 14, dpi = 300
)

#Supp: Full species pairwise heatmap (NC / UC / CD all pairs)
#NC / UC / CD
supp_sp3_full <- ggplot(data = pvals_sp3 %>% filter(!is.na(Taxon))) +
  geom_tile(aes(
    x    = region2,
    y    = factor(Taxon, levels = rev(as.character(taxa_means_sp3$taxon))),
    fill = log10(AdjP_BH)
  )) +
  scale_fill_gradient2(low = "red", high = "white") +
  theme(
    axis.text.x   = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.title    = element_blank(), axis.ticks = element_blank(),
    legend.position = "top", text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title = "Adjusted P Value (all pairwise)",
    title.position = "top", title.hjust = 0.5
  )) +
  labs(title = "Supplementary: Species-level pairwise DA (NC/UC/CD)")

ggsave(
  file.path(output_dir, "Supp_Figure_species_pairwise_heatmap_3groups.pdf"),
  plot = supp_sp3_full, width = 8, height = 14, dpi = 300
)

#Supp: All genus ridges NC / UC / CD (extended beyond top 10)
#NC / UC / CD10
genus_long_g3_all <- genus_rel3[, c(as.character(taxa_means_g3$taxon), "disease_subtype")] %>%
  pivot_longer(-disease_subtype, names_to = "taxon", values_to = "rel_abundance")
genus_long_g3_all$taxon <- factor(genus_long_g3_all$taxon,
                                   levels = as.character(taxa_means_g3$taxon))

supp_ridges_g3 <- ggplot(data = genus_long_g3_all) +
  geom_density_ridges(
    aes(x    = log10(rel_abundance),
        y    = factor(disease_subtype, levels = c("NC", "UC", "CD")),
        fill = disease_subtype),
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    text            = element_text(color = "black")
  ) +
  labs(x = "Relative Abundance", y = "Density",
       title = "Supplementary: All genus ridges (NC/UC/CD)") +
  scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
  facet_wrap(~taxon, ncol = 5)

ggsave(
  file.path(output_dir, "Supp_Figure_genus_all_ridges_3groups.pdf"),
  plot = supp_ridges_g3, width = 14, height = 20, dpi = 300
)

#Supp: All species ridges NC / UC / CD (extended beyond top 10)
#NC / UC / CD10
species_long_sp3_all <- species_rel3[, c(as.character(taxa_means_sp3$taxon), "disease_subtype")] %>%
  pivot_longer(-disease_subtype, names_to = "taxon", values_to = "rel_abundance")
species_long_sp3_all$taxon <- factor(species_long_sp3_all$taxon,
                                      levels = as.character(taxa_means_sp3$taxon))

supp_ridges_sp3 <- ggplot(data = species_long_sp3_all) +
  geom_density_ridges(
    aes(x    = log10(rel_abundance),
        y    = factor(disease_subtype, levels = c("NC", "UC", "CD")),
        fill = disease_subtype),
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_fill_manual(values = diseasecolors_3, aesthetics = c("colour", "fill")) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    text            = element_text(color = "black")
  ) +
  labs(x = "Relative Abundance", y = "Density",
       title = "Supplementary: All species ridges (NC/UC/CD)") +
  scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
  facet_wrap(~taxon, ncol = 5)

ggsave(
  file.path(output_dir, "Supp_Figure_species_all_ridges_3groups.pdf"),
  plot = supp_ridges_sp3, width = 14, height = 20, dpi = 300
)

#Supp: All genus ridges NC / IBD (extended)
#NC / IBD
genus_long_g2_all <- genus_rel2[, c(as.character(taxa_means_g2$taxon), "disease")] %>%
  pivot_longer(-disease, names_to = "taxon", values_to = "rel_abundance")
genus_long_g2_all$taxon   <- factor(genus_long_g2_all$taxon,
                                     levels = as.character(taxa_means_g2$taxon))
genus_long_g2_all$disease <- factor(genus_long_g2_all$disease,
                                     levels = c("healthy", "IBD"))

supp_ridges_g2 <- ggplot(data = genus_long_g2_all) +
  geom_density_ridges(
    aes(x    = log10(rel_abundance),
        y    = factor(disease, levels = c("healthy", "IBD")),
        fill = disease),
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_fill_manual(values = diseasecolors_2, aesthetics = c("colour", "fill")) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    text            = element_text(color = "black")
  ) +
  labs(x = "Relative Abundance", y = "Density",
       title = "Supplementary: All genus ridges (NC/IBD)") +
  scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
  facet_wrap(~taxon, ncol = 5)

ggsave(
  file.path(output_dir, "Supp_Figure_genus_all_ridges_2groups.pdf"),
  plot = supp_ridges_g2, width = 14, height = 18, dpi = 300
)

#Supp: All species ridges NC / IBD (extended)
#NC / IBD
species_long_sp2_all <- species_rel2[, c(as.character(taxa_means_sp2$taxon), "disease")] %>%
  pivot_longer(-disease, names_to = "taxon", values_to = "rel_abundance")
species_long_sp2_all$taxon   <- factor(species_long_sp2_all$taxon,
                                        levels = as.character(taxa_means_sp2$taxon))
species_long_sp2_all$disease <- factor(species_long_sp2_all$disease,
                                        levels = c("healthy", "IBD"))

supp_ridges_sp2 <- ggplot(data = species_long_sp2_all) +
  geom_density_ridges(
    aes(x    = log10(rel_abundance),
        y    = factor(disease, levels = c("healthy", "IBD")),
        fill = disease),
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_fill_manual(values = diseasecolors_2, aesthetics = c("colour", "fill")) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    text            = element_text(color = "black")
  ) +
  labs(x = "Relative Abundance", y = "Density",
       title = "Supplementary: All species ridges (NC/IBD)") +
  scale_x_continuous(labels = scales::math_format(), breaks = c(-5, -1)) +
  facet_wrap(~taxon, ncol = 5)

ggsave(
  file.path(output_dir, "Supp_Figure_species_all_ridges_2groups.pdf"),
  plot = supp_ridges_sp2, width = 14, height = 18, dpi = 300
)

message("All Figure 4 panels (A–J) and supplementary figures saved to: ", output_dir)
