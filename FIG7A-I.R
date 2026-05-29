#Figure 7
#FigA-I
#Supplementary tables

required_pkgs <- c(
  "dplyr", "readr", "readxl", "tidyr", "tibble", "ggplot2", "stringr",
  "ape", "phytools", "ggtree", "ggtreeExtra", "ggnewscale",
  "colorspace", "cowplot", "scales", "lmerTest",
  "fgsea", "ggrepel", "ggforce", "ggbeeswarm", "magick"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- if (!is.na(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg))) else normalizePath(getwd())
output_dir <- script_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

find_file <- function(filename) {
  candidates <- c(
    file.path(script_dir, filename),
    file.path(dirname(script_dir), filename),
    file.path("D:/R代码", filename),
    file.path("D:/R代码/发育树代码-cell/Phylogenies_AndreSanchezetal", filename),
    file.path("D:/R代码/发育树代码-cell/Results", filename),
    file.path("E:/2TE盘备份/IBS/iHMP_IBDMDB_2019", filename),
    file.path("E:/2TE盘备份/Results3", filename)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Required file not found: ", filename)
  hit
}

tree_dir <- "D:/R代码/发育树代码-cell/Phylogenies_AndreSanchezetal"
stage_colors <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")
ibd_colors <- c("NC" = "#bdd7e7", "IBD" = "#FF7F7F")

iso3_to_iso2 <- tibble::tribble(
  ~Country3, ~iso,
  "ITA","IT","USA","US","GBR","GB","SWE","SE","DEU","DE","FJI","FJ","CHN","CN","KAZ","KZ",
  "BGD","BD","IND","IN","AUT","AT","CAN","CA","DNK","DK","LUX","LU","FRA","FR","NOR","NO",
  "SVK","SK","HUN","HU","EST","EE","FIN","FI","ISL","IS","ESP","ES","IRL","IE","NLD","NL",
  "MNG","MN","CMR","CM","JPN","JP","PER","PE","MDG","MG","SLV","SV","CHE","CH","TZA","TZ",
  "IDN","ID","LBR","LR","GHA","GH","ETH","ET","ARG","AR","COL","CO","GNB","GW","RUS","RU",
  "KOR","KR","ISR","IL"
)

clean_tip_names <- function(x) {
  x <- as.character(x)
  x <- vapply(strsplit(x, "_metaphlan4", fixed = TRUE), `[`, character(1), 1)
  x <- vapply(strsplit(x, "_"), function(parts) {
    if (any(grepl("^HV", parts))) parts[1] else paste(parts, collapse = "_")
  }, character(1))
  x
}

make_country_legends <- function(iso_values) {
  iso_values <- sort(unique(stats::na.omit(iso_values)))
  country_colors <- colorspace::qualitative_hcl(
    n = length(iso_values), palette = "Dark 3", h = c(0, 360), c = 60, l = 65
  )
  names(country_colors) <- iso_values
  df_leg <- data.frame(iso = factor(iso_values, levels = iso_values), x = 1, y = 1)

  p_leg_h <- ggplot(df_leg, aes(x, y, fill = iso)) +
    geom_tile() +
    scale_fill_manual(values = country_colors, name = "Country") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"))

  p_leg_v <- ggplot(df_leg, aes(x, y, fill = iso)) +
    geom_tile() +
    scale_fill_manual(values = country_colors, name = "Country") +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE, title.position = "top")) +
    theme_void() +
    theme(legend.position = "right", legend.title = element_text(size = 12, face = "bold"))

  leg_h <- cowplot::ggdraw(cowplot::get_legend(p_leg_h))
  leg_v <- cowplot::get_legend(p_leg_v)
  ggsave(file.path(output_dir, "Figure7H_country_legend_horizontal.pdf"), leg_h, width = 8, height = 2)
  ggsave(file.path(output_dir, "Figure7H_country_legend_vertical.pdf"), leg_v, width = 2.5, height = 5)
  country_colors
}

#FigA-D
iso_list <- c(
  "AT","BD","CA","CH","CM","CN","CO","DE","DK","EE","ES","FI","FR","GB","IE","IL","IN",
  "IT","JP","KR","KZ","MG","NL","PE","RU","SE","TZ","US","NO","HU","IS","GH"
)
expanded_colors <- colorspace::qualitative_hcl(
  n = length(iso_list),
  palette = "Dark 3",
  h = c(0, 360),
  c = 60,
  l = 65
)
names(expanded_colors) <- iso_list

Data_fixed <- readr::read_tsv(find_file("Phenotypes_merged.tsv"), show_col_types = FALSE) %>%
  dplyr::filter(!study_name %in% c("ThomasAM_2019_c", "500FG_FSK")) %>%
  dplyr::filter(!is.na(stage), !is.na(iso))

run_tree_panel <- function(tree_file, data_use = Data_fixed, perm = 999) {
  tree <- ape::read.tree(tree_file)
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  tree <- phytools::midpoint.root(tree)
  tree$tip.label <- vapply(tree$tip.label, clean_tip_names, character(1))

  data_local <- data_use %>%
    dplyr::select(ID_anal, subject_id, iso, stage) %>%
    dplyr::distinct(subject_id, .keep_all = TRUE) %>%
    dplyr::filter(ID_anal %in% tree$tip.label) %>%
    dplyr::distinct(ID_anal, .keep_all = TRUE)

  tree2 <- ape::keep.tip(tree, data_local$ID_anal)
  data_local <- data_local[match(tree2$tip.label, data_local$ID_anal), ]

  distance_tree <- ape::cophenetic.phylo(tree2)
  stage_dist <- stats::dist(as.numeric(data_local$stage), method = "euclidean")
  test_stage <- vegan::mantel(as.dist(distance_tree), stage_dist, method = "pearson", permutations = perm)

  result <- tibble::tibble(
    Tree_file = basename(tree_file),
    Rho = unname(test_stage$statistic),
    P = test_stage$signif,
    Permutations = perm
  )

  plot_obj <- ggtree::ggtree(tree2, layout = "circular") %<+%
    (data_local %>% dplyr::rename(label = ID_anal)) +
    ggtree::geom_tippoint(ggplot2::aes(color = factor(stage)), size = 1) +
    ggplot2::scale_color_manual(values = stage_colors, name = "stage") +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      geom = "geom_tile",
      mapping = ggplot2::aes(fill = iso),
      width = 0.015,
      offset = 0.10
    ) +
    ggplot2::scale_fill_manual(values = expanded_colors, name = "country")

  list(result = result, plot = plot_obj, data_local = data_local)
}

fig_a <- run_tree_panel(find_file("IQtree.t__SGB14546_group.TreeShrink.tre"))
p_a <- fig_a$plot
p_a
ggsave(file.path(output_dir, "Figure7A_tree.pdf"), p_a, width = 10, height = 10)

fig_b <- run_tree_panel(find_file("RAxML_bestTree.t__SGB15323.TreeShrink.tre"))
p_b <- fig_b$plot
p_b
ggsave(file.path(output_dir, "Figure7B_tree.pdf"), p_b, width = 10, height = 10)

fig_c <- run_tree_panel(find_file("RAxML_bestTree.t__SGB4837_group.TreeShrink.tre"))
p_c <- fig_c$plot
p_c
ggsave(file.path(output_dir, "Figure7C_tree.pdf"), p_c, width = 10, height = 10)

fig_d <- run_tree_panel(find_file("RAxML_bestTree.t__SGB2301.TreeShrink.tre"))
p_d <- fig_d$plot
p_d
ggsave(file.path(output_dir, "Figure7D_tree.pdf"), p_d, width = 10, height = 10)

readr::write_tsv(
  dplyr::bind_rows(
    fig_a$result %>% dplyr::mutate(panel = "A"),
    fig_b$result %>% dplyr::mutate(panel = "B"),
    fig_c$result %>% dplyr::mutate(panel = "C"),
    fig_d$result %>% dplyr::mutate(panel = "D")
  ),
  file.path(output_dir, "Figure7A-D_source_mapping.tsv")
)

#FigE
sgb_taxonomy <- readr::read_csv(find_file("species.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    SGB = sub("_group$", "", SGB),
    Taxonomy = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = "|")
  ) %>%
  tidyr::separate(
    Taxonomy,
    into = c("domain", "phylum", "class", "order", "family", "genera", "species"),
    sep = "\\|",
    remove = FALSE
  )

geography_assocations <- readr::read_csv(
  find_file("Supplementary_Table_32_Mantel_stage_ALL_combined_1659.csv"),
  show_col_types = FALSE
) %>%
  dplyr::left_join(sgb_taxonomy, by = "SGB")

taxa_list <- list()
for (level in c("phylum", "class", "order", "family", "genera")) {
  taxa_list <- c(taxa_list, split(geography_assocations$SGB, geography_assocations[[level]]))
}

r1_clean <- geography_assocations %>%
  dplyr::filter(is.finite(Rho)) %>%
  dplyr::arrange(dplyr::desc(Rho))
ranks <- r1_clean$Rho
names(ranks) <- as.character(r1_clean$SGB)

fgsea_res <- fgsea::fgsea(taxa_list, ranks, eps = 0)
assoc_rho_taxa <- tibble::as_tibble(fgsea_res) %>% dplyr::arrange(padj)
assoc_rho_taxa2 <- assoc_rho_taxa %>%
  dplyr::mutate(leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ",")))
readr::write_tsv(assoc_rho_taxa2, file.path(output_dir, "Supplementary_Table_34_Geography_enrichment.tsv"))

p_e <- fgsea::plotEnrichment(taxa_list[["g__Lachnospiraceae_unclassified"]], ranks) +
  ggplot2::labs(title = "g__Lachnospiraceae_unclassified") +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
assign("p_e", p_e, envir = .GlobalEnv)
p_e
ggplot2::ggsave(file.path(output_dir, "Figure7E_Lachnospiraceae_enrichment.pdf"), p_e, width = 5, height = 4, dpi = 300)

#FigF-G
traitar <- read.table(
  find_file("traitar_merged_1655SGBs_traits.txt"),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
traitar <- as.data.frame(lapply(traitar, function(col) if (is.character(col)) enc2utf8(col) else col))
traitar <- traitar %>%
  dplyr::rename(SGB_clean = SGB) %>%
  dplyr::mutate(SGB_clean = stringr::str_replace(SGB_clean, "_group$", ""))

geo_effect <- readxl::read_excel(find_file("Supplementary_Table_32_Mantel_stage_ALL_combined_1659.xls")) %>%
  dplyr::mutate(
    SGB_clean = gsub("^t__", "", SGB),
    SGB_clean = sub("_group$", "", SGB_clean),
    Rho = as.numeric(Rho),
    P = as.numeric(P),
    FDR = as.numeric(FDR)
  )

sgb_taxonomy2 <- readr::read_csv(find_file("species.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(SGB = sub("_group$", "", SGB))

traitar_merged <- dplyr::left_join(geo_effect, traitar, by = "SGB_clean") %>%
  tidyr::drop_na() %>%
  dplyr::filter(!is.na(Rho), is.finite(Rho), N > 50) %>%
  dplyr::left_join(
    sgb_taxonomy2 %>% dplyr::select(SGB, Genus, Family),
    by = c("SGB_clean" = "SGB")
  )

traitar_results <- tibble::tibble()
for (ann in colnames(traitar)) {
  if (ann == "SGB_clean") next
  merged2 <- traitar_merged
  suppressWarnings(merged2[merged2[[ann]] == 0.5, ann] <- NA)
  tbl <- table(as.vector(merged2[[ann]]))
  if (length(tbl) < 2 || min(tbl) < 2) next
  model <- tryCatch(
    stats::lm(stats::as.formula(paste0("Rho ~ `", ann, "`")), data = tidyr::drop_na(merged2)) %>% summary(),
    error = function(e) NULL
  )
  if (is.null(model) || nrow(model$coefficients) < 2) next
  traitar_results <- dplyr::bind_rows(
    traitar_results,
    tibble::as_tibble(t(model$coefficients[2, ])) %>% dplyr::mutate(Annotation = ann)
  )
}
traitar_results <- traitar_results %>%
  dplyr::arrange(`Pr(>|t|)`) %>%
  dplyr::mutate(FDR = p.adjust(`Pr(>|t|)`, method = "BH"))
readr::write_tsv(traitar_results, file.path(output_dir, "Supplementary_Table_35_Traitar_results.tsv"))
readr::write_tsv(traitar_merged, file.path(output_dir, "Traitar_merged_Figure7FG.tsv"))

p_f <- traitar_results %>%
  ggplot2::ggplot(ggplot2::aes(x = Estimate, y = -log10(`Pr(>|t|)`), color = FDR < 0.05)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_bw() +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "#d62728")) +
  ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(FDR < 0.05, Annotation, "")), box.padding = 0.4, point.padding = 0.4, size = 4) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 16), axis.text = ggplot2::element_text(size = 14)) +
  ggplot2::labs(y = "-log10(p)", x = "Effect size (Estimate)")
assign("p_f", p_f, envir = .GlobalEnv)
p_f
ggplot2::ggsave(file.path(output_dir, "Figure7F_Traitar_volcano.pdf"), p_f, width = 5, height = 4, dpi = 300)

to_plot <- traitar_results %>%
  dplyr::arrange(`Pr(>|t|)`) %>%
  dplyr::slice(1:3) %>%
  dplyr::pull(Annotation)

p_g <- traitar_merged %>%
  dplyr::select(Rho, SGB_clean, Genus, dplyr::all_of(to_plot)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(to_plot), names_to = "Trait", values_to = "Available") %>%
  dplyr::filter(Available != 0.5) %>%
  dplyr::mutate(
    Available = ifelse(Available == 1, "Yes", "No"),
    Lachno = ifelse(Genus == "g__Lachnospiraceae_unclassified", "TRUE", "FALSE")
  ) %>%
  tidyr::drop_na(Rho, Available, Lachno) %>%
  ggplot2::ggplot(ggplot2::aes(x = Available, y = Rho, color = Lachno, fill = Lachno)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.35, position = ggplot2::position_dodge(width = 0.8)) +
  ggforce::geom_sina(alpha = 0.7, size = 1.6, position = ggplot2::position_dodge(width = 0.8)) +
  ggplot2::facet_wrap(~Trait, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Predicted trait", y = "Rho") +
  ggplot2::scale_color_manual(values = c("FALSE" = "#6A3D9A", "TRUE" = "#FF7F00")) +
  ggplot2::scale_fill_manual(values = c("FALSE" = "#6A3D9A50", "TRUE" = "#FF7F0050")) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16, face = "bold"),
    axis.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11)
  )
assign("p_g", p_g, envir = .GlobalEnv)
ggplot2::ggsave(file.path(output_dir, "Figure7G_Traitar_boxplots.pdf"), p_g, width = 7, height = 4, dpi = 300)
p_g
#FigH
metadata_metabolites <- readr::read_tsv(find_file("Phenotypes_with_metabolites.tsv"), show_col_types = FALSE) %>%
  dplyr::filter(!study_name %in% c("ThomasAM_2019_c", "500FG_FSK")) %>%
  dplyr::group_by(ID_anal) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(stage)) %>%
  dplyr::transmute(
    sample_id = ID_anal,
    stage = as.character(stage),
    IBD_group = dplyr::case_when(IBD == 1 ~ "IBD", IBD == 0 ~ "NC", TRUE ~ NA_character_),
    age = age,
    Country = Country,
    study = study_name
  ) %>%
  dplyr::left_join(iso3_to_iso2, by = c("Country" = "Country3"))

phylo_4993 <- readr::read_tsv(find_file("SGB4993_IBD_phylo_effect_per_sample.tsv"), show_col_types = FALSE) %>%
  dplyr::select(sample_id, phylo_effect_median)

tree_candidates <- list.files(tree_dir, pattern = "SGB4993.*\\.tre$", full.names = TRUE)
tree_file <- if (any(grepl("TreeShrink", tree_candidates))) tree_candidates[grepl("TreeShrink", tree_candidates)][1] else tree_candidates[1]
tree_h <- ape::read.tree(tree_file)
if (inherits(tree_h, "multiPhylo")) tree_h <- tree_h[[1]]
tree_h <- phytools::midpoint.root(tree_h)
tree_h$tip.label <- clean_tip_names(tree_h$tip.label)

meta_h <- metadata_metabolites %>%
  dplyr::filter(sample_id %in% tree_h$tip.label) %>%
  dplyr::filter(!is.na(IBD_group), !is.na(stage), !is.na(age), !is.na(iso)) %>%
  dplyr::left_join(phylo_4993, by = "sample_id")

tree_h <- ape::keep.tip(tree_h, meta_h$sample_id)
meta_h <- meta_h %>% dplyr::arrange(match(sample_id, tree_h$tip.label))
country_colors <- make_country_legends(meta_h$iso)

p_h <- ggtree::ggtree(tree_h, layout = "fan", open.angle = 15, size = 0.1) %<+% meta_h
p_h <- p_h +
  ggtree::geom_tippoint(aes(color = IBD_group), size = 0.25, alpha = 0.9) +
  scale_color_manual(values = ibd_colors, name = "IBD")
p_h <- p_h +
  ggtreeExtra::geom_fruit(geom = geom_tile, mapping = aes(fill = IBD_group), width = 0.05, offset = 0.05) +
  scale_fill_manual(values = ibd_colors, name = "IBD / NC")
p_h <- p_h + ggnewscale::new_scale_fill()
p_h <- p_h +
  ggtreeExtra::geom_fruit(geom = geom_tile, mapping = aes(fill = phylo_effect_median), width = 0.05, offset = 0.05) +
  scale_fill_viridis_c(option = "magma", name = "Phylo. effect\nmedian", na.value = "white")
p_h <- p_h + ggnewscale::new_scale_fill()
p_h <- p_h +
  ggtreeExtra::geom_fruit(geom = geom_tile, mapping = aes(fill = age), width = 0.05, offset = 0.05) +
  scale_fill_viridis_c(option = "viridis", name = "Age", na.value = "white")
p_h <- p_h + ggnewscale::new_scale_fill()
p_h <- p_h +
  ggtreeExtra::geom_fruit(geom = geom_tile, mapping = aes(fill = iso), width = 0.05, offset = 0.05) +
  scale_fill_manual(values = country_colors, name = "Country")
p_h <- p_h + ggnewscale::new_scale_fill()
p_h <- p_h +
  ggtreeExtra::geom_fruit(geom = geom_tile, mapping = aes(fill = stage), width = 0.05, offset = 0.05) +
  scale_fill_manual(values = stage_colors, name = "Stage")
assign("p_h", p_h, envir = .GlobalEnv)
ggsave(file.path(output_dir, "Figure7H_tree.pdf"), plot = p_h, width = 8, height = 15)
save(tree_h, meta_h, tree_file, file = file.path(output_dir, "Figure7H_panel.RData"))
p_h
#FigI
ba_vars <- c("CA","CDCA","DCA","LCA","GCA","TCA","TDCA","GUDCA","TUDCA")
cluster_xlsx <- list.files("E:/2TE盘备份/Results3", pattern = "\\.xlsx$", full.names = TRUE)[1]
df_cluster <- readxl::read_xlsx(cluster_xlsx) %>%
  dplyr::mutate(Group = as.numeric(as.character(Group)), study = as.character(study)) %>%
  dplyr::select(sample_id, Group, age, study, dplyr::all_of(ba_vars))

fit_one <- function(outcome) {
  df <- df_cluster %>%
    dplyr::select(sample_id, Group, age, study, dplyr::all_of(outcome)) %>%
    tidyr::drop_na()
  if (nrow(df) == 0) {
    return(tibble::tibble(BileAcid = outcome, n = 0, estimate = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  df <- df %>% dplyr::mutate(Group = factor(Group), study = factor(study))
  if (dplyr::n_distinct(df$study) > 1) {
    fit <- lmerTest::lmer(stats::as.formula(paste0(outcome, " ~ Group + age + (1|study)")), data = df, REML = FALSE)
    co <- summary(fit)$coefficients
  } else {
    fit <- stats::lm(stats::as.formula(paste0(outcome, " ~ Group + age")), data = df)
    co <- summary(fit)$coefficients
  }
  group_row <- grep("^Group", rownames(co))[1]
  tibble::tibble(
    BileAcid = outcome,
    n = nrow(df),
    estimate = co[group_row, "Estimate"],
    se = co[group_row, "Std. Error"],
    t = co[group_row, "t value"],
    p = co[group_row, "Pr(>|t|)"]
  )
}

lmm_res <- purrr::map_dfr(ba_vars, fit_one) %>%
  dplyr::mutate(FDR = p.adjust(p, method = "fdr")) %>%
  dplyr::arrange(FDR)
readr::write_tsv(lmm_res, file.path(output_dir, "Supplementary_Table_37_SGB4993_bileAcid_GroupCompare.tsv"))

df_long <- df_cluster %>%
  tidyr::pivot_longer(cols = dplyr::all_of(ba_vars), names_to = "BA", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(Group = factor(Group, levels = c(0, 1), labels = c("Outside clade", "In clade")), log_value = log10(value + 1))

sig_marks <- lmm_res %>%
  dplyr::filter(!is.na(FDR), FDR < 0.05) %>%
  dplyr::select(BA = BileAcid) %>%
  dplyr::left_join(df_long %>% dplyr::group_by(BA) %>% dplyr::summarise(y = max(log_value, na.rm = TRUE), .groups = "drop"), by = "BA") %>%
  dplyr::mutate(x = 1.5, label = "*")

p_i <- ggplot(df_long, aes(x = Group, y = log_value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  facet_wrap(~BA, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  ylab("log10(BA level)") +
  xlab("")
if (nrow(sig_marks) > 0) {
  p_i <- p_i + geom_text(data = sig_marks, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 8)
}
assign("p_i", p_i, envir = .GlobalEnv)
ggsave(file.path(output_dir, "Figure7I_bile_acid_boxplots.pdf"), p_i, width = 9, height = 7)
p_i
writeLines(
  c(
    "Figure7 A-D: generated directly in this master script",
    "Figure7 E-G: generated directly in this master script",
    "Figure7 H-I: generated directly in this master script"
  ),
  con = file.path(output_dir, "Figure7_run_notes.txt")
)
