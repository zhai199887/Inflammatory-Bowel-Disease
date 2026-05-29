#Figure 5
#Figure 5 complete master script (A-D + supplementary figures & tables)


#Package loading

pkgs <- c(
  "data.table", "readxl", "stringr", "dplyr", "tidyr",
  "coin", "ggplot2", "patchwork", "RColorBrewer",
  "Hmisc", "igraph", "ggraph", "scales", "ggrepel"
)
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE))


#Helper: file finder

script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd()
)
output_dir <- script_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

find_file <- function(filename) {
  candidates <- c(file.path(script_dir, filename), filename, file.path("D:/R代码", filename))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("File not found: ", filename)
  hit
}



#Figure 5A-C Data Preparation
#merged_6cohorts_common.tsv, Supplementary_Table_21__metadata_common.tsv, Table_S6.xlsx
#Input files: merged_6cohorts_common.tsv, Supplementary_Table_21__metadata_common.tsv, Table_S6.xlsx



#metadata / Read pathway wide table and sample metadata
merged_6cohorts <- fread(find_file("merged_6cohorts_common.tsv"), sep = "\t", header = TRUE)
metadata_common  <- fread(find_file("Supplementary_Table_21__metadata_common.tsv"), sep = "\t", header = TRUE)

#metadata / Intersect pathway sample columns with metadata sample column
sample_cols     <- setdiff(names(merged_6cohorts), "pathway")
common_samples  <- intersect(sample_cols, metadata_common$sample)

#Keep only intersection samples
merged_6cohorts <- merged_6cohorts[, c("pathway", common_samples), with = FALSE]
metadata_common <- metadata_common[sample %in% common_samples]

#Convert to matrix and re-close pathway subset to sum=1 per sample
mat <- as.matrix(merged_6cohorts[, ..common_samples])
colsum2 <- colSums(mat, na.rm = TRUE)
colsum2[colsum2 == 0] <- NA_real_
mat_rel_common <- sweep(mat, 2, colsum2, "/")
merged_rel <- data.table(pathway = merged_6cohorts$pathway, mat_rel_common)
setnames(merged_rel, c("pathway", common_samples))

#metadata / Pivot to long + merge metadata
pw_long <- melt(merged_rel, id.vars = "pathway", variable.name = "sample", value.name = "abund")
pw_long <- merge(pw_long, metadata_common, by = "sample", all.x = TRUE)

#NA projectabund HUMAnN /
#Basic cleaning: remove NA project; numeric abund; remove UNMAPPED / UNINTEGRATED; re-normalize within sample
pw_long <- pw_long[!is.na(project)]
pw_long[, abund := as.numeric(abund)]
pw_long[is.na(abund), abund := 0]
pw_long <- pw_long[!pathway %in% c("UNMAPPED", "UNINTEGRATED")]
pw_long[, tot := sum(abund), by = sample]
pw_long[tot > 0, abund := abund / tot]
pw_long[, tot := NULL]

#Table_S6 → / Pathway → metabolic function mapping ----------
map_raw <- as.data.table(read_excel(find_file("Table_S6.xlsx")))
setnames(
  map_raw,
  old = c("Parent of metabolic functions", "Metabolic functions", "Metaolic pathways"),
  new = c("parent_func", "met_func", "pathway_raw"),
  skip_absent = TRUE
)
map_raw[, pathway_id := sub(":.*$", "", pathway_raw)]
map_raw <- map_raw[!is.na(met_func) & met_func != ""]
map_unique <- map_raw[, .SD[1], by = pathway_id]

pw_long[, pathway_id := sub(":.*$", "", pathway)]
pw_long1 <- merge(pw_long, map_unique[, .(pathway_id, met_func)], by = "pathway_id", all.x = TRUE)
pw_long1[is.na(met_func) | met_func == "", met_func := "Others"]

#Manual extra mapping (author-defined categories) ----------
map_extra <- fread(text = "
pathway_id\tmet_func
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY\tPhenolic-Compounds-Degradation
AEROBACTINSYN-PWY\tSiderophores-Biosynthesis
ALLANTOINDEG-PWY\tAllantoin-degradation
ANAEROFRUCAT-PWY\tFermentation-of-Pyruvate
ARG+POLYAMINE-SYN\tPolyamine-Biosynthesis
ARGDEG-PWY\tAMINE-DEG
ASPASN-PWY\tAmino-Acid-Biosynthesis
BRANCHED-CHAIN-AA-SYN-PWY\tAmino-Acid-Biosynthesis
CATECHOL-ORTHO-CLEAVAGE-PWY\tCatechol-Degradation
CODH-PWY\tAcetyl-CoA-Biosynthesis
DENITRIFICATION-PWY\tANAEROBIC-RESPIRATION
DHGLUCONATE-PYR-CAT-PWY\tSugars-Degradation
FOLSYN-PWY\tVitamin-Biosynthesis
GLCMANNANAUT-PWY\tSugar-Derivatives-Degradation
HOMOSER-METSYN-PWY\tProteinogenic-Amino-Acid-Biosynthesis
LEU-DEG2-PWY\tProteinogenic-Amino-Acids-Degradation
LPSSYN-PWY\tLipopolysaccharide-Biosynthesis
LYSINE-DEG1-PWY\tProteinogenic-Amino-Acids-Degradation
METH-ACETATE-PWY\tANAEROBIC-RESPIRATION
ORNARGDEG-PWY\tProteinogenic-Amino-Acids-Degradation
P161-PWY\tEnergy-Metabolism-others
P165-PWY\tPurine-Degradation
P221-PWY\tAROMATIC-COMPOUNDS-DEGRADATION
P241-PWY\tCofactor-Biosynthesis
P562-PWY\tSugar-Derivatives-Degradation
P621-PWY\tPolysaccharide Degradation
PHOSLIPSYN-PWY\tPhospholipid-Biosynthesis
PHOTOALL-PWY\tPhotosynthesis
PROPFERM-PWY\tFermentation-to-Short-Chain-Fatty-Acids
PWY-101\tPhotosynthesis
PWY-1269\tLipopolysaccharide-Biosynthesis
PWY-1361\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-1501\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-1622\tFormaldehyde-Assimilation
PWY-2221\tGLYCOLYSIS-VARIANTS
PWY-3001\tProteinogenic-Amino-Acid-Biosynthesis
PWY-3841\tVitamin-Biosynthesis
PWY-4702\tPhosphorus-Compounds
PWY-4722\tCreatinine-Degradation
PWY-5005\tVitamin-Biosynthesis
PWY-5028\tProteinogenic-Amino-Acids-Degradation
PWY-5055\tNAD-Metabolism
PWY-5079\tProteinogenic-Amino-Acids-Degradation
PWY-5103\tProteinogenic-Amino-Acid-Biosynthesis
PWY-5156\tFatty-acid-biosynthesis
PWY-5178\tToluene-Degradation
PWY-5181\tToluene-Degradation
PWY-5183\tToluene-Degradation
PWY-5392\tCO2-Fixation
PWY-5417\tCatechol-Degradation
PWY-5419\tCatechol-Degradation
PWY-5420\tCatechol-Degradation
PWY-5431\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-5509\tVitamin-Biosynthesis
PWY-5514\tSugar-Biosynthesis
PWY-5531\tTetrapyrrole-Biosynthesis
PWY-5532\tPurine-Degradation
PWY-5534\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-561\tSuper-Pathways
PWY-5647\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-5654\tProteinogenic-Amino-Acids-Degradation
PWY-5751\tSECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5754\t4-Hydroxybenzoate-Biosynthesis
PWY-5870\tQuinone-Biosynthesis
PWY-5994\tFatty-acid-biosynthesis
PWY-6071\tAMINE-DEG
PWY-6107\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-6126\tPurine-Nucleotide-Biosynthesis
PWY-6143\tSugar-Biosynthesis
PWY-6148\tCofactor-Biosynthesis
PWY-6185\tCatechol-Degradation
PWY-6210\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-6215\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-6269\tVitamin-Biosynthesis
PWY-6281\tOther-Amino-Acid-Biosynthesis
PWY-6307\tProteinogenic-Amino-Acids-Degradation
PWY-6309\tProteinogenic-Amino-Acids-Degradation
PWY-6313\tAMINE-DEG
PWY-6318\tProteinogenic-Amino-Acids-Degradation
PWY-6344\tProteinogenic-Amino-Acids-Degradation
PWY-6349\tPhospholipid-Biosynthesis
PWY-6396\tButanediol-Biosynthesis
PWY-6467\tLipopolysaccharide-Biosynthesis
PWY-6470\tCell-Wall-Biosynthesis
PWY-6478\tSugar-Biosynthesis
PWY-6538\tPurine-Degradation
PWY-6562\tPolyamine-Biosynthesis
PWY-6565\tPolyamine-Biosynthesis
PWY-6572\tPolysaccharide Degradation
PWY-6612\tVitamin-Biosynthesis
PWY-6641\tSulfur-Metabolism
PWY-6660\tSECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6662\tSECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6748\tANAEROBIC-RESPIRATION
PWY-6749\tSugar-Biosynthesis
PWY-6785\tEnergy-Metabolism-others
PWY-6797\tVitamin-Biosynthesis
PWY-6834\tPolyamine-Biosynthesis
PWY-6837\tFatty-Acid-Degradation
PWY-6895\tVitamin-Biosynthesis
PWY-6901\tSugars-Degradation
PWY-6906\tPolysaccharide Degradation
PWY-6953\tSugar-Biosynthesis
PWY-6957\tAROMATIC-COMPOUNDS-DEGRADATION
PWY-6992\tSugars-Degradation
PWY-7031\tGlycan-Biosynthesis
PWY-7039\tLipid-Biosynthesis
PWY-7046\tPhenolic-Compounds-Degradation
PWY-7090\tSugar-Biosynthesis
PWY-7094\tFatty-acid-biosynthesis
PWY-7118\tFermentation-to-Alcohols
PWY-7124\tSECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7159\tTetrapyrrole-Biosynthesis
PWY-7165\tVitamin-Biosynthesis
PWY-7184\tDeoxyribonucleotide-Biosynthesis
PWY-7197\tDeoxyribonucleotide-Biosynthesis
PWY-7200\tPyrimidine-Nucleotide-Biosynthesis
PWY-7204\tVitamin-Biosynthesis
PWY-7208\tPyrimidine-Nucleotide-Biosynthesis
PWY-7218\tPhotosynthesis
PWY-722\tNAD-Metabolism
PWY-7245\tNAD-Metabolism
PWY-7268\tNAD-Metabolism
PWY-7288\tFatty-Acid-Degradation
PWY-7290\tSugar-Biosynthesis
PWY-7312\tSugar-Biosynthesis
PWY-7316\tSugar-Biosynthesis
PWY-7317\tSugar-Biosynthesis
PWY-7332\tSugar-Biosynthesis
PWY-7373\tQuinone-Biosynthesis
PWY-7385\tFermentation-to-Alcohols
PWY-7399\tPhosphorus-Compounds
PWY-7413\tSugar-Biosynthesis
PWY-7420\tLipid-Biosynthesis
PWY-7431\tAMINE-DEG
PWY-7528\tProteinogenic-Amino-Acid-Biosynthesis
PWY-7616\tFormaldehyde-Oxidation
PWY-841\tPurine-Nucleotide-Biosynthesis
PWY0-166\tDeoxyribonucleotide-Biosynthesis
PWY0-321\tAROMATIC-COMPOUNDS-DEGRADATION
PWY0-881\tFatty-acid-biosynthesis
PWY3O-1109\t4-Hydroxybenzoate-Biosynthesis
PWY3O-355\tFatty-acid-biosynthesis
PWY5F9-12\tAROMATIC-COMPOUNDS-DEGRADATION
PWY66-388\tFatty-Acid-Degradation
PWY66-391\tFatty-Acid-Degradation
REDCITCYC\tTCA-VARIANTS
RHAMCAT-PWY\tSugars-Degradation
THRESYN-PWY\tProteinogenic-Amino-Acid-Biosynthesis
TYRFUMCAT-PWY\tProteinogenic-Amino-Acids-Degradation
UDPNACETYLGALSYN-PWY\tSugar-Biosynthesis
URSIN-PWY\tNitrogen-Compound-Metabolism
")

pw_long1[, pathway_id := as.character(pathway_id)]
map_extra[, pathway_id := as.character(pathway_id)]
#Others / Replace only entries still labeled Others
pw_long1[map_extra, on = "pathway_id",
         met_func := fifelse(met_func == "Others" & !is.na(i.met_func), i.met_func, met_func)]

pw_long <- pw_long1#pw_long / Use pw_long uniformly from here on

#Aggregate to functional level ----------
func_long <- pw_long[, .(abund = sum(abund)), by = .(sample, project, disease_subtype, disease, met_func)]

#Statistical helper functions ----------
#blocked Wilcoxon project / Blocked Wilcoxon test (stratified by project)
blocked_wilcox_p <- function(dt_one_feature) {
  as.numeric(coin::pvalue(wilcox_test(abund ~ group2 | block,
                                      data = dt_one_feature,
                                      distribution = "approximate")))
}

#project / Prepare test data (keep only projects with both groups)
prepare_blocked_data <- function(dt, group_col, levels_keep, block_col = "project") {
  out <- copy(dt)[get(group_col) %in% levels_keep]
  out[, group2 := factor(get(group_col), levels = levels_keep)]
  out[, block  := factor(get(block_col))]
  out <- out[!is.na(block)]
  valid_blocks <- out[, uniqueN(group2), by = block][V1 == length(levels_keep), block]
  out[block %in% valid_blocks]
}

#Color ramp helper
make_pal <- function(n, base_colors) {
  if (n == 0) return(character(0))
  grDevices::colorRampPalette(base_colors)(n)
}



#Figure 5A
#NC / CD / UC /
#Stacked bar plot: metabolic function relative abundance in NC / CD / UC
#p_a / Variable: p_a



#NC / CD / UC NC vs IBD(CD+UC) blocked Wilcoxon /
#Keep NC / CD / UC; use NC vs IBD(CD+UC) blocked Wilcoxon to screen significant functions
func_subtype <- func_long[disease_subtype %in% c("NC", "CD", "UC")]
func_subtype[, group_for_test := ifelse(disease_subtype == "NC", "NC", "IBD")]

dt_test <- prepare_blocked_data(
  dt         = func_subtype,
  group_col  = "group_for_test",
  levels_keep = c("NC", "IBD"),
  block_col  = "project"
)

stat_4a <- dt_test[, .(
  mean_all = mean(abund),
  prev_all = mean(abund > 0),
  mean_NC  = mean(abund[group2 == "NC"]),
  mean_IBD = mean(abund[group2 == "IBD"]),
  pval = blocked_wilcox_p(.SD)
), by = met_func]
stat_4a[, fdr := p.adjust(pval, method = "BH")]

#FDR<0.05 & >1% & >=5% /
#Significance threshold: FDR<0.05 & mean abundance >1% & prevalence >=5%
keep_funcs <- stat_4a[fdr < 0.05 & mean_all > 0.01 & prev_all >= 0.05, met_func]

#Others / Assemble plot data: non-significant merged into Others
plot_4a_subtype <- func_subtype[, .(mean_abund = mean(abund)), by = .(disease_subtype, met_func)]
plot_4a_subtype[, met_func_plot := ifelse(met_func %in% keep_funcs, met_func, "Others")]
plot_4a_subtype <- plot_4a_subtype[, .(mean_abund = sum(mean_abund)), by = .(disease_subtype, met_func_plot)]

#NC→IBD→ / Function direction (NC-higher=blue, IBD-higher=red)
dir_tbl <- stat_4a[met_func %in% keep_funcs, .(met_func, direction = ifelse(mean_IBD > mean_NC, "IBD_higher", "NC_higher"))]
dir_tbl_no_others <- dir_tbl[met_func != "Others"]

blue_funcs <- stat_4a[
  met_func %in% dir_tbl_no_others[direction == "NC_higher", met_func],
  .(met_func, diff = mean_NC - mean_IBD)][order(diff)][, met_func]

red_funcs <- stat_4a[
  met_func %in% dir_tbl_no_others[direction == "IBD_higher", met_func],
  .(met_func, diff = mean_IBD - mean_NC)][order(diff)][, met_func]

#NC) + (IBD) + (Others) /
#Blue ramp (NC-higher) + Red ramp (IBD-higher) + Grey (Others)
cols <- c(
  setNames(make_pal(length(blue_funcs), c("#dbe9f6", "#08306b")), blue_funcs),
  setNames(make_pal(length(red_funcs),  c("#fde0dd", "#67000d")), red_funcs),
  Others = "grey80"
)
func_levels <- c(rev(red_funcs), blue_funcs, "Others")

plot_4a_subtype[, disease_subtype := factor(disease_subtype, levels = c("NC", "CD", "UC"))]
plot_4a_subtype[, met_func_plot   := factor(met_func_plot, levels = func_levels)]

#Figure 5A
p_a <- ggplot(plot_4a_subtype, aes(x = disease_subtype, y = mean_abund, fill = met_func_plot)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols, name = "Metabolic functions") +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(colour = "black", linewidth = 0.6))) +
  labs(y = "Relative abundance", x = NULL) +
  theme_classic() +
  theme(
    legend.key.size = unit(0.45, "cm"),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 10),
    legend.key      = element_rect(fill = NA, colour = NA)
  )

p_a#View



#Figure 5B
#NC vs IBD CD+UC /
#Stacked bar plot: NC vs IBD (pooled CD+UC) metabolic function relative abundance
#p_b / Variable: p_b



#healthy / IBD healthy NC /
#Keep healthy / IBD groups; relabel healthy as NC
func_disease <- func_long[disease %in% c("healthy", "IBD")]
func_disease[, disease := ifelse(disease == "healthy", "NC", "IBD")]

dt_test2 <- prepare_blocked_data(
  dt         = func_disease,
  group_col  = "disease",
  levels_keep = c("NC", "IBD"),
  block_col  = "project"
)

stat_4a2 <- dt_test2[, .(
  mean_all = mean(abund),
  prev_all = mean(abund > 0),
  mean_H   = mean(abund[group2 == "NC"]),
  mean_I   = mean(abund[group2 == "IBD"]),
  pval = blocked_wilcox_p(.SD)
), by = met_func]
stat_4a2[, fdr := p.adjust(pval, method = "BH")]

keep_funcs2 <- stat_4a2[fdr < 0.05 & mean_all > 0.01 & prev_all >= 0.05, met_func]

plot_4a_disease <- func_disease[, .(mean_abund = mean(abund)), by = .(disease, met_func)]
plot_4a_disease[, met_func_plot := ifelse(met_func %in% keep_funcs2, met_func, "Others")]
plot_4a_disease <- plot_4a_disease[, .(mean_abund = sum(mean_abund)), by = .(disease, met_func_plot)]

dir_tbl2 <- stat_4a2[met_func %in% keep_funcs2,
                     .(met_func, direction = ifelse(mean_I > mean_H, "IBD_higher", "NC_higher"))]
dir_tbl2_no_others <- dir_tbl2[met_func != "Others"]

blue_funcs2 <- stat_4a2[
  met_func %in% dir_tbl2_no_others[direction == "NC_higher", met_func],
  .(met_func, diff = mean_H - mean_I)][order(diff)][, met_func]

red_funcs2 <- stat_4a2[
  met_func %in% dir_tbl2_no_others[direction == "IBD_higher", met_func],
  .(met_func, diff = mean_I - mean_H)][order(diff)][, met_func]

cols2 <- c(
  setNames(make_pal(length(blue_funcs2), c("#dbe9f6", "#08306b")), blue_funcs2),
  setNames(make_pal(length(red_funcs2),  c("#fde0dd", "#67000d")), red_funcs2),
  Others = "grey80"
)
func_levels2 <- c(rev(red_funcs2), blue_funcs2, "Others")

plot_4a_disease[, met_func_plot := factor(met_func_plot, levels = func_levels2)]
plot_4a_disease[, disease       := factor(disease, levels = c("NC", "IBD"))]

#Figure 5B
p_b <- ggplot(plot_4a_disease, aes(x = disease, y = mean_abund, fill = met_func_plot)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols2, name = "Metabolic functions") +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(colour = "black", linewidth = 0.6))) +
  labs(y = "Relative abundance", x = NULL) +
  theme_classic() +
  theme(
    legend.key.size = unit(0.45, "cm"),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 10),
    legend.key      = element_rect(fill = NA, colour = NA)
  )

p_b#View



#Figure 5C
#CD vs NC / UC vs NC /
#Differential pathway plot (left: mean relative abundance %; mid / right: CD vs NC / UC vs NC diff)
#p_c / Variable: p_c



#blocked WilcoxonCD vs NC; UC vs NC /
#Run blocked Wilcoxon per pathway (CD vs NC; UC vs NC)
analyze_pathways_blocked <- function(pw_dt, group_col, control, case, block_col = "project") {
  dt <- copy(pw_dt)[get(group_col) %in% c(control, case)]
  dt[, group2 := factor(get(group_col), levels = c(control, case))]
  dt[, block  := factor(get(block_col))]
  dt <- dt[!is.na(block)]
  valid_blocks <- dt[, uniqueN(group2), by = block][V1 == 2, block]
  dt <- dt[block %in% valid_blocks]

  setkey(dt, pathway)
  pws <- unique(dt$pathway)

  res <- rbindlist(lapply(pws, function(pw) {
    x <- dt[J(pw)]
    if (length(unique(x$block)) < 2)      return(NULL)
    if (any(table(x$group2) < 2))         return(NULL)
    if (sd(x$abund, na.rm = TRUE) == 0)   return(NULL)
    x <- x[!is.na(abund)]
    if (nrow(x) < 4)                       return(NULL)

    m_ctrl <- mean(x[group2 == control, abund])
    m_case <- mean(x[group2 == case,    abund])

    pv <- tryCatch(
      as.numeric(coin::pvalue(wilcox_test(abund ~ group2 | block, data = x, distribution = "approximate"))),
      error = function(e) NA_real_
    )
    data.table(pathway = pw, case = case, control = control,
               mean_case = m_case, mean_control = m_ctrl,
               diff = (m_case - m_ctrl), pval = pv)
  }), fill = TRUE)

  res[, fdr := p.adjust(pval, method = "BH")]
  res
}

pw_subtype <- pw_long[disease_subtype %in% c("NC", "CD", "UC")]
res_cd <- analyze_pathways_blocked(pw_subtype, "disease_subtype", control = "NC", case = "CD")
res_uc <- analyze_pathways_blocked(pw_subtype, "disease_subtype", control = "NC", case = "UC")
res_all <- rbind(res_cd, res_uc, fill = TRUE)

#CD+UC / Filter using case-side (CD+UC) stats
case_stats <- pw_subtype[disease_subtype %in% c("CD", "UC"), .(
  mean_case_all = mean(abund),
  prev_case_all = mean(abund > 0)
), by = pathway]


#At least one significant + mean_case >0.5% + prevalence >=5%
sig_any <- res_all[, .(sig_any = any(fdr < 0.05)), by = pathway]
keep_pw <- merge(case_stats, sig_any, by = "pathway")
keep_pw <- keep_pw[sig_any == TRUE & mean_case_all > 0.005 & prev_case_all >= 0.05, pathway]

res_plot <- res_all[pathway %in% keep_pw]
bar_plot <- case_stats[pathway %in% keep_pw]

#Sample sizes
n_nc <- uniqueN(pw_subtype[disease_subtype == "NC", sample])
n_cd <- uniqueN(pw_subtype[disease_subtype == "CD", sample])
n_uc <- uniqueN(pw_subtype[disease_subtype == "UC", sample])

#ID / Pathway label: strip ID prefix before colon
res_plot[, pathway_label := str_replace(pathway, "^.*?:\\s*", "")]
bar_plot[, pathway_label := str_replace(pathway, "^.*?:\\s*", "")]

#met_func / / Reuse met_func classification (for color / order)
tmp_map <- unique(pw_subtype[, .(pathway, met_func)])
res_plot <- merge(res_plot, tmp_map, by = "pathway", all.x = TRUE)
bar_plot <- merge(bar_plot,  tmp_map, by = "pathway", all.x = TRUE)

#Sort: by function category, then case mean descending
ord <- bar_plot[order(met_func, -mean_case_all)]
lvl <- ord$pathway_label
res_plot[, pathway_label := factor(pathway_label, levels = rev(lvl))]
bar_plot[, pathway_label := factor(pathway_label, levels = rev(lvl))]

#facet / facet labels with sample counts
facet_lab <- c(
  "NC" = sprintf("NC\n(n=%d)", n_nc),
  "CD" = sprintf("CD\n(n=%d)", n_cd),
  "UC" = sprintf("UC\n(n=%d)", n_uc)
)
res_plot[, case_lab := factor(facet_lab[case], levels = facet_lab[c("NC", "CD", "UC")])]

#Metabolic function color map (fixed order, as in paper)
met_levels <- c(
  "Amino-Acid-Biosynthesis", "Aminoacyl-tRNAs-Charging", "Cell-Wall-Biosynthesis",
  "Chorismate-Biosynthesis", "CO2-Fixation", "Coenzyme-A-Biosynthesis",
  "Deoxyribonucleotide-Biosynthesis", "Fatty-acid-biosynthesis", "Fermentation-to-Alcohols",
  "GLYCOLYSIS-VARIANTS", "NAD-Metabolism", "Nucleic-Acid-Processing",
  "Other-Amino-Acid-Biosynthesis", "Pentose-Phosphate-Cycle", "Phospholipid-Biosynthesis",
  "Polysaccharide Degradation", "Polysaccharides-Biosynthesis",
  "Proteinogenic-Amino-Acid-Biosynthesis", "Purine-Degradation", "Purine-Nucleotide-Biosynthesis",
  "Pyrimidine-Nucleotide-Biosynthesis", "SECONDARY-METABOLITE-BIOSYNTHESIS",
  "Sugar-Biosynthesis", "Sugars-Degradation", "Vitamin-Biosynthesis"
)
met_cols <- c(
  "Amino-Acid-Biosynthesis"                = "#666666",
  "Aminoacyl-tRNAs-Charging"               = "#999999",
  "Cell-Wall-Biosynthesis"                 = "#377EB8",
  "Chorismate-Biosynthesis"                = "#8DA0CB",
  "CO2-Fixation"                           = "#B3B3B3",
  "Coenzyme-A-Biosynthesis"                = "#984EA3",
  "Deoxyribonucleotide-Biosynthesis"       = "#7570B3",
  "Fatty-acid-biosynthesis"                = "#E41A1C",
  "Fermentation-to-Alcohols"               = "#E7298A",
  "GLYCOLYSIS-VARIANTS"                    = "#F781BF",
  "NAD-Metabolism"                         = "#E78AC3",
  "Nucleic-Acid-Processing"                = "#A65628",
  "Other-Amino-Acid-Biosynthesis"          = "#A6761D",
  "Pentose-Phosphate-Cycle"                = "#E6AB02",
  "Phospholipid-Biosynthesis"              = "#E5C494",
  "Polysaccharide Degradation"             = "#D95F02",
  "Polysaccharides-Biosynthesis"           = "#FF7F00",
  "Proteinogenic-Amino-Acid-Biosynthesis"  = "#FC8D62",
  "Purine-Degradation"                     = "#FFD92F",
  "Purine-Nucleotide-Biosynthesis"         = "#FFFF33",
  "Pyrimidine-Nucleotide-Biosynthesis"     = "#A6D854",
  "SECONDARY-METABOLITE-BIOSYNTHESIS"      = "#66A61E",
  "Sugar-Biosynthesis"                     = "#4DAF4A",
  "Sugars-Degradation"                     = "#66C2A5",
  "Vitamin-Biosynthesis"                   = "#1B9E77"
)

#Figure 5C left: mean relative abundance %
p_left <- ggplot(
  bar_plot,
  aes(y = pathway_label, x = mean_case_all * 100,
      fill = factor(met_func, levels = met_levels))
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = met_cols, breaks = met_levels, drop = FALSE) +
  guides(fill = guide_legend(title = "Metabolic functions", ncol = 2, byrow = FALSE)) +
  labs(x = "Relative abundance, %", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y       = element_text(size = 7),
    legend.position   = "left",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width  = unit(0.28, "cm")
  )

#Only significant pathways get asterisks
res_plot_only_CD_UC <- res_plot[case %in% c("CD", "UC")]
sig_dt <- copy(res_plot_only_CD_UC[!is.na(fdr) & fdr < 0.05])
sig_dt[, x_star := ifelse(diff >= 0, diff * 100 + 0.15, diff * 100 - 0.15)]

#Figure 5C right: case-control difference %
p_right <- ggplot(
  res_plot_only_CD_UC,
  aes(y = pathway_label, x = diff * 100,
      fill = factor(met_func, levels = met_levels))
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.75) +
  geom_text(
    data    = sig_dt,
    aes(x = diff * 100, label = "*", hjust = ifelse(diff >= 0, -0.15, 1.15)),
    vjust   = 0.5, size = 3, colour = "black"
  ) +
  facet_grid(. ~ case_lab, scales = "free_x") +
  scale_fill_manual(values = met_cols, breaks = met_levels, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  coord_cartesian(clip = "off") +
  labs(x = "Difference of relative abundance (case - control), %", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
    strip.text       = element_text(size = 8),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    panel.spacing.x  = unit(0.25, "cm"),
    axis.text.x      = element_text(size = 7, angle = 45, hjust = 1),
    plot.margin      = margin(5.5, 14, 5.5, 0)
  )

#Figure 5C combined
p_c <- p_left + p_right + plot_layout(widths = c(1.25, 2.25))
p_c#View



#Figure 5D
#IBD NC /
#Species interaction network (IBD-discriminatory taxa; NC-only correlations)
#p_d / Variable: p_d



#LMM out_tbl dat / Load LMM workspace (contains out_tbl and dat)
load(find_file("lmm-IBD.RData"))

#Parameters ----------
FDR_CUT  <- 0.05#BH-FDR / Edge significance threshold
COR_CUT  <- 0.20#Absolute correlation threshold
PSEUDO   <- 1e-6#log / log pseudo count
USE_GROUP <- c("Health-associated", "IBD-associated")#Draw only these two node types

#LMM Health / IBD / Classify Health / IBD taxa from LMM ----------
da0 <- out_tbl %>%
  mutate(Comparison = str_squish(Comparison)) %>%
  tidyr::separate(Comparison, into = c("grp1", "grp2"), sep = " - ", remove = FALSE)

#disease - NC / Standardize direction to disease - NC
da_nc <- da0 %>%
  filter(
    (grp1 %in% c("UC", "CD") & grp2 == "NC") |
    (grp1 == "NC" & grp2 %in% c("UC", "CD"))
  ) %>%
  mutate(
    disease = ifelse(grp1 == "NC", grp2, grp1),
    est_disease_vs_NC = ifelse(grp1 == "NC", -Estimate, Estimate),
    fdr = AdjP_BH
  ) %>%
  dplyr::select(Taxon, disease, est_disease_vs_NC, fdr)

#Classify: pos-only / neg-only / mixed / not significant
tax_groups <- da_nc %>%
  group_by(Taxon) %>%
  summarise(
    any_pos = any(fdr < FDR_CUT & est_disease_vs_NC > 0, na.rm = TRUE),
    any_neg = any(fdr < FDR_CUT & est_disease_vs_NC < 0, na.rm = TRUE),
    min_fdr = suppressWarnings(min(fdr, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(group = case_when(
    any_pos & !any_neg ~ "IBD-associated",
    any_neg & !any_pos ~ "Health-associated",
    any_pos & any_neg  ~ "Mixed",
    TRUE               ~ "NotSig"
  )) %>%
  filter(group %in% USE_GROUP)

#NC / Build correlation matrix (NC samples only) ----------
dat_fixed <- dat
#make.names / Normalize column names for taxa (cols 21+)
colnames(dat_fixed)[21:ncol(dat_fixed)] <- make.names(colnames(dat_fixed)[21:ncol(dat_fixed)], unique = TRUE)

taxa_in_dat <- intersect(tax_groups$Taxon, colnames(dat_fixed))

dat_fixed_NC <- dat_fixed %>% dplyr::filter(disease_subtype == "NC")
X <- dat_fixed_NC[, taxa_in_dat, drop = FALSE]

#log /
#Force numeric; divide by 100 if percentage; log-transform
X <- X %>% mutate(across(everything(), ~as.numeric(as.character(.x)))) %>% as.matrix()
if (max(X, na.rm = TRUE) > 1.5) X <- X / 100
X_log <- log(X + PSEUDO)

#Spearman + Pearson BH-FDR /
#Compute Spearman + Pearson correlation, dual BH-FDR filtering ----------
rc_s <- Hmisc::rcorr(X_log, type = "spearman")
rc_p <- Hmisc::rcorr(X_log, type = "pearson")

get_edges <- function(R, P, suffix) {
  ut <- upper.tri(R, diag = FALSE)
  data.frame(
    from = colnames(R)[row(R)[ut]],
    to   = colnames(R)[col(R)[ut]],
    cor  = as.numeric(R[ut]),
    p    = as.numeric(P[ut]),
    stringsAsFactors = FALSE
  ) %>% rename_with(~paste0(.x, "_", suffix), c("cor", "p"))
}

e_s <- get_edges(rc_s$r, rc_s$P, "s")
e_p <- get_edges(rc_p$r, rc_p$P, "p")

edges <- dplyr::left_join(e_s, e_p, by = c("from", "to")) %>%
  mutate(
    FDR_s = p.adjust(p_s, method = "BH"),
    FDR_p = p.adjust(p_p, method = "BH")
  )

#Filter: both significant + same sign + both above threshold
edges_f <- edges %>%
  filter(FDR_s < FDR_CUT, FDR_p < FDR_CUT) %>%
  filter(abs(cor_s) >= COR_CUT, abs(cor_p) >= COR_CUT) %>%
  filter(sign(cor_s) == sign(cor_p)) %>%
  mutate(cor_mean = (cor_s + cor_p) / 2,
         sign_lab = ifelse(cor_mean > 0, "Positive", "Negative"))

if (nrow(edges_f) == 0) {
  stop("没有任何边通过筛选：请放宽 FDR_CUT 或 COR_CUT。/ No edges passed filter: relax FDR_CUT or COR_CUT.")
}

#Node info: parse genus and species names ----------
parse_genus <- function(x) {
  x <- as.character(x)
  g <- ifelse(grepl("g__", x), sub(".*g__([^|;\\.]+).*", "\\1", x), x)
  g <- gsub("^g__", "", g)
  g <- gsub("\\.+", " ", g)#make.names / Restore make.names dots to spaces
  g
}
parse_species <- function(x) {
  x <- as.character(x)
  s <- ifelse(grepl("s__", x), sub(".*s__([^\\.]+).*", "\\1", x), x)
  s <- gsub("_", " ", s)
  s
}

vertices <- tax_groups %>%
  mutate(
    name    = Taxon,
    genus   = parse_genus(Taxon),
    species = parse_species(Taxon)
  ) %>%
  dplyr::select(name, group, genus, species)

#Build graph and remove isolated nodes ----------
g <- igraph::graph_from_data_frame(
  d        = edges_f %>% dplyr::select(from, to, sign_lab, cor_mean, FDR_s, FDR_p),
  directed = FALSE,
  vertices = vertices
)
g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
if (igraph::vcount(g) < 3) stop("去孤立点后节点太少。/ Too few nodes after removing isolates.")

#Health above / IBD below / Structured layout ----------
set.seed(1)
edge_df <- as.data.frame(igraph::as_data_frame(g, what = "edges"))
vert_df <- as.data.frame(igraph::as_data_frame(g, what = "vertices"))

edge_df <- edge_df %>%
  left_join(vert_df %>% dplyr::select(name, group), by = c("from" = "name")) %>%
  rename(group_from = group) %>%
  left_join(vert_df %>% dplyr::select(name, group), by = c("to"   = "name")) %>%
  rename(group_to   = group)

vdf <- data.frame(
  name    = V(g)$name,
  group   = V(g)$group,
  genus   = V(g)$genus,
  species = V(g)$species,
  x = NA_real_, y = NA_real_,
  stringsAsFactors = FALSE
)

#Use within-group edges for subgraph layout
edges_within <- edge_df %>% dplyr::filter(group_from == group_to)
g_within <- igraph::graph_from_data_frame(d = edges_within, directed = FALSE, vertices = vert_df)

y_shift <- 3
for (grp in unique(vdf$group)) {
  sub_nodes <- vdf$name[vdf$group == grp]
  subg <- igraph::induced_subgraph(g_within, vids = sub_nodes)
  L <- if (igraph::ecount(subg) > 0) igraph::layout_with_fr(subg) else igraph::layout_with_kk(subg)
  L   <- scale(L)
  idx <- match(V(subg)$name, vdf$name)
  vdf$x[idx] <- L[, 1]
  vdf$y[idx] <- L[, 2] + ifelse(grp == "Health-associated", y_shift, -y_shift)
}
vdf <- vdf[match(V(g)$name, vdf$name), ]

#Node degree + genus color map ----------
V(g)$degree <- igraph::degree(g)
vdf$degree  <- V(g)$degree

genus_levels <- sort(unique(vdf$genus))
genus_cols   <- setNames(scales::hue_pal()(length(genus_levels)), genus_levels)

#Figure 5D
p_d <- ggraph(g, layout = "manual", x = vdf$x, y = vdf$y) +
#Edges: colored by correlation direction
  geom_edge_link(aes(edge_colour = sign_lab, edge_width = abs(cor_mean))) +
  scale_edge_width(range = c(0.3, 1.5)) +
  scale_edge_colour_manual(
    values = c("Positive" = "black", "Negative" = "grey70"),
    name   = "Correlation"
  ) +
#Nodes
  geom_node_point(aes(fill = genus, size = degree), shape = 21, colour = "black", stroke = 0.4) +
  scale_size(range = c(5, 12)) +
#degree>=1 / Labels (species with degree >=1)
  ggrepel::geom_text_repel(
    data          = vdf %>% dplyr::filter(degree >= 1),
    aes(x = x, y = y, label = species),
    size          = 3, max.overlaps = 30, segment.color = "grey50"
  ) +
#Divider line and group labels
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = min(vdf$x, na.rm = TRUE), y = max(vdf$y, na.rm = TRUE) + 0.8,
           label = "Health-associated", hjust = 0, size = 5) +
  annotate("text", x = max(vdf$x, na.rm = TRUE), y = min(vdf$y, na.rm = TRUE) - 0.8,
           label = "IBD-associated",    hjust = 1, size = 5) +
  scale_fill_manual(values = genus_cols, name = "Genus") +
  theme_void(base_size = 12) +
  theme(legend.position = "right", legend.box = "vertical")

p_d#View




#Output: main figures + supplementary figures + supplementary tables



#Figure 5A main ------
ggsave(
  filename = file.path(output_dir, "Figure5A_pathway_subtype.pdf"),
  plot     = p_a,
  width    = 6, height = 4, units = "in", device = cairo_pdf
)

#Figure 5B main ------
ggsave(
  filename = file.path(output_dir, "Figure5B_pathway_disease.pdf"),
  plot     = p_b,
  width    = 5.5, height = 4, units = "in", device = cairo_pdf
)

#Figure 5C main ------
ggsave(
  filename = file.path(output_dir, "Figure5C_pathway_differential.pdf"),
  plot     = p_c,
  width    = 13, height = 9, units = "in", device = cairo_pdf
)

#Figure 5D main ------
ggsave(
  filename = file.path(output_dir, "Figure5D_IBD_interaction_network.pdf"),
  plot     = p_d,
  width    = 12, height = 9, units = "in", device = cairo_pdf
)

#Figure 10A (subtype) ------
#Figure 5A (subtype version)
ggsave(
  filename = file.path(output_dir, "Supplementary_Figure10A_subtype.pdf"),
  plot     = p_a,
  width    = 6, height = 4, units = "in", device = cairo_pdf
)

#Figure 10A (disease) ------
#Figure 5B (disease version)
ggsave(
  filename = file.path(output_dir, "Supplementary_Figure10A_disease.pdf"),
  plot     = p_b,
  width    = 5.5, height = 4, units = "in", device = cairo_pdf
)

#Figure 10C (subtype differential) ------
#Figure 5C
ggsave(
  filename = file.path(output_dir, "Supplementary_Figure10C_subtype.pdf"),
  plot     = p_c,
  width    = 13, height = 9, units = "in", device = cairo_pdf
)

#Supplementary Table 40: pathway long table ------
fwrite(pw_long, file.path(output_dir, "Supplementary_Table_40_pathway.csv"))

#NC vs IBD / Supplementary Table 41: NC vs IBD function stats ------
fwrite(stat_4a2, file.path(output_dir, "Supplementary_Table_41_pathway.csv"))

#CD / UC vs NC / Supplementary Table 42: pathway DA results (CD / UC vs NC) ------
fwrite(res_all, file.path(output_dir, "Supplementary_Table_42_pathway.csv"))

#NC / Supplementary Table 43: NC-only species interaction edges ------
write.csv(
  edges_f,
  file      = file.path(output_dir, "Supplementary_Table_43_NC_only_species_interaction_edges.csv"),
  row.names = FALSE
)


#Network statistics summary

N_edges  <- nrow(edges_f)
N_pos    <- sum(edges_f$sign_lab == "Positive")
N_neg    <- sum(edges_f$sign_lab == "Negative")
N_nodes  <- length(unique(c(edges_f$from, edges_f$to)))
prop_neg <- round(N_neg / N_edges, 3)

cat(sprintf(
  "Network: %d nodes, %d edges (%d positive, %d negative; %.1f%% negative)\n",
  N_nodes, N_edges, N_pos, N_neg, prop_neg * 100
))

#Hub taxa (by degree) / Hub
deg     <- igraph::degree(g)
top_hub <- head(sort(deg, decreasing = TRUE), 10)
print(top_hub)


#Figure 5 A-D /
#Figure 5 combined (A-D)
#A+B C D /
#Layout: A+B top row, C middle, D bottom

library(patchwork)

fig5_combined <- (
  wrap_elements(full = p_a) | wrap_elements(full = p_b)
) /
  wrap_elements(full = p_c) /
  wrap_elements(full = p_d) +
  plot_layout(heights = c(1.2, 2.0, 2.0)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 16, face = "bold"))
  )

ggsave(
  filename = file.path(output_dir, "Figure5_A-D_combined.pdf"),
  plot     = fig5_combined,
  width    = 20,
  height   = 30,
  units    = "in",
  device   = cairo_pdf
)
cat("Figure5_A-D_combined.pdf saved.\n")
