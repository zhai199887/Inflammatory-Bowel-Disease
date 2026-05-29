
#Figure 2 A–D: Global gut microbiome phylum composition
#Figures 12–14 and Supplementary Tables 9–12



#Packages & helper functions

required_packages <- c(
  "dplyr", "data.table", "tidyr", "tibble", "ggplot2",
  "patchwork", "cowplot", "scales", "colorspace", "readr",
  "reshape2", "boot", "purrr", "rstatix"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

#Locate the directory of this script (or working directory when run interactively)
get_script_dir <- function() {
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(cmd_file) > 0 && file.exists(cmd_file[1])) {
    return(dirname(normalizePath(cmd_file[1], winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
script_dir <- get_script_dir()

#Save a ggplot to PDF in the script directory
save_pdf <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(script_dir, filename),
    plot     = plot,
    device   = "pdf",
    width    = width,
    height   = height,
    units    = "in"
  )
}

#Convert a row vector to relative abundance (sums to 1)
rel <- function(x) {
  x     <- as.numeric(x)
  total <- sum(x, na.rm = TRUE)
  if (is.na(total) || total == 0) return(rep(0, length(x)))
  x / total
}

#Extract phylum name (second field of dot-separated taxon)
getPhylum <- function(taxon_name) {
  temp <- unlist(strsplit(taxon_name, split = "[.]"))[2]
  if (identical(temp, integer(0))) return(taxon_name)
  temp
}

#Convert count table to relative-abundance table (rows = samples)
#Note: apply(table, 1, rel) loses colnames because rel() calls as.numeric(),
#so we restore them explicitly after the transformation.
make_rel <- function(table) {
  mat <- t(apply(table, 1, rel))
  colnames(mat) <- colnames(table)
  data.frame(mat, check.names = FALSE)
}

#Collapse duplicate taxon columns by summing
combine_taxa <- function(dataset) {
  dataset$srr <- rownames(dataset)
  dt <- dataset %>%
    tidyr::pivot_longer(!srr, names_to = "taxon", values_to = "count") %>%
    data.table::as.data.table()
  final <- dt[, lapply(.SD, sum), by = list(srr, taxon)] %>%
    tidyr::pivot_wider(names_from = taxon, values_from = count, values_fill = 0) %>%
    data.frame()
  rownames(final) <- final$srr
  final$srr <- NULL
  final
}

#Rename new NCBI phylum names to classic names used in figures
rename_phyla_cols <- function(df) {
  colnames(df) <- gsub("Bacillota",               "Firmicutes",      colnames(df))
  colnames(df) <- gsub("Pseudomonadota",           "Proteobacteria",  colnames(df))
  colnames(df) <- gsub("Actinomycetota",           "Actinobacteria",  colnames(df))
  colnames(df) <- gsub("Bacteroidota",             "Bacteroidetes",   colnames(df))
  colnames(df) <- gsub("Thermodesulfobacteriota",  "Desulfobacterota",colnames(df))
  df
}

#Build a long-format myData object from a wide relative-abundance table,
#grouped by the `region` column already attached to final.rel.
#do_sample10: if TRUE, take a random 10 % of samples per region (seed 42)
#max_n: hard cap on samples per region (NULL = no cap).
#Set to avoid rendering 100k+ bars in a single PDF.
build_myData <- function(final_rel, regionDF, topTaxa,
                         do_sample10 = FALSE, max_n = NULL) {
  myData <- NULL
  for (i in seq_len(nrow(regionDF))) {
    currRegion <- regionDF$region[i]
    idx        <- which(final_rel$region == currRegion)
    if (length(idx) == 0) next

    df_sub <- final_rel[idx, , drop = FALSE]
    df_sub <- rename_phyla_cols(df_sub)

#Optional 10 % down-sampling
    if (do_sample10) {
      set.seed(42)
      n_take <- max(1, round(0.1 * nrow(df_sub)))
      df_sub <- df_sub[sample(nrow(df_sub), n_take), , drop = FALSE]
    }

#Hard cap: keep at most max_n samples (randomly selected, seed 42)
    if (!is.null(max_n) && nrow(df_sub) > max_n) {
      set.seed(42)
      df_sub <- df_sub[sample(nrow(df_sub), max_n), , drop = FALSE]
    }

    df_sub$sample <- rownames(df_sub)
    df_sub$region <- NULL

#Ensure every top taxon column is present
    for (tx in topTaxa) {
      if (!tx %in% colnames(df_sub)) df_sub[[tx]] <- 0
    }

#Sum all remaining (non-top) numeric columns into "other"
    other_cols <- setdiff(colnames(df_sub), c("sample", topTaxa))
    num_cols   <- other_cols[sapply(df_sub[, other_cols, drop = FALSE], is.numeric)]
    df_sub$other <- if (length(num_cols) == 0) 0 else
      rowSums(df_sub[, num_cols, drop = FALSE], na.rm = TRUE)

#Pivot to long format
    final.long <- tidyr::pivot_longer(df_sub, cols = -sample,
                                      names_to = "taxon", values_to = "rel")
    final.long$taxon[!(final.long$taxon %in% c(topTaxa, "other"))] <- "other"
    final.long$taxon <- factor(final.long$taxon, levels = c(topTaxa, "other"))

#Sort samples by cumulative Firmicutes → Proteobacteria → order
    ord_df  <- df_sub[, c("sample", topTaxa), drop = FALSE]
    ord_key <- with(ord_df,
      order(Firmicutes,
            Firmicutes + Proteobacteria,
            Firmicutes + Proteobacteria + Actinobacteria,
            Firmicutes + Proteobacteria + Actinobacteria + Bacteroidetes,
            Firmicutes + Proteobacteria + Actinobacteria + Bacteroidetes + Desulfobacterota))
    final.long$sample  <- factor(final.long$sample, levels = ord_df$sample[ord_key])
    final.long$region  <- currRegion

    myData <- if (is.null(myData)) final.long else rbind(myData, final.long)
  }

#Re-normalise each sample's phylum vector so rows sum to 1
  myData <- myData %>%
    group_by(sample) %>%
    mutate(rel = rel / sum(rel, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(rel))

#Re-sort samples within each region by ascending Firmicutes abundance
  myData <- myData %>%
    group_by(region, sample) %>%
    summarise(Firmicutes = sum(rel[taxon == "Firmicutes"], na.rm = TRUE),
              .groups = "drop") %>%
    arrange(region, Firmicutes) %>%
    mutate(sample = factor(sample, levels = unique(sample))) %>%
    dplyr::select(region, sample) %>%
    right_join(myData, by = c("region", "sample"))

  myData
}

#Build a stacked-bar plot + colour-strip inset from myData
build_stacked_bar <- function(myData, regionDF, colour_scale, region_order,
                               strip_height = 64) {
  myData$region <- factor(myData$region, levels = region_order)
  myData <- myData %>%
    arrange(region) %>%
    mutate(sample = factor(sample, levels = unique(sample))) %>%
    filter(!is.na(region))

  regionDF <- regionDF %>%
    mutate(region = factor(region, levels = region_order)) %>%
    arrange(region)

  phylum_scale <- c(
    Firmicutes       = "#BCD2EE",
    Proteobacteria   = "#832161",
    Actinobacteria   = "#06D6A0",
    Bacteroidetes    = "#E88873",
    Desulfobacterota = "#6153CC",
    other            = "gray"
  )

  bar <- ggplot(myData, aes(x = sample, y = rel, fill = taxon)) +
    geom_bar(stat = "identity", width = 1) +
    theme_bw() +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      axis.text        = element_text(size = 11),
      axis.title       = element_text(size = 11)
    ) +
    scale_fill_manual(values = phylum_scale,
                      aesthetics = c("colour", "fill")) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0, 0), limits = c(0, 1)) +
    labs(y = "Relative abundance", fill = "Phylum")

#Bottom colour strip: one tile per sample coloured by region
  breaks <- vapply(seq_len(nrow(regionDF)), function(i) {
    min(which(myData$region == regionDF$region[i]))
  }, numeric(1))

  colour_strip <- ggplot(myData) +
    geom_tile(aes(x = seq_len(nrow(myData)), y = 1, fill = region)) +
    scale_x_discrete(breaks = breaks, labels = regionDF$region, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colour_scale, aesthetics = c("colour", "fill")) +
    xlab("Region") +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks       = element_blank(),
      axis.title.y     = element_blank(),
      panel.spacing    = unit(c(0, 0, 0, 0), "pt"),
      plot.margin      = unit(c(0, 0, 0, 0), "pt")
    )

  bar + inset_element(colour_strip, 0, -0.2, 1, 0,
                      align_to = "panel", ignore_tag = TRUE)
}



#Load shared data

regions_df <- read.csv("D:/R代码/regions.csv", stringsAsFactors = FALSE)
metadata   <- readRDS("D:/R代码/metadata_for_diffAbundance.rds")
metadata   <- metadata %>% left_join(regions_df, by = c("country" = "country"))

meta <- data.table::fread("D:/R代码/Supplementary_Table_44_sample_metadata.tsv",
                           sep = "\t", header = TRUE) %>%
  mutate(sample_id = paste(project, srr, sep = "_")) %>%
  left_join(regions_df, by = "iso")

x_filtered <- readRDS("D:/R代码/filtered.rds")



#Figure 2A: Rarefaction / phylum-discovery curve by country


#Load pre-computed rarefaction results (one row per replicate × sample-size)
res <- readRDS("D:/R代码/rarefaction_phylum_by_iso_with_sampleid.rds")

res$millions     <- res$total_reads / 1e6
res$taxaPerRead  <- res$observed / res$total_reads
res$taxaPerMillion <- res$observed / res$millions

#Summarise to mean per (ISO × sample-size)
excluded_iso <- c("UM", "TW", "HK", "BW", "unknown")
output <- res %>%
  filter(!(iso %in% excluded_iso)) %>%
  group_by(currRegion = iso, currSampleSize = scount) %>%
  summarise(
    mean.reads    = mean(total_reads,    na.rm = TRUE),
    mean.taxa     = mean(observed,       na.rm = TRUE),
    taxaPerMillion = mean(taxaPerMillion, na.rm = TRUE),
    taxaPerRead   = mean(taxaPerRead,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(currSampleSize, mean.reads, mean.taxa,
                  taxaPerMillion, taxaPerRead), as.numeric))

#Label positions: keep only the row with the largest sample size per country
textPos <- output %>%
  group_by(currRegion) %>%
  slice_max(currSampleSize, n = 1) %>%
  ungroup() %>%
  mutate(
    mean.taxa = case_when(
      currRegion == "US" ~ mean.taxa + 100,
      currRegion == "CN" ~ mean.taxa - 50,
      TRUE               ~ mean.taxa
    )
  )

#Colour palette (64 countries)
iso_list <- c(
  "JO","SD","AZ","EE","MZ","GT","MX","CM","TZ","ML","RO","NG","KZ","EC","TH","RS","VE",
  "RU","CF","MG","MW","PE","SG","BR","PL","TR","LT","HR","PT","UG","AE","IN","IT","ZW",
  "AU","IE","DK","ES","BD","CO","KE","HU","NL","MY","NO","GR","AM","CH","CA","SE","AT",
  "IS","FR","GB","NZ","FI","DE","IL","BE","CN","KR","US","GH","JP"
)
expanded_colors <- qualitative_hcl(n = 64, palette = "Dark 3",
                                   h = c(0, 360), c = 60, l = 65)
names(expanded_colors) <- iso_list

#Main panel: mean unique phyla vs sample size
p_a <- ggplot(output, aes(y = mean.taxa, x = currSampleSize,
                           colour = currRegion, group = currRegion)) +
  geom_point() +
  geom_line() +
  geom_text(data = textPos,
            aes(x = currSampleSize + 1000, label = currRegion),
            hjust = 0, size = 3, show.legend = FALSE) +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = expanded_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Sample size", y = "Phylum richness") +
  theme_bw()

#Inset: phyla per million reads (log-scaled y)
p_a_inset <- ggplot(output, aes(y = taxaPerMillion, x = currSampleSize,
                                 colour = currRegion, group = currRegion)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_log10(labels = scales::comma_format()) +
  scale_colour_manual(values = expanded_colors, guide = "none") +
  labs(x = "Sample size", y = "Phyla per million reads") +
  theme(legend.position = "none", axis.title = element_text(size = 8))

p_a_combined <- p_a + inset_element(p_a_inset, 0.4, 0.05, 0.98, 0.65,
                                     align_to = "panel", ignore_tag = TRUE)
p_a_combined
save_pdf(p_a,          "Figure2A_scree.pdf",                           width = 7, height = 4)
save_pdf(p_a_combined, "Figure2A_microbiome_composition_phylum.pdf",   width = 7, height = 4)



#Figure 12: Phyla per read (not per sample) by country

p_supp12 <- ggplot(output, aes(y = taxaPerRead, x = currSampleSize,
                                colour = currRegion, group = currRegion)) +
  geom_point() +
  geom_line() +
  geom_text(data = textPos,
            aes(x = currSampleSize + max(output$currSampleSize) * 0.05,
                label = currRegion),
            hjust = 0, size = 3, show.legend = FALSE) +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::scientific_format(digits = 2),
                     limits = c(0, 6e-6)) +
  scale_colour_manual(values = expanded_colors, guide = "none") +
  labs(x = "Sample size", y = "Phyla per read") +
  theme_bw()
p_supp12
save_pdf(p_supp12, "Supplementary_Figure12_microbiome_composition_phylumPerRead.pdf",
         width = 7, height = 8)

#Export Supplementary Table 9
write.csv(output, file.path(script_dir, "Supplementary_Table_9_output.csv"),
          row.names = FALSE)



#Figure 14: Bootstrap CI – mean phyla vs mean reads per country

boot_mean_ci <- function(v, R = 1000) {
  v <- v[!is.na(v)]
  if (length(v) < 2 || length(unique(v)) == 1) {
    m <- mean(v, na.rm = TRUE)
    return(c(m, m))
  }
  b  <- boot(v, statistic = function(d, idx) mean(d[idx], na.rm = TRUE), R = R)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5],
                 error = function(e) c(mean(v, na.rm = TRUE), mean(v, na.rm = TRUE)))
  ci
}

ci_df <- res %>%
  filter(!(iso %in% excluded_iso)) %>%
  group_by(currRegion = iso, currSampleSize = scount) %>%
  summarise(
    mean_reads = mean(total_reads, na.rm = TRUE),
    mean_taxa  = mean(observed,    na.rm = TRUE),
    ci         = list(boot_mean_ci(observed, R = 1000)),
    .groups    = "drop"
  ) %>%
  mutate(ci_lower = map_dbl(ci, 1), ci_upper = map_dbl(ci, 2)) %>%
  dplyr::select(-ci)

p_supp14 <- ggplot(ci_df, aes(x = mean_reads / 1e6, y = mean_taxa,
                               colour = currRegion, group = currRegion)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, alpha = 0.3) +
  geom_text(
    data = ci_df %>% group_by(currRegion) %>% slice_max(currSampleSize, n = 1),
    aes(label = currRegion), hjust = -0.1, size = 2.8, show.legend = FALSE
  ) +
  scale_colour_manual(values = expanded_colors, guide = "none") +
  scale_x_continuous(name = "Mean reads (millions)", labels = scales::comma) +
  scale_y_continuous(name = "Mean unique phyla observed", labels = scales::comma) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 9),
        panel.grid.minor = element_blank())
p_supp14
save_pdf(p_supp14, "supplementary_Figure14__ci_df.pdf", width = 7.5, height = 5.2)



#Figure 2B: Phylum relative-abundance distribution by IBD stage
#freq-poly, log y-axis)


#Prepare phylum-level relative abundance table joined with metadata
x_phylum <- x_filtered
colnames(x_phylum) <- gsub("-", "", colnames(x_phylum))
colnames(x_phylum) <- gsub(" ", "", colnames(x_phylum))
colnames(x_phylum) <- gsub("^(\\w+\\.\\w+)\\.\\w+\\..+$", "\\1", colnames(x_phylum))

data_phylum <- combine_taxa(x_phylum)
colnames(data_phylum) <- gsub("^\\w+\\.(\\w+)$", "\\1", colnames(data_phylum))

#Merge duplicate columns
data_phylum <- as.data.frame(t(rowsum(t(data_phylum), group = colnames(data_phylum))))
rel_phylum  <- as.data.frame(t(apply(data_phylum, 1, function(x) {
  s <- sum(x, na.rm = TRUE); if (s > 0) x / s else rep(0, length(x))
})))
colnames(rel_phylum) <- colnames(data_phylum)

meta_fig2 <- read.delim("D:/R代码/Supplementary_Table_44_sample_metadata.tsv",
                         sep = "\t") %>%
  mutate(sample = paste(project, srr, sep = "_"))

toplot <- rel_phylum %>%
  mutate(sample = sub("^X", "", rownames(rel_phylum))) %>%
  tidyr::pivot_longer(!sample, names_to = "taxon", values_to = "rel") %>%
  left_join(meta_fig2, by = "sample") %>%
  left_join(regions_df, by = "iso") %>%
  filter(rel > 0)

#Shorten region labels
toplot$region <- toplot$region %>%
  gsub("South-Eastern", "S.E.", .) %>%
  gsub("Eastern",       "E.",   .) %>%
  gsub("Southern",      "S.",   .) %>%
  gsub("Northern",      "N.",   .) %>%
  gsub("Western",       "W.",   .) %>%
  gsub("Latin America", "Latin Amer.", .)

#Rename new phylum names to classic names
toplot$taxon <- toplot$taxon %>%
  gsub("Bacillota",               "Firmicutes",      .) %>%
  gsub("Actinomycetota",          "Actinobacteria",  .) %>%
  gsub("Pseudomonadota",          "Proteobacteria",  .) %>%
  gsub("Thermodesulfobacteriota", "Desulfobacterota",.) %>%
  gsub("Bacteroidota",            "Bacteroidetes",   .)

target_phyla <- c("Firmicutes", "Actinobacteria", "Proteobacteria",
                  "Desulfobacterota", "Bacteroidetes")
toplot <- toplot %>%
  mutate(taxon = ifelse(taxon %in% target_phyla, taxon, "other"))

phylum_palette <- c(
  Firmicutes       = "#BCD2EE",
  Actinobacteria   = "#06D6A0",
  Proteobacteria   = "#832161",
  Desulfobacterota = "#6153CC",
  Bacteroidetes    = "#E88873",
  other            = "gray70"
)

#Figure 2B: faceted by IBD stage
toplot_stage <- toplot %>%
  filter(!is.na(stage), !(iso %in% c("BW","TW","HK","UM","unknown"))) %>%
  mutate(stage = as.factor(stage))

p_b <- ggplot(toplot_stage, aes(x = rel, y = after_stat(count), colour = taxon)) +
  geom_freqpoly(linewidth = 1.0, binwidth = 0.02, na.rm = TRUE) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                limits = c(1, 5000)) +
  scale_x_continuous(
    labels = c("", scales::percent(0.5), "", scales::percent(1)),
    limits = c(0, 1), expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)
  ) +
  facet_wrap(~stage, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = phylum_palette) +
  labs(x = "Relative abundance", y = "Samples",
       title = "Distribution of major bacterial phyla across IBD stages",
       colour = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(
    plot.margin   = unit(c(20, 5, 5, 5), "pt"),
    axis.title    = element_text(size = 11),
    strip.text    = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )
p_b
save_pdf(p_b, "Figure2B_phylum_stage_plot.pdf", width = 10, height = 4)

#Export Supplementary Table 10 (toplot.sum)
write.csv(toplot, file.path(script_dir, "Supplementary_Table_10_toplot.sum.csv"),
          row.names = FALSE)



#Figure 13: Phylum distribution by country (64 ISO codes)

toplot_iso <- toplot %>%
  filter(!is.na(iso), !(iso %in% c("BW","TW","HK","UM","unknown"))) %>%
  mutate(iso = as.factor(iso))

p_supp13 <- ggplot(toplot_iso, aes(x = rel, y = after_stat(count), colour = taxon)) +
  geom_freqpoly(linewidth = 0.8, binwidth = 0.02, na.rm = TRUE) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                limits = c(1, 5000)) +
  scale_x_continuous(
    labels = c("", scales::percent(0.5), "", scales::percent(1)),
    limits = c(0, 1), expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)
  ) +
  facet_wrap(~iso, ncol = 8, scales = "free_y") +
  scale_colour_manual(values = phylum_palette) +
  labs(x = "Relative abundance", y = "Samples",
       title = "Distribution of major bacterial phyla across 64 countries",
       colour = "Phylum") +
  theme_minimal(base_size = 10) +
  theme(
    plot.margin  = unit(c(20, 5, 5, 5), "pt"),
    strip.text   = element_text(size = 7, face = "bold"),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 7),
    legend.position = "bottom"
  )
p_supp13
save_pdf(p_supp13, "Supplementary_Figure13_phylum_iso_plot.pdf",
         width = 15, height = 20)



#Figures 2C and 2D (stacked bar panels)


#Re-read regions (used for both C and D)
regions_df <- read.csv("D:/R代码/regions.csv", stringsAsFactors = FALSE)

#Reload metadata and filtered data for stacked-bar panels
metadata_bar <- readRDS("D:/R代码/metadata_for_diffAbundance.rds") %>%
  left_join(regions_df %>% dplyr::select(iso, country), by = "iso")

x_bar <- readRDS("D:/R代码/filtered.rds")

topTaxa <- c("Firmicutes", "Proteobacteria", "Actinobacteria",
             "Bacteroidetes", "Desulfobacterota")



#Figure 2C: Stacked bar – phylum by IBD stage (10 % down-sampled per stage)


#Assign region = stage
metadata_bar$region <- metadata_bar$stage
regionDF_stage <- regions_df %>%
  dplyr::select(stage) %>%
  distinct() %>%
  dplyr::rename(region = stage)

keep_stage  <- metadata_bar$region %in% regionDF_stage$region
meta_stage  <- metadata_bar[keep_stage, ]
x_stage     <- x_bar[keep_stage, ]

taxphylum_stage <- x_stage
colnames(taxphylum_stage) <- gsub("^(\\w+\\.\\w+)\\.\\w+\\..+$", "\\1",
                                   colnames(taxphylum_stage))
colnames(taxphylum_stage) <- sapply(colnames(taxphylum_stage), getPhylum)
final_rel_stage        <- combine_taxa(taxphylum_stage) %>% make_rel()
final_rel_stage$region <- meta_stage$region

stage_order  <- c("1", "2", "3")
stage_colors <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

myData_stage <- build_myData(final_rel_stage, regionDF_stage, topTaxa,
                              do_sample10 = TRUE)
myData_stage$region <- factor(myData_stage$region, levels = stage_order)

p_c <- build_stacked_bar(myData_stage, regionDF_stage,
                          stage_colors, stage_order)

save_pdf(p_c, "Figure2C_phylum_stage_stackedbar.pdf", width = 8.7, height = 3)



#Figure 2D: Stacked bar – phylum by country (10 % down-sampled per country, seed 42)


#Assign region = ISO code
metadata_bar$region <- metadata_bar$iso
regionDF_iso <- regions_df %>%
  dplyr::select(iso, country) %>%
  distinct() %>%
  dplyr::rename(region = iso)

keep_iso  <- metadata_bar$region %in% regionDF_iso$region
meta_iso  <- metadata_bar[keep_iso, ]
x_iso     <- x_bar[keep_iso, ]

taxphylum_iso <- x_iso
colnames(taxphylum_iso) <- gsub("^(\\w+\\.\\w+)\\.\\w+\\..+$", "\\1",
                                 colnames(taxphylum_iso))
colnames(taxphylum_iso) <- sapply(colnames(taxphylum_iso), getPhylum)
final_rel_iso        <- combine_taxa(taxphylum_iso) %>% make_rel()
final_rel_iso$region <- meta_iso$region

iso_order <- c(
  "US","CN","DK","JP","GB","FI","DE","CA","IN","NL","SE","NZ","IT","ES","FR","IL",
  "BD","BE","MW","RU","CM","MX","SG","MG","NO","TZ","AT","CH","AU","CF","EE","MZ",
  "KR","UG","ZW","PL","CO","RO","BR","EC","TR","TH","KE","GR","NG","GH","VE","PT",
  "ML","IE","GT","RS","AZ","SD","KZ","JO","LT","HU","HR","MY","AE","AM","IS","PE"
)

myData_iso <- build_myData(final_rel_iso, regionDF_iso, topTaxa,
                            do_sample10 = TRUE)#sampling, seed 42 (set inside build_myData)
myData_iso$region <- factor(myData_iso$region, levels = iso_order)

p_d <- build_stacked_bar(myData_iso, regionDF_iso,
                          expanded_colors, iso_order)

save_pdf(p_d, "Figure2D_phylum_country_stackedbar.pdf", width = 8.7, height = 3)



#Supplementary Tables 11 & 12: Pairwise Wilcoxon tests on phylum abundance


#Helper: run Kruskal + pairwise Wilcoxon, return tidy data frame
#max_per_group caps samples per region to keep pairwise tests tractable
run_phylum_wilcox <- function(myData, topTaxa, min_samples = 5,
                               max_per_group = 200) {
  df <- myData %>%
    filter(taxon %in% topTaxa) %>%
    dplyr::select(sample, region, taxon, rel) %>%
    distinct() %>%
    mutate(rel = as.numeric(rel)) %>%
    filter(!is.na(rel), !is.na(region), !is.na(taxon),
           !taxon %in% c("NA", "NA.", "Unclassified", "unknown"))

  valid_regions <- df %>%
    group_by(region) %>%
    summarise(n = n_distinct(sample), .groups = "drop") %>%
    filter(n >= min_samples) %>%
    pull(region)

  df_sub <- df %>% filter(region %in% valid_regions)

#Cap samples per region so pairwise tests finish in reasonable time
  if (!is.null(max_per_group) && max_per_group > 0) {
    set.seed(42)
    keep_samples <- df_sub %>%
      dplyr::select(sample, region) %>%
      distinct() %>%
      group_by(region) %>%
      slice_sample(n = max_per_group) %>%
      ungroup() %>%
      pull(sample)
    df_sub <- df_sub %>% filter(sample %in% keep_samples)
  }

#Kruskal test to identify taxa with significant regional variation
  sig_taxa <- df_sub %>%
    group_by(taxon) %>%
    kruskal_test(rel ~ region) %>%
    filter(p < 0.05) %>%
    pull(taxon)

#Pairwise Wilcoxon on significant taxa only
  df_sub %>%
    filter(taxon %in% sig_taxa) %>%
    group_by(taxon) %>%
    summarise(
      test = list({
        d           <- pick(everything())
        valid_grps  <- d %>%
          group_by(region) %>%
          summarise(n_valid = sum(rel > 0, na.rm = TRUE), .groups = "drop") %>%
          filter(n_valid >= 5)
        if (nrow(valid_grps) < 2) {
          tibble(group1=NA, group2=NA, n1=NA, n2=NA,
                 statistic=NA, p=NA, p.adj=NA, p.adj.signif=NA)
        } else {
          pairwise_wilcox_test(
            data      = d %>% filter(region %in% valid_grps$region),
            formula   = rel ~ region,
            p.adjust.method = "BH"
          )
        }
      }),
      .groups = "drop"
    ) %>%
    tidyr::unnest(test) %>%
    mutate(signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )) %>%
    arrange(taxon, p.adj)
}

#Supplementary Table 11: phylum differences across countries (ISO)
wilcox_iso <- run_phylum_wilcox(myData_iso, topTaxa)
write.csv(wilcox_iso,
          file.path(script_dir, "Supplementary_Table_11_wilcox_results by iso.csv"),
          row.names = FALSE)

#Supplementary Table 12: phylum differences across IBD stages
wilcox_stage <- run_phylum_wilcox(myData_stage, topTaxa)
write.csv(wilcox_stage,
          file.path(script_dir, "Supplementary_Table_12_wilcox_results by stage.csv"),
          row.names = FALSE)
