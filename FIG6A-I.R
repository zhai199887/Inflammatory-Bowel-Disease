required_pkgs <- c(
  "dplyr", "forcats", "ggplot2", "ggrepel", "pROC",
  "randomForest", "RColorBrewer", "patchwork",
  "pdftools", "png", "grid", "ggplotify"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

get_script_dir <- function() {
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(cmd_file) > 0) {
    cmd_dir <- dirname(cmd_file[1])
    if (dir.exists(cmd_dir)) {
      return(normalizePath(cmd_dir, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
output_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
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

save_panel <- function(plot_obj, filename, width, height) {
  outfile <- file.path(output_dir, filename)
  if (inherits(plot_obj, "recordedplot")) {
    grDevices::pdf(outfile, width = width, height = height, onefile = FALSE)
    replayPlot(plot_obj)
    grDevices::dev.off()
  } else {
    ggplot2::ggsave(outfile, plot = plot_obj, device = cairo_pdf, width = width, height = height, units = "in")
  }
}

save_base_panel <- function(draw_fun, filename, width, height) {
  outfile <- file.path(output_dir, filename)
  grDevices::pdf(outfile, width = width, height = height, onefile = FALSE)
  on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
  draw_fun()
  grDevices::dev.off()
}

base_panel_to_plot <- function(filename, draw_fun, width, height) {
  save_base_panel(draw_fun, filename, width, height)
  pdf_to_plot(file.path(output_dir, filename))
}

make_placeholder_plot <- function(title_text) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title_text, size = 6) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

draw_roc_manual <- function(roc_obj, col, lwd = 2, add = FALSE,
                            xlab = "1-Specificity", ylab = "Sensitivity",
                            xlim = c(0, 1), ylim = c(0, 1), main = NULL) {
  x <- 1 - roc_obj$specificities
  y <- roc_obj$sensitivities
  ord <- order(x, y)
  if (isTRUE(add)) {
    graphics::lines(x[ord], y[ord], col = col, lwd = lwd)
  } else {
    graphics::plot(
      x[ord], y[ord],
      type = "l",
      col = col,
      lwd = lwd,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      main = main
    )
  }
}

pdf_to_plot <- function(path, dpi = 200) {
  png_file <- pdftools::pdf_convert(path, format = "png", dpi = dpi, pages = 1)
  on.exit(unlink(png_file[1]), add = TRUE)
  bitmap <- png::readPNG(png_file[1])
  ggplotify::as.ggplot(grid::rasterGrob(bitmap, interpolate = FALSE))
}

copy_into_output <- function(src_name) {
  src <- find_file(src_name)
  dst <- file.path(output_dir, basename(src_name))
  if (normalizePath(src, winslash = "/", mustWork = TRUE) != normalizePath(dst, winslash = "/", mustWork = FALSE)) {
    file.copy(src, dst, overwrite = TRUE)
  }
}

load(find_file("Figure6_tenfold_workspace.RData"))
load(find_file("MIRS.RData"))
load(find_file("lasso.RData"))

stage_colors <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

#figA
p_a <- rf_importance_top19 %>%
  dplyr::mutate(protein = forcats::fct_reorder(protein, MeanDecreaseAccuracy)) %>%
  ggplot2::ggplot() +
  ggplot2::geom_col(
    ggplot2::aes(x = protein, y = MeanDecreaseAccuracy, fill = MeanDecreaseAccuracy),
    width = 0.8,
    color = NA
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  ggplot2::scale_fill_gradient(low = "#1B9E77", high = "#7570B3") +
  ggplot2::labs(x = NULL, y = "Feature Importance") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "none",
    axis.text.y = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(size = 11)
  )

#figB
roc_b_110 <- roc_allmicro
roc_b_21 <- roc_micro21
roc_b_19 <- roc_core
draw_fig_b <- function() {
  draw_roc_manual(roc_b_110, col = "#fd9200", xlim = c(0, 1), ylim = c(0, 1))
  graphics::abline(0, 1, col = "grey80")
  draw_roc_manual(roc_b_110, col = "#fd9200", add = FALSE, xlim = c(0, 1), ylim = c(0, 1))
  draw_roc_manual(roc_b_21, col = "#d33524", add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
  draw_roc_manual(roc_b_19, col = "#1576ae", add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
  legend("bottomright",
         legend = c("AUC=0.87(110 genera)", "AUC=0.92 (21 genera)", "AUC=0.92(19 genera)"),
         col = c("#fd9200", "#d33524", "#1576ae"), lty = 1, cex = 0.82)
}
p_b <- make_placeholder_plot("Figure 6B will be generated below")

#figC
draw_fig_c <- function() {
  draw_roc_manual(roc_random_list[[1]], col = color_list[1], lwd = 2,
                  xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  graphics::abline(0, 1, col = "grey80")
  draw_roc_manual(roc_random_list[[1]], col = color_list[1], lwd = 2,
                  xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  for (i in 2:10) {
    draw_roc_manual(roc_random_list[[i]], add = TRUE, col = color_list[i], lwd = 2,
                    xlim = c(0, 1), ylim = c(0, 1))
  }
  draw_roc_manual(roc_core, add = TRUE, col = "#1576ae", lwd = 2.5,
                  xlim = c(0, 1), ylim = c(0, 1))
  legend_labels <- c(
    paste0("Random subset ", 1:10, " (AUC = ", sprintf("%.2f", auc_random), ")"),
    "Core 19 genera (AUC = 0.92)"
  )
  legend("bottomright", legend = legend_labels,
         col = c(color_list, "#1576ae"), lty = 1, lwd = 2, cex = 0.78, bty = "n")
}
p_c <- make_placeholder_plot("Figure 6C will be generated below")

#figD
predict_test_prob <- predict(rf_top, newdata = test_top, type = "prob")
roc_micro19 <- pROC::roc(
  response = test_top$disease,
  predictor = as.numeric(predict_test_prob[, "IBD"]),
  levels = c("healthy", "IBD"),
  direction = "<"
)

rf_age_country <- randomForest::randomForest(
  disease ~ age_category + country,
  data = dat,
  ntree = 1000,
  importance = TRUE,
  na.action = na.omit
)
pred_age_country <- predict(rf_age_country, newdata = test_data, type = "prob")
roc_age_country <- pROC::roc(
  response = test_top$disease,
  predictor = as.numeric(pred_age_country[, "IBD"]),
  levels = c("healthy", "IBD"),
  direction = "<"
)

train_full <- data.frame(
  disease = train_data$disease,
  age_category = train_data$age_category,
  country = train_data$country,
  train_data[, selectvar, drop = FALSE]
)
test_full <- data.frame(
  disease = test_data$disease,
  age_category = test_data$age_category,
  country = test_data$country,
  test_data[, selectvar, drop = FALSE]
)
rf_full <- randomForest::randomForest(
  disease ~ .,
  data = train_full,
  ntree = 1000,
  importance = TRUE,
  na.action = na.omit
)
pred_full <- predict(rf_full, newdata = test_full, type = "prob")
roc_full <- pROC::roc(
  response = test_full$disease,
  predictor = as.numeric(pred_full[, "IBD"]),
  levels = c("healthy", "IBD"),
  direction = "<"
)

draw_fig_d <- function() {
  graphics::plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "1-Specificity",
    ylab = "Sensitivity",
    type = "n"
  )
  draw_roc_manual(roc_full, col = "#1576ae", lwd = 2, add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
  draw_roc_manual(roc_micro19, add = TRUE, col = "#d33524", lwd = 2, xlim = c(0, 1), ylim = c(0, 1))
  draw_roc_manual(roc_age_country, add = TRUE, col = "#fd9200", lwd = 2, xlim = c(0, 1), ylim = c(0, 1))
  legend(
    "bottomright",
    legend = c(
      paste0("AUC=", sprintf("%.2f", as.numeric(pROC::auc(roc_age_country))), "(Age + Country)"),
      paste0("AUC=", sprintf("%.2f", as.numeric(pROC::auc(roc_micro19))), " (19 genera)"),
      paste0("AUC=", sprintf("%.2f", as.numeric(pROC::auc(roc_full))), "(Age + Country +\n19 genera)")
    ),
    col = c("#fd9200", "#d33524", "#1576ae"),
    lty = 1,
    cex = 0.9,
    bty = "n"
  )
}
p_d <- make_placeholder_plot("Figure 6D will be generated below")

#figE
cor_stage <- out %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    r_prev = stats::cor(prob, prevalence, method = "pearson"),
    r_inci = stats::cor(prob, incidence, method = "pearson"),
    .groups = "drop"
  )
cor_long <- tidyr::pivot_longer(cor_stage, cols = c(r_prev, r_inci),
                                names_to = "metric", values_to = "r")
cor_long$metric <- ifelse(cor_long$metric == "r_prev", "prevalence", "incidence")

p_e <- ggplot2::ggplot(out_long, ggplot2::aes(x = prob, y = value, color = factor(stage))) +
  ggplot2::geom_point(size = 3, alpha = 0.8) +
  ggplot2::geom_smooth(ggplot2::aes(color = factor(stage)), method = "lm", se = FALSE, linetype = "dashed") +
  ggplot2::facet_wrap(~ metric, scales = "free_y") +
  ggplot2::scale_color_manual(values = stage_colors, name = "Stage") +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  ggplot2::labs(x = "MIRS score", y = "Epidemiological metric")

#figF
draw_fig_f <- function() {
  roc_first <- pROC::roc(y_binary, df[[top_feats[1]]], quiet = TRUE)
  draw_roc_manual(roc_first, col = cols[1], lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  graphics::abline(0, 1, col = "grey80")
  draw_roc_manual(roc_first, col = cols[1], lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  auc_labels <- paste0(top_feats[1], " (AUC=", sprintf("%.2f", auc(roc_first)), ")")
  if (topN >= 2) {
    for (i in 2:topN) {
      rci <- pROC::roc(y_binary, df[[top_feats[i]]], quiet = TRUE)
      draw_roc_manual(rci, add = TRUE, col = cols[i], lwd = 2, xlim = c(0, 1), ylim = c(0, 1))
      auc_labels <- c(auc_labels, paste0(top_feats[i], " (AUC=", sprintf("%.2f", auc(rci)), ")"))
    }
  }
  draw_roc_manual(roc_combined, add = TRUE, col = "#1576ae", lwd = 2.5, xlim = c(0, 1), ylim = c(0, 1))
  legend("bottomright",
         legend = c(auc_labels, paste0("Combined model (AUC=", sprintf("%.2f", auc(roc_combined)), ")")),
         col = c(cols[1:topN], "#1576ae"), lty = 1, lwd = 2, cex = 0.65, bty = "n")
}
p_f <- make_placeholder_plot("Figure 6F will be generated below")

#figG
draw_fig_g <- function() {
  roc_first1 <- pROC::roc(y_binary1, df[[top_feats1[1]]], quiet = TRUE)
  draw_roc_manual(roc_first1, col = cols1[1], lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  graphics::abline(0, 1, col = "grey80")
  draw_roc_manual(roc_first1, col = cols1[1], lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
                  main = "")
  auc_labels1 <- paste0(top_feats1[1], " (AUC=", sprintf("%.2f", auc(roc_first1)), ")")
  if (topN1 >= 2) {
    for (i in 2:topN1) {
      rci1 <- pROC::roc(y_binary1, df[[top_feats1[i]]], quiet = TRUE)
      draw_roc_manual(rci1, add = TRUE, col = cols1[i], lwd = 2, xlim = c(0, 1), ylim = c(0, 1))
      auc_labels1 <- c(auc_labels1, paste0(top_feats1[i], " (AUC=", sprintf("%.2f", auc(rci1)), ")"))
    }
  }
  draw_roc_manual(roc_combined1, add = TRUE, col = "#1576ae", lwd = 2.5, xlim = c(0, 1), ylim = c(0, 1))
  legend("bottomright",
         legend = c(auc_labels1, paste0("Combined model (AUC=", sprintf("%.2f", auc(roc_combined1)), ")")),
         col = c(cols1[1:topN1], "#1576ae"), lty = 1, lwd = 2, cex = 0.65, bty = "n")
}
p_g <- make_placeholder_plot("Figure 6G will be generated below")

#figH
p_h <- ggplot2::ggplot(df, ggplot2::aes(x = log_incidence, y = predicted_incidence)) +
  ggplot2::geom_point(size = 3, color = "#1576ae", alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", color = "#d33524", linetype = "dashed", se = TRUE) +
  ggplot2::labs(x = "Observed log(IBD incidence)", y = "Predicted log(IBD incidence)",
                title = "Predicted vs. Observed IBD incidence") +
  ggplot2::annotate("text", x = min(df$log_incidence), y = max(df$predicted_incidence),
                    label = "r = 0.45, p < 0.001", hjust = 0, size = 5, color = "black") +
  ggplot2::theme_classic(base_size = 14)

#figI
p_i <- ggplot2::ggplot(df, ggplot2::aes(x = log_prevalence, y = predicted_prevalence)) +
  ggplot2::geom_point(size = 3, color = "#fd9200", alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#d33524", linetype = "dashed") +
  ggplot2::labs(x = "Observed log(IBD prevalence)", y = "Predicted log(IBD prevalence)",
                title = "Predicted vs. Observed IBD prevalence") +
  ggplot2::annotate(
    "text",
    x = min(df$log_prevalence),
    y = max(df$predicted_prevalence),
    label = paste0("r = ", round(cor_test_pre$estimate, 2),
                   ", p = ", ifelse(cor_test_pre$p.value < 0.001, "<0.001", signif(cor_test_pre$p.value, 2))),
    hjust = 0,
    size = 5
  ) +
  ggplot2::theme_classic(base_size = 14)

#supplementary figure 20
p_s20 <- ggplot2::ggplot(
  data.frame(Genus = names(sorted_importance1), Coefficient = as.numeric(sorted_importance1)),
  ggplot2::aes(x = Coefficient, y = reorder(Genus, Coefficient), fill = Coefficient > 0)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"), guide = "none") +
  ggplot2::labs(x = "Coefficient Value", y = NULL, title = "Genera Impact Direction on log(IBD incidence)") +
  ggplot2::theme_bw(base_size = 12)

#supplementary figure 21
p_s21 <- ggplot2::ggplot(
  data.frame(Genus = names(sorted_importance), Coefficient = as.numeric(sorted_importance)),
  ggplot2::aes(x = Coefficient, y = reorder(Genus, Coefficient), fill = Coefficient > 0)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"), guide = "none") +
  ggplot2::labs(x = "Coefficient Value", y = NULL, title = "Genera Impact Direction on log(IBD prevalence)") +
  ggplot2::theme_bw(base_size = 12)

save_panel(p_a, "Figure6A_feature_importance.pdf", 3.0, 2.6)
p_b <- base_panel_to_plot("Figure6B_all_vs_21_vs_19_genera_ROC.pdf", draw_fig_b, 4.8, 4.8)
p_c <- base_panel_to_plot("Figure6C_core19_vs_random19_ROC.pdf", draw_fig_c, 4.8, 4.8)
p_d <- base_panel_to_plot("Figure6D_age_country_plus_19genera_ROC.pdf", draw_fig_d, 4.8, 4.8)
save_panel(p_e, "Figure6E_MIRS_vs_epidemiology.pdf", 8.4, 4.8)
p_f <- base_panel_to_plot("Figure6F_combined_vs_single_genera_prevalence_ROC.pdf", draw_fig_f, 5.2, 4.8)
p_g <- base_panel_to_plot("Figure6G_combined_vs_single_genera_incidence_ROC.pdf", draw_fig_g, 5.2, 4.8)
save_panel(p_h, "Figure6H_log_incidence_correlation.pdf", 4.6, 3.8)
save_panel(p_i, "Figure6I_log_prevalence_correlation.pdf", 4.6, 3.8)
save_panel(p_s20, "Supplementary_Figure20_feature_importance_incidence.pdf", 7.0, 5.0)
save_panel(p_s21, "Supplementary_Figure21_feature_importance_prevalence.pdf", 7.0, 5.0)

copy_into_output("Supplementary_Table_29_selected_genera_LASSO.csv")
copy_into_output("Supplementary_Table_30_LASSO.csv")
copy_into_output("Supplementary_Table_31_selected_genera_LASSO.csv")

panel_files <- c(
  "Figure6A_feature_importance.pdf",
  "Figure6B_all_vs_21_vs_19_genera_ROC.pdf",
  "Figure6C_core19_vs_random19_ROC.pdf",
  "Figure6D_age_country_plus_19genera_ROC.pdf",
  "Figure6E_MIRS_vs_epidemiology.pdf",
  "Figure6F_combined_vs_single_genera_prevalence_ROC.pdf",
  "Figure6G_combined_vs_single_genera_incidence_ROC.pdf",
  "Figure6H_log_incidence_correlation.pdf",
  "Figure6I_log_prevalence_correlation.pdf"
)
panel_plots <- lapply(panel_files, function(x) pdf_to_plot(file.path(output_dir, x)))
names(panel_plots) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

figure6_design <- "
AAABBBCCC
AAABBBCCC
DDDEEEEEE
DDDEEEEEE
DDDEEEEEE
FFFFGGGGG
FFFFGGGGG
FFFFGGGGG
HHHHIIIII
HHHHIIIII
"

p_figure6 <- patchwork::wrap_plots(
  A = panel_plots$A,
  B = panel_plots$B,
  C = panel_plots$C,
  D = panel_plots$D,
  E = panel_plots$E,
  F = panel_plots$F,
  G = panel_plots$G,
  H = panel_plots$H,
  I = panel_plots$I,
  design = figure6_design
)

ggplot2::ggsave(file.path(output_dir, "Figure6_A-I_composite.pdf"), plot = p_figure6,
                device = grDevices::pdf, width = 14.0, height = 19.0, units = "in")
