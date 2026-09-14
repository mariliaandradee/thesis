library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ragg)

# =========================================================
# 1. INPUT
# =========================================================

# pasta com copia de todas as figuras usadas na tese, nomeadas pelo numero
# de figura do capitulo de Resultados (facil de saber qual grafico e qual,
# sobretudo neste script onde as comparacoes envolvem tempo/dpi)
figuras_tese_dir <- "/Users/mariliaandrade/Desktop/Tese/resultados_final2/figuras_tese"
dir.create(figuras_tese_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- "/Users/mariliaandrade/Desktop/Tese/featurecounts_output.txt"
coldata_file <- "/Users/mariliaandrade/Desktop/Tese/updated_coldata.csv"

count_data <- read.table(counts_file, header = TRUE, row.names = 1)
count_data <- count_data[, 7:ncol(count_data)]
count_data <- as.matrix(count_data)
mode(count_data) <- "numeric"

coldata <- read.csv2(coldata_file, row.names = 1)

sample_columns <- gsub(".*\\.([0-9A-Z]+)_S[0-9]+_L[0-9]+.*", "\\1", colnames(count_data))
sample_columns <- sub("^X\\.", "", sample_columns)
sample_columns <- trimws(tolower(sample_columns))

colnames(count_data) <- sample_columns
rownames(coldata) <- trimws(tolower(rownames(coldata)))

# remover duplicadas no count_data
count_data <- count_data[, !duplicated(colnames(count_data))]

# alinhar coldata ao count_data
coldata <- coldata[colnames(count_data), , drop = FALSE]

# confirmar alinhamento
all(rownames(coldata) == colnames(count_data))

# padronizar fatores
coldata$Sampling <- factor(
  toupper(trimws(coldata$Sampling)),
  levels = c("4DPI", "14DPI")
)

coldata$Condition <- factor(
  toupper(trimws(coldata$Condition)),
  levels = c("CONTROL", "RCV", "RHDV")
)

coldata$Tissue <- factor(trimws(coldata$Tissue))

# =========================================================
# 2. OUTPUT PATHS
# =========================================================

base_path_tempo <- "/Users/mariliaandrade/Desktop/Tese/resultados_final/tempo"

dir.create(base_path_tempo, recursive = TRUE, showWarnings = FALSE)

dir_tempo_rcv_4dpi_comloc        <- file.path(base_path_tempo, "RCV_4dpi", "com_locs")
dir_tempo_rcv_4dpi_semloc        <- file.path(base_path_tempo, "RCV_4dpi", "sem_locs")

dir_tempo_rhdv_14dpi_comloc      <- file.path(base_path_tempo, "RHDV_14dpi", "com_locs")
dir_tempo_rhdv_14dpi_semloc      <- file.path(base_path_tempo, "RHDV_14dpi", "sem_locs")

dir_tempo_rcv_ctrl_4dpi_comloc   <- file.path(base_path_tempo, "RCV_CONTROL_4DPI", "com_locs")
dir_tempo_rcv_ctrl_4dpi_semloc   <- file.path(base_path_tempo, "RCV_CONTROL_4DPI", "sem_locs")

dir_tempo_rcv_ctrl_14dpi_comloc  <- file.path(base_path_tempo, "RCV_CONTROL_14DPI", "com_locs")
dir_tempo_rcv_ctrl_14dpi_semloc  <- file.path(base_path_tempo, "RCV_CONTROL_14DPI", "sem_locs")

dir_tempo_rhdv_ctrl_4dpi_comloc  <- file.path(base_path_tempo, "RHDV_CONTROL_4DPI", "com_locs")
dir_tempo_rhdv_ctrl_4dpi_semloc  <- file.path(base_path_tempo, "RHDV_CONTROL_4DPI", "sem_locs")

dir_tempo_rhdv_ctrl_14dpi_comloc <- file.path(base_path_tempo, "RHDV_CONTROL_14DPI", "com_locs")
dir_tempo_rhdv_ctrl_14dpi_semloc <- file.path(base_path_tempo, "RHDV_CONTROL_14DPI", "sem_locs")

dir_tempo_rhdv_rcv_comloc        <- file.path(base_path_tempo, "RHDV_RCV", "com_locs")
dir_tempo_rhdv_rcv_semloc        <- file.path(base_path_tempo, "RHDV_RCV", "sem_locs")

all_dirs <- c(
  dir_tempo_rcv_4dpi_comloc, dir_tempo_rcv_4dpi_semloc,
  dir_tempo_rhdv_14dpi_comloc, dir_tempo_rhdv_14dpi_semloc,
  dir_tempo_rcv_ctrl_4dpi_comloc, dir_tempo_rcv_ctrl_4dpi_semloc,
  dir_tempo_rcv_ctrl_14dpi_comloc, dir_tempo_rcv_ctrl_14dpi_semloc,
  dir_tempo_rhdv_ctrl_4dpi_comloc, dir_tempo_rhdv_ctrl_4dpi_semloc,
  dir_tempo_rhdv_ctrl_14dpi_comloc, dir_tempo_rhdv_ctrl_14dpi_semloc,
  dir_tempo_rhdv_rcv_comloc, dir_tempo_rhdv_rcv_semloc
)

for (d in all_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =========================================================
# 3. FUNCTIONS
# =========================================================

prepare_res_df <- function(res, remove_loc = FALSE) {
  df <- as.data.frame(res)
  df$Gene <- rownames(df)
  df <- df[, c("Gene", setdiff(colnames(df), "Gene")), drop = FALSE]
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), , drop = FALSE]
  
  if (remove_loc) {
    df <- df[!grepl("^LOC", df$Gene), , drop = FALSE]
  }
  
  df$negLogPadj <- -log10(pmax(df$padj, .Machine$double.xmin))
  df
}

save_result_tables <- function(res, out_dir, prefix, remove_loc = FALSE,
                               alpha = 0.05, lfc_cut = 1) {
  
  df <- prepare_res_df(res, remove_loc = remove_loc)
  df <- df[order(df$padj, -abs(df$log2FoldChange)), , drop = FALSE]
  
  write.csv(
    df,
    file.path(out_dir, paste0(prefix, "_all_genes.csv")),
    row.names = FALSE
  )
  
  sig_df <- subset(df, padj < alpha & abs(log2FoldChange) > lfc_cut)
  
  write.csv(
    sig_df,
    file.path(out_dir, paste0(prefix, "_sig.csv")),
    row.names = FALSE
  )
  
  up_df <- subset(sig_df, log2FoldChange > lfc_cut)
  down_df <- subset(sig_df, log2FoldChange < -lfc_cut)
  
  write.csv(
    up_df,
    file.path(out_dir, paste0(prefix, "_UP.csv")),
    row.names = FALSE
  )
  
  write.csv(
    down_df,
    file.path(out_dir, paste0(prefix, "_DOWN.csv")),
    row.names = FALSE
  )
  
  deg_counts <- data.frame(
    Comparison = prefix,
    Up = nrow(up_df),
    Down = nrow(down_df)
  )
  
  write.csv(
    deg_counts,
    file.path(out_dir, paste0(prefix, "_DEG_counts.csv")),
    row.names = FALSE
  )
  
  invisible(list(full = df, sig = sig_df, up = up_df, down = down_df, counts = deg_counts))
}

plot_volcano <- function(res, title, path_plot,
                         right_label, left_label,
                         remove_loc = TRUE,
                         lfc_min = -3, lfc_max = 3,
                         x_limit = 4,
                         y_limit = NULL,
                         padj_min = 0.05,
                         point_size = 2.5,
                         max.overlaps = 20,
                         highlight_genes = NULL,
                         top_n_label = 10,
                         fig_num = NULL,
                         fig_label = NULL) {

  df <- as.data.frame(res)
  df$Gene <- rownames(df)
  df <- df[!grepl("^LOC", df$Gene), ]
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]

  df$negLogPadj <- -log10(pmax(df$padj, .Machine$double.xmin))

  df$Significance <- "Not Sig"
  df$Significance[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$Significance[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"

  # genes citados no texto dos Resultados: sempre rotulados e destacados
  df$Highlight <- df$Gene %in% highlight_genes
  df_highlight <- df[df$Highlight, ]

  # rotular APENAS os genes destacados no texto; se nao houver nenhum,
  # rotular só os top_n_label mais significativos - nunca a nuvem toda
  if (nrow(df_highlight) > 0) {
    ord <- match(highlight_genes, df_highlight$Gene)
    ord <- ord[!is.na(ord)]
    df_label <- df_highlight[ord, ]
  } else {
    sig_df <- df[df$Significance != "Not Sig", ]
    sig_df <- sig_df[order(sig_df$padj), ]
    df_label <- head(sig_df, top_n_label)
  }

  n_up   <- sum(df$Significance == "Up")
  n_down <- sum(df$Significance == "Down")

  p <- ggplot(df, aes(x = log2FoldChange, y = negLogPadj)) +
    geom_point(aes(color = Significance), alpha = 0.55, size = point_size) +
    geom_point(
      data = df_highlight,
      shape = 21, fill = NA, color = "black",
      stroke = 1, size = point_size + 1.5
    ) +
    geom_text_repel(
      data = df_label,
      aes(label = Gene),
      size = 3.4,
      fontface = "bold",
      color = "black",
      bg.color = "white",
      bg.r = 0.15,
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.3,
      force = 3,
      force_pull = 0.6,
      max.iter = 20000,
      max.time = 2,
      seed = 42,
      min.segment.length = 0,
      segment.size = 0.3,
      segment.color = "grey30",
      segment.alpha = 0.7
    ) +
    scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F5FB2", "Not Sig" = "grey82")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(title = title,
         subtitle = paste0("← Higher in ", left_label,
                            "                              Higher in ", right_label, " →"),
         caption = paste0("Total DEGs: ", n_up + n_down,
                           " (Up: ", n_up, ", Down: ", n_down,
                           "). Labelled genes are those discussed in the text; ",
                           "the full list is provided in the corresponding sig.csv / Supplementary Table."),
         x = "Log2 Fold Change",
         y = "-Log10(adjusted p-value)",
         color = "Expression") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 11),
      plot.caption = element_text(hjust = 0.5, size = 8.5, color = "grey30",
                                   margin = margin(t = 10)),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(path_plot, plot = p, width = 11, height = 8.5, dpi = 320, bg = "white",
         device = ragg::agg_png)

  if (!is.null(fig_num)) {
    ggsave(
      file.path(figuras_tese_dir, paste0(fig_num, "_", fig_label, ".png")),
      plot = p, width = 11, height = 8.5, dpi = 320, bg = "white",
      device = ragg::agg_png
    )
  }
}

save_all_outputs <- function(res, prefix, title,
                             right_label, left_label,
                             out_dir_comloc, out_dir_semloc,
                             volcano_name = "volcano.png",
                             highlight_genes = NULL,
                             fig_num = NULL,
                             fig_label = NULL) {

  # com LOCs (Material Suplementar - sem copia em figuras_tese)
  save_result_tables(
    res = res,
    out_dir = out_dir_comloc,
    prefix = prefix,
    remove_loc = FALSE
  )

  plot_volcano(
    res = res,
    title = title,
    path_plot = file.path(out_dir_comloc, volcano_name),
    right_label = right_label,
    left_label = left_label,
    remove_loc = FALSE,
    highlight_genes = highlight_genes
  )

  # sem LOCs (usado no corpo da tese - com copia em figuras_tese)
  save_result_tables(
    res = res,
    out_dir = out_dir_semloc,
    prefix = prefix,
    remove_loc = TRUE
  )

  plot_volcano(
    res = res,
    title = title,
    path_plot = file.path(out_dir_semloc, volcano_name),
    right_label = right_label,
    left_label = left_label,
    remove_loc = TRUE,
    highlight_genes = highlight_genes,
    fig_num = fig_num,
    fig_label = fig_label
  )
}

# numero da figura no capitulo de Resultados, por comparacao temporal
fig_lookup_tempo <- list(
  "RHDV_4DPI_vs_CONTROL_14DPI"  = "Figure23",
  "RCV_4DPI_vs_CONTROL_14DPI"   = "Figure24",
  "RHDV_vs_RCV_4DPI"            = "Figure25",
  "RHDV_14DPI_vs_CONTROL_14DPI" = "Figure26",
  "RCV_14DPI_vs_CONTROL_14DPI"  = "Figure27",
  "RHDV_vs_RCV_14DPI"           = "Figure28",
  "RHDV_14DPI_vs_4DPI"          = "Figure29",
  "RCV_14DPI_vs_4DPI"           = "Figure30"
)

# genes citados no texto dos Resultados, por comparacao temporal — usados
# para destacar os mesmos genes nos volcano plots do time-course
highlight_tempo <- list(
  "RHDV_14DPI_vs_4DPI" = c("TMEM86A", "FPR1", "FCGR2A", "IL1R2", "CLEC1A", "NREP",
                           "MYO7A", "FCER1G", "RTP4", "ESAM", "XAF1", "OAS1",
                           "TTC39B", "APLNR", "DDX4"),
  "RCV_14DPI_vs_4DPI"  = c("DHRS3", "RAD9B", "PLEKHB2", "SLC4A1", "SPATA9", "DDX4",
                           "CHIT1", "SUCNR1", "SFTPD", "GNG10_1", "S100A9", "POLI"),
  "RHDV_14DPI_vs_CONTROL_14DPI" = c("NOXRED1", "C16H1orf131", "GRO-A", "NLGN3"),
  "RCV_14DPI_vs_CONTROL_14DPI"  = c("RPL34", "DHRS3", "LIPJ", "KLK6", "GNAT3"),
  "RHDV_4DPI_vs_CONTROL_14DPI"  = c("RTP4", "IFIT3", "IFIT2", "RSAD2", "FCGR2A",
                                     "IFIT5", "PARP9", "OAS1", "XAF1", "USP18",
                                     "OAS2", "UBE2L6",
                                     "HSPB11", "NREP", "RIBC1", "MRPS33", "ADRA1A", "AXIN2"),
  "RCV_4DPI_vs_CONTROL_14DPI"   = c("RPL34", "MRPL17", "C16H1orf131", "PLEKHB2",
                                     "RAD9B", "RIBC1", "SPATA9", "CASP4", "IRF7",
                                     "XAF1", "ETV7", "XCR1"),
  "RHDV_vs_RCV_4DPI"  = c("NEU1", "IL1R2", "SLC43A2", "GLMP", "ATP6V0C", "CD68",
                          "SLC7A7", "NREP", "PLVAP", "CTSB", "PFKFB4", "CLEC1A",
                          "FADS2", "AXIN2", "APLNR", "SLC17A8", "SERPINF1", "EFNB3"),
  "RHDV_vs_RCV_14DPI" = c("DHRS3", "ADGRE2", "PLCD4", "TMOD1", "NLGN3", "EPHB1",
                          "FANK1", "CNTN2", "PGC", "PDZRN4")
)

run_deseq_subset <- function(count_data, coldata, subset_idx,
                             design_formula, contrast_vec) {
  
  col_sub <- coldata[subset_idx, , drop = FALSE]
  col_sub <- droplevels(col_sub)
  
  counts_sub <- count_data[, rownames(col_sub), drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = col_sub,
    design = design_formula
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(dds, contrast = contrast_vec)
  res <- res[order(res$padj), ]
  
  return(res)
}

# =========================================================
# 4. TEMPORAL ANALYSIS
# =========================================================

# ---------------------------------------------------------
# 4.1 RHDV 14DPI vs 4DPI
# ---------------------------------------------------------

idx_rhdv_time <- rownames(coldata) %in%
  rownames(subset(coldata, Condition == "RHDV" & Sampling %in% c("4DPI", "14DPI")))

coldata_rhdv_time <- coldata
coldata_rhdv_time$Sampling <- factor(coldata_rhdv_time$Sampling, levels = c("4DPI", "14DPI"))
coldata_rhdv_time$Tissue <- factor(coldata_rhdv_time$Tissue)

res_rhdv_14vs4 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_rhdv_time,
  subset_idx = idx_rhdv_time,
  design_formula = ~ Tissue + Sampling,
  contrast_vec = c("Sampling", "14DPI", "4DPI")
)

save_all_outputs(
  res = res_rhdv_14vs4,
  prefix = "RHDV_14DPI_vs_4DPI",
  title = "RHDV 14DPI vs 4DPI",
  right_label = "14DPI",
  left_label = "4DPI",
  out_dir_comloc = dir_tempo_rhdv_14dpi_comloc,
  out_dir_semloc = dir_tempo_rhdv_14dpi_semloc,
  volcano_name = "volcano_RHDV_14DPI_vs_4DPI.png",
  highlight_genes = highlight_tempo[["RHDV_14DPI_vs_4DPI"]],
  fig_num = fig_lookup_tempo[["RHDV_14DPI_vs_4DPI"]],
  fig_label = "RHDV_14dpi_vs_4dpi_withinRHDV"
)

# ---------------------------------------------------------
# 4.2 RCV 14DPI vs 4DPI
# ---------------------------------------------------------

idx_rcv_time <- rownames(coldata) %in%
  rownames(subset(coldata, Condition == "RCV" & Sampling %in% c("4DPI", "14DPI")))

coldata_rcv_time <- coldata
coldata_rcv_time$Sampling <- factor(coldata_rcv_time$Sampling, levels = c("4DPI", "14DPI"))
coldata_rcv_time$Tissue <- factor(coldata_rcv_time$Tissue)

res_rcv_14vs4 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_rcv_time,
  subset_idx = idx_rcv_time,
  design_formula = ~ Tissue + Sampling,
  contrast_vec = c("Sampling", "14DPI", "4DPI")
)

save_all_outputs(
  res = res_rcv_14vs4,
  prefix = "RCV_14DPI_vs_4DPI",
  title = "RCV 14DPI vs 4DPI",
  right_label = "14DPI",
  left_label = "4DPI",
  out_dir_comloc = dir_tempo_rcv_4dpi_comloc,
  out_dir_semloc = dir_tempo_rcv_4dpi_semloc,
  volcano_name = "volcano_RCV_14DPI_vs_4DPI.png",
  highlight_genes = highlight_tempo[["RCV_14DPI_vs_4DPI"]],
  fig_num = fig_lookup_tempo[["RCV_14DPI_vs_4DPI"]],
  fig_label = "RCV_14dpi_vs_4dpi_withinRCV"
)

# ---------------------------------------------------------
# 4.3 RCV 14DPI vs CONTROL 14DPI
# ---------------------------------------------------------

col_14_rcv_ctrl <- subset(coldata, Sampling == "14DPI" & Condition %in% c("CONTROL", "RCV"))
idx_14_rcv_ctrl <- rownames(coldata) %in% rownames(col_14_rcv_ctrl)

coldata_14_rcv_ctrl <- coldata
coldata_14_rcv_ctrl$Group <- ifelse(
  coldata_14_rcv_ctrl$Condition == "RCV" & coldata_14_rcv_ctrl$Sampling == "14DPI",
  "RCV_14DPI", "CONTROL_14DPI"
)
coldata_14_rcv_ctrl$Group <- factor(
  coldata_14_rcv_ctrl$Group,
  levels = c("CONTROL_14DPI", "RCV_14DPI")
)

res_rcv14_vs_ctrl14 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_14_rcv_ctrl,
  subset_idx = idx_14_rcv_ctrl,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RCV_14DPI", "CONTROL_14DPI")
)

save_all_outputs(
  res = res_rcv14_vs_ctrl14,
  prefix = "RCV_14DPI_vs_CONTROL_14DPI",
  title = "RCV 14DPI vs CONTROL",
  right_label = "RCV 14DPI",
  left_label = "CONTROL",
  out_dir_comloc = dir_tempo_rcv_ctrl_14dpi_comloc,
  out_dir_semloc = dir_tempo_rcv_ctrl_14dpi_semloc,
  volcano_name = "volcano_RCV_14DPI_vs_CONTROL_14DPI.png",
  highlight_genes = highlight_tempo[["RCV_14DPI_vs_CONTROL_14DPI"]],
  fig_num = fig_lookup_tempo[["RCV_14DPI_vs_CONTROL_14DPI"]],
  fig_label = "RCV_vs_Control_14dpi"
)

# ---------------------------------------------------------
# 4.4 RHDV 14DPI vs CONTROL 14DPI
# ---------------------------------------------------------

col_14_rhdv_ctrl <- subset(coldata, Sampling == "14DPI" & Condition %in% c("CONTROL", "RHDV"))
idx_14_rhdv_ctrl <- rownames(coldata) %in% rownames(col_14_rhdv_ctrl)

coldata_14_rhdv_ctrl <- coldata
coldata_14_rhdv_ctrl$Group <- ifelse(
  coldata_14_rhdv_ctrl$Condition == "RHDV" & coldata_14_rhdv_ctrl$Sampling == "14DPI",
  "RHDV_14DPI", "CONTROL_14DPI"
)
coldata_14_rhdv_ctrl$Group <- factor(
  coldata_14_rhdv_ctrl$Group,
  levels = c("CONTROL_14DPI", "RHDV_14DPI")
)

res_rhdv14_vs_ctrl14 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_14_rhdv_ctrl,
  subset_idx = idx_14_rhdv_ctrl,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RHDV_14DPI", "CONTROL_14DPI")
)

save_all_outputs(
  res = res_rhdv14_vs_ctrl14,
  prefix = "RHDV_14DPI_vs_CONTROL_14DPI",
  title = "RHDV 14DPI vs CONTROL",
  right_label = "RHDV 14DPI",
  left_label = "CONTROL",
  out_dir_comloc = dir_tempo_rhdv_ctrl_14dpi_comloc,
  out_dir_semloc = dir_tempo_rhdv_ctrl_14dpi_semloc,
  volcano_name = "volcano_RHDV_14DPI_vs_CONTROL_14DPI.png",
  highlight_genes = highlight_tempo[["RHDV_14DPI_vs_CONTROL_14DPI"]],
  fig_num = fig_lookup_tempo[["RHDV_14DPI_vs_CONTROL_14DPI"]],
  fig_label = "RHDV_vs_Control_14dpi"
)

# ---------------------------------------------------------
# 4.5 RCV 4DPI vs CONTROL 14DPI
# ---------------------------------------------------------

col_4_rcv_ctrl14 <- subset(
  coldata,
  (Condition == "RCV" & Sampling == "4DPI") |
    (Condition == "CONTROL" & Sampling == "14DPI")
)
idx_4_rcv_ctrl14 <- rownames(coldata) %in% rownames(col_4_rcv_ctrl14)

coldata_4_rcv_ctrl14 <- coldata
coldata_4_rcv_ctrl14$Group <- NA
coldata_4_rcv_ctrl14$Group[coldata_4_rcv_ctrl14$Condition == "RCV" &
                             coldata_4_rcv_ctrl14$Sampling == "4DPI"] <- "RCV_4DPI"
coldata_4_rcv_ctrl14$Group[coldata_4_rcv_ctrl14$Condition == "CONTROL" &
                             coldata_4_rcv_ctrl14$Sampling == "14DPI"] <- "CONTROL_14DPI"
coldata_4_rcv_ctrl14$Group <- factor(
  coldata_4_rcv_ctrl14$Group,
  levels = c("CONTROL_14DPI", "RCV_4DPI")
)

res_rcv4_vs_ctrl14 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_4_rcv_ctrl14,
  subset_idx = idx_4_rcv_ctrl14,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RCV_4DPI", "CONTROL_14DPI")
)

save_all_outputs(
  res = res_rcv4_vs_ctrl14,
  prefix = "RCV_4DPI_vs_CONTROL_14DPI",
  title = "RCV 4DPI vs CONTROL",
  right_label = "RCV 4DPI",
  left_label = "CONTROL",
  out_dir_comloc = dir_tempo_rcv_ctrl_4dpi_comloc,
  out_dir_semloc = dir_tempo_rcv_ctrl_4dpi_semloc,
  volcano_name = "volcano_RCV_4DPI_vs_CONTROL_14DPI.png",
  highlight_genes = highlight_tempo[["RCV_4DPI_vs_CONTROL_14DPI"]],
  fig_num = fig_lookup_tempo[["RCV_4DPI_vs_CONTROL_14DPI"]],
  fig_label = "RCV_vs_Control_4dpi"
)

# ---------------------------------------------------------
# 4.6 RHDV 4DPI vs CONTROL 14DPI
# ---------------------------------------------------------

col_4_rhdv_ctrl14 <- subset(
  coldata,
  (Condition == "RHDV" & Sampling == "4DPI") |
    (Condition == "CONTROL" & Sampling == "14DPI")
)
idx_4_rhdv_ctrl14 <- rownames(coldata) %in% rownames(col_4_rhdv_ctrl14)

coldata_4_rhdv_ctrl14 <- coldata
coldata_4_rhdv_ctrl14$Group <- NA
coldata_4_rhdv_ctrl14$Group[coldata_4_rhdv_ctrl14$Condition == "RHDV" &
                              coldata_4_rhdv_ctrl14$Sampling == "4DPI"] <- "RHDV_4DPI"
coldata_4_rhdv_ctrl14$Group[coldata_4_rhdv_ctrl14$Condition == "CONTROL" &
                              coldata_4_rhdv_ctrl14$Sampling == "14DPI"] <- "CONTROL_14DPI"
coldata_4_rhdv_ctrl14$Group <- factor(
  coldata_4_rhdv_ctrl14$Group,
  levels = c("CONTROL_14DPI", "RHDV_4DPI")
)

res_rhdv4_vs_ctrl14 <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_4_rhdv_ctrl14,
  subset_idx = idx_4_rhdv_ctrl14,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RHDV_4DPI", "CONTROL_14DPI")
)

save_all_outputs(
  res = res_rhdv4_vs_ctrl14,
  prefix = "RHDV_4DPI_vs_CONTROL_14DPI",
  title = "RHDV 4DPI vs CONTROL",
  right_label = "RHDV 4DPI",
  left_label = "CONTROL",
  out_dir_comloc = dir_tempo_rhdv_ctrl_4dpi_comloc,
  out_dir_semloc = dir_tempo_rhdv_ctrl_4dpi_semloc,
  volcano_name = "volcano_RHDV_4DPI_vs_CONTROL_14DPI.png",
  highlight_genes = highlight_tempo[["RHDV_4DPI_vs_CONTROL_14DPI"]],
  fig_num = fig_lookup_tempo[["RHDV_4DPI_vs_CONTROL_14DPI"]],
  fig_label = "RHDV_vs_Control_4dpi"
)

# ---------------------------------------------------------
# 4.7 RHDV vs RCV em 4DPI
# ---------------------------------------------------------

col_4_virus <- subset(coldata, Sampling == "4DPI" & Condition %in% c("RCV", "RHDV"))
idx_4_virus <- rownames(coldata) %in% rownames(col_4_virus)

coldata_4_virus <- coldata
coldata_4_virus$Group <- NA
coldata_4_virus$Group[coldata_4_virus$Condition == "RCV" &
                        coldata_4_virus$Sampling == "4DPI"] <- "RCV_4DPI"
coldata_4_virus$Group[coldata_4_virus$Condition == "RHDV" &
                        coldata_4_virus$Sampling == "4DPI"] <- "RHDV_4DPI"
coldata_4_virus$Group <- factor(
  coldata_4_virus$Group,
  levels = c("RCV_4DPI", "RHDV_4DPI")
)

res_rhdv_vs_rcv_4dpi <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_4_virus,
  subset_idx = idx_4_virus,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RHDV_4DPI", "RCV_4DPI")
)

save_all_outputs(
  res = res_rhdv_vs_rcv_4dpi,
  prefix = "RHDV_vs_RCV_4DPI",
  title = "RHDV vs RCV at 4DPI",
  right_label = "RHDV 4DPI",
  left_label = "RCV 4DPI",
  out_dir_comloc = dir_tempo_rhdv_rcv_comloc,
  out_dir_semloc = dir_tempo_rhdv_rcv_semloc,
  volcano_name = "volcano_RHDV_vs_RCV_4DPI.png",
  highlight_genes = highlight_tempo[["RHDV_vs_RCV_4DPI"]],
  fig_num = fig_lookup_tempo[["RHDV_vs_RCV_4DPI"]],
  fig_label = "RHDV_vs_RCV_4dpi"
)

# ---------------------------------------------------------
# 4.8 RHDV vs RCV em 14DPI
# ---------------------------------------------------------

col_14_virus <- subset(coldata, Sampling == "14DPI" & Condition %in% c("RCV", "RHDV"))
idx_14_virus <- rownames(coldata) %in% rownames(col_14_virus)

coldata_14_virus <- coldata
coldata_14_virus$Group <- NA
coldata_14_virus$Group[coldata_14_virus$Condition == "RCV" &
                         coldata_14_virus$Sampling == "14DPI"] <- "RCV_14DPI"
coldata_14_virus$Group[coldata_14_virus$Condition == "RHDV" &
                         coldata_14_virus$Sampling == "14DPI"] <- "RHDV_14DPI"
coldata_14_virus$Group <- factor(
  coldata_14_virus$Group,
  levels = c("RCV_14DPI", "RHDV_14DPI")
)

res_rhdv_vs_rcv_14dpi <- run_deseq_subset(
  count_data = count_data,
  coldata = coldata_14_virus,
  subset_idx = idx_14_virus,
  design_formula = ~ Tissue + Group,
  contrast_vec = c("Group", "RHDV_14DPI", "RCV_14DPI")
)

save_all_outputs(
  res = res_rhdv_vs_rcv_14dpi,
  prefix = "RHDV_vs_RCV_14DPI",
  title = "RHDV vs RCV at 14DPI",
  right_label = "RHDV 14DPI",
  left_label = "RCV 14DPI",
  out_dir_comloc = dir_tempo_rhdv_rcv_comloc,
  out_dir_semloc = dir_tempo_rhdv_rcv_semloc,
  volcano_name = "volcano_RHDV_vs_RCV_14DPI.png",
  highlight_genes = highlight_tempo[["RHDV_vs_RCV_14DPI"]],
  fig_num = fig_lookup_tempo[["RHDV_vs_RCV_14DPI"]],
  fig_label = "RHDV_vs_RCV_14dpi"
)