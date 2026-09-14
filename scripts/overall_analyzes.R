library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# pasta com copia de todas as figuras usadas na tese, nomeadas pelo numero
# de figura do capitulo de Resultados (facil de saber qual grafico e qual)
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

coldata$Condition <- factor(
  toupper(trimws(coldata$Condition)),
  levels = c("CONTROL", "RCV", "RHDV")
)
coldata$Tissue <- factor(trimws(coldata$Tissue))

# covariavel de batch: sequencing_run (2 niveis, A/B) em vez de Lane (6
# niveis) - Lane fica confundida com Condition dentro de cada tecido (por
# ex. no Duodeno a lane A nunca aparece em CONTROL), o que faz o DESeq2
# recusar o modelo (model matrix not full rank). sequencing_run esta bem
# mais equilibrada entre tecidos e condicoes.
coldata$sequencing_run <- factor(trimws(coldata$sequencing_run))

dds_all <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~ sequencing_run + Tissue + Condition
)

keep <- rowSums(counts(dds_all) >= 10) >= 3
dds_all <- dds_all[keep, ]
dds_all <- DESeq(dds_all)

res_rhdv_vs_control <- results(dds_all, contrast = c("Condition", "RHDV", "CONTROL"))
res_rcv_vs_control  <- results(dds_all, contrast = c("Condition", "RCV", "CONTROL"))
res_rcv_vs_rhdv     <- results(dds_all, contrast = c("Condition", "RCV", "RHDV"))

# ordenar por padj ascendente, seguindo a convenção de Love et al. (2014);
# esta ordenação propaga-se automaticamente para as tabelas completas e para
# os subconjuntos de significativos (sig_*) derivados abaixo via subset()
res_rhdv_vs_control <- res_rhdv_vs_control[order(res_rhdv_vs_control$padj), ]
res_rcv_vs_control  <- res_rcv_vs_control[order(res_rcv_vs_control$padj), ]
res_rcv_vs_rhdv     <- res_rcv_vs_rhdv[order(res_rcv_vs_rhdv$padj), ]

base_path <- "/Users/mariliaandrade/Desktop/Tese/resultados_final/geral"

dir_rhdv_ctrl <- file.path(base_path, "RHDV_CONTROL")
dir_rcv_ctrl  <- file.path(base_path, "RCV_CONTROL")
dir_rhdv_rcv  <- file.path(base_path, "RHDV_RCV")

dir.create(dir_rhdv_ctrl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rcv_ctrl,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rhdv_rcv,  recursive = TRUE, showWarnings = FALSE)

# salvar resultados completos
write.csv(
  as.data.frame(res_rhdv_vs_control),
  file.path(dir_rhdv_ctrl, "RHDV_vs_CONTROL.csv")
)

write.csv(
  as.data.frame(res_rcv_vs_control),
  file.path(dir_rcv_ctrl, "RCV_vs_CONTROL.csv")
)

write.csv(
  as.data.frame(res_rcv_vs_rhdv),
  file.path(dir_rhdv_rcv, "RCV_vs_RHDV.csv")
)

# significativos
sig_rhdv <- subset(
  as.data.frame(res_rhdv_vs_control),
  padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1
)

sig_rcv <- subset(
  as.data.frame(res_rcv_vs_control),
  padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1
)

sig_rcv_rhdv <- subset(
  as.data.frame(res_rcv_vs_rhdv),
  padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1
)

# NAO remover genes LOC aqui: estes sig_* completos (com_locs) sao a fonte
# do Material Suplementar. A versao usada no corpo da tese (sem_locs) e
# criada no bloco seguinte.
# sig_rhdv <- sig_rhdv[!grepl("^LOC", rownames(sig_rhdv)), ]
# sig_rcv <- sig_rcv[!grepl("^LOC", rownames(sig_rcv)), ]
# sig_rcv_rhdv <- sig_rcv_rhdv[!grepl("^LOC", rownames(sig_rcv_rhdv)), ]

# salvar significativos (com_locs, para Material Suplementar)
write.csv(
  sig_rhdv,
  file.path(dir_rhdv_ctrl, "RHDV_vs_CONTROL_sig.csv")
)

write.csv(
  sig_rcv,
  file.path(dir_rcv_ctrl, "RCV_vs_CONTROL_sig.csv")
)

write.csv(
  sig_rcv_rhdv,
  file.path(dir_rhdv_rcv, "RCV_vs_RHDV_sig.csv")
)

# ---- versao sem_locs: criterio oficial usado no corpo da tese ----
# (excluindo loci LOC nao caracterizados; excecao: LOCs entretanto
# renomeados no NCBI devem ser reportados com o simbolo atual, ver
# LOC_genes_renamed_NCBI.csv)
dir_rhdv_ctrl_sl <- file.path(dir_rhdv_ctrl, "sem_locs")
dir_rcv_ctrl_sl  <- file.path(dir_rcv_ctrl, "sem_locs")
dir_rhdv_rcv_sl  <- file.path(dir_rhdv_rcv, "sem_locs")
dir.create(dir_rhdv_ctrl_sl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rcv_ctrl_sl,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rhdv_rcv_sl,  recursive = TRUE, showWarnings = FALSE)

sig_rhdv_sem_loc     <- sig_rhdv[!grepl("^LOC", rownames(sig_rhdv)), ]
sig_rcv_sem_loc      <- sig_rcv[!grepl("^LOC", rownames(sig_rcv)), ]
sig_rcv_rhdv_sem_loc <- sig_rcv_rhdv[!grepl("^LOC", rownames(sig_rcv_rhdv)), ]

write.csv(
  sig_rhdv_sem_loc,
  file.path(dir_rhdv_ctrl_sl, "RHDV_vs_CONTROL_sig.csv")
)

write.csv(
  sig_rcv_sem_loc,
  file.path(dir_rcv_ctrl_sl, "RCV_vs_CONTROL_sig.csv")
)

write.csv(
  sig_rcv_rhdv_sem_loc,
  file.path(dir_rhdv_rcv_sl, "RCV_vs_RHDV_sig.csv")
)

# tabela com contagem de genes up/down (com_locs, mantida por compatibilidade)
deg_counts <- data.frame(
  Comparison = c("RHDV vs CONTROL", "RCV vs CONTROL", "RCV vs RHDV"),
  Up = c(
    sum(sig_rhdv$log2FoldChange > 1, na.rm = TRUE),
    sum(sig_rcv$log2FoldChange > 1, na.rm = TRUE),
    sum(sig_rcv_rhdv$log2FoldChange > 1, na.rm = TRUE)
  ),
  Down = c(
    sum(sig_rhdv$log2FoldChange < -1, na.rm = TRUE),
    sum(sig_rcv$log2FoldChange < -1, na.rm = TRUE),
    sum(sig_rcv_rhdv$log2FoldChange < -1, na.rm = TRUE)
  )
)

# tabela oficial (sem_locs) usada no corpo da tese (Tabela 1 - nivel global)
deg_counts_sem_locs <- data.frame(
  Comparison = c("RHDV vs CONTROL", "RCV vs CONTROL", "RCV vs RHDV"),
  Total = c(nrow(sig_rhdv_sem_loc), nrow(sig_rcv_sem_loc), nrow(sig_rcv_rhdv_sem_loc)),
  Up = c(
    sum(sig_rhdv_sem_loc$log2FoldChange > 1, na.rm = TRUE),
    sum(sig_rcv_sem_loc$log2FoldChange > 1, na.rm = TRUE),
    sum(sig_rcv_rhdv_sem_loc$log2FoldChange > 1, na.rm = TRUE)
  ),
  Down = c(
    sum(sig_rhdv_sem_loc$log2FoldChange < -1, na.rm = TRUE),
    sum(sig_rcv_sem_loc$log2FoldChange < -1, na.rm = TRUE),
    sum(sig_rcv_rhdv_sem_loc$log2FoldChange < -1, na.rm = TRUE)
  )
)

write.csv(
  deg_counts_sem_locs,
  file.path(base_path, "DEG_counts_sem_locs.csv"),
  row.names = FALSE
)

write.csv(
  deg_counts,
  file.path(base_path, "DEG_counts.csv"),
  row.names = FALSE
)

# Figura de abertura do capitulo: PCA das amostras (nao repete a Tabela 1,
# mostra antes a estrutura dos dados - normalmente o tecido domina a
# variacao, o que justifica visualmente a decisao metodologica de analisar
# cada comparacao tambem por tecido, e nao so no modelo global).
vsd_all <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd_all, intgroup = c("Condition", "Tissue"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# nomes das amostras (coluna "name" devolvida por plotPCA) para conseguires
# identificar exatamente quais sao os pontos afastados do respetivo grupo
# de tecido - sem isto nao da para saber quais animais sao os outliers
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Tissue)) +
  geom_point(size = 3.5, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = name),
    size = 2.6,
    color = "black",
    bg.color = "white",
    bg.r = 0.12,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.15,
    force = 2,
    seed = 42,
    segment.size = 0.25,
    segment.color = "grey40",
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("CONTROL" = "grey40", "RCV" = "#1F5FB2", "RHDV" = "#D62728")) +
  labs(
    title = "Sample Clustering (PCA, all tissues and conditions)",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Condition",
    shape = "Tissue"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(base_path, "Sample_PCA.png"),
  plot = p_pca,
  width = 10,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

# copia com o numero da figura usado no capitulo de Resultados
ggsave(
  file.path(figuras_tese_dir, "Figure08_Overview_Sample_PCA.png"),
  plot = p_pca,
  width = 10,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

plot_volcano <- function(res, title, path_plot,
                         right_label, left_label,
                         remove_loc = TRUE,
                         lfc_min = -3, lfc_max = 3,
                         padj_min = 0.05,
                         point_size = 2.5,
                         highlight_genes = NULL,
                         top_n_label = 10,
                         fig_num = NULL,
                         fig_label = NULL) {

  df <- as.data.frame(res)
  df$Gene <- rownames(df)

  if (remove_loc) {
    df <- df[!grepl("^LOC", df$Gene), ]
  }

  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]
  df$negLogPadj <- -log10(pmax(df$padj, .Machine$double.xmin))

  df$Significance <- "Not Sig"
  df$Significance[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$Significance[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"

  # genes citados no texto dos Resultados: sempre rotulados e destacados
  # com contorno preto
  df$Highlight <- df$Gene %in% highlight_genes
  df_highlight <- df[df$Highlight, ]

  # rotular APENAS os genes destacados no texto; se nao houver nenhum,
  # rotular só os top_n_label mais significativos - nunca a nuvem toda,
  # para o grafico manter-se legivel
  if (nrow(df_highlight) > 0) {
    # preservar a ordem do vetor highlight_genes (a ordem em que os genes
    # sao discutidos no texto), nao a ordem do data frame
    ord <- match(highlight_genes, df_highlight$Gene)
    ord <- ord[!is.na(ord)]
    df_label <- df_highlight[ord, ]
  } else {
    sig_df <- df[df$Significance != "Not Sig", ]
    sig_df <- sig_df[order(sig_df$padj), ]
    df_label <- head(sig_df, top_n_label)
  }

  # contagem total de significativos (a lista completa continua disponivel
  # no ficheiro *_sig.csv correspondente)
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
    labs(
      title = title,
      subtitle = paste0("← Higher in ", left_label,
                         "                              Higher in ", right_label, " →"),
      caption = paste0("Total DEGs: ", n_up + n_down,
                        " (Up: ", n_up, ", Down: ", n_down,
                        "). Labelled genes are those discussed in the text; ",
                        "the full list is provided in the corresponding sig.csv / Supplementary Table."),
      x = "Log2 Fold Change",
      y = "-Log10(adjusted p-value)",
      color = "Expression"
    ) +
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

  # bg = "white" força fundo opaco no PNG (o tema minimal deixa o fundo
  # transparente, o que fazia os graficos aparecerem pretos/ilegiveis em
  # visualizadores com modo escuro)
  ggsave(path_plot, plot = p, width = 11, height = 8.5, dpi = 320, bg = "white")

  # copia com o numero da figura usado no capitulo de Resultados, para
  # facilitar identificar qual grafico e qual (sobretudo no time-course)
  if (!is.null(fig_num)) {
    ggsave(
      file.path(figuras_tese_dir, paste0(fig_num, "_", fig_label, ".png")),
      plot = p, width = 11, height = 8.5, dpi = 320, bg = "white"
    )
  }

  return(p)
}

# genes citados no texto dos Resultados para cada comparação global —
# mantidos aqui para que o destaque no gráfico acompanhe sempre o texto
# NOTA: os 8 genes "Down" (ADIPOQ, KCNQ3, NLGN3, KLB, SLC26A9, COL6A6, MC5R,
# THRSP) foram adicionados para que o lado azul do vulcano também fique
# identificado — confirmados a partir de RHDV_CONTROL/sem_locs/
# RHDV_vs_CONTROL_sig.csv (execução de 18/08, 135 total = 127 up + 8 down).
# Se ainda não estiverem mencionados no texto de Resultados, acrescentar uma
# frase equivalente à usada para outros genes sem função antiviral descrita.
highlight_rhdv_control <- c("CD163", "NOXRED1", "ISG15", "RSAD2", "CCL5", "S100A9",
                             "VSIG4", "CX3CL1", "IRF7", "LGALS9", "CASP4",
                             "ADIPOQ", "KCNQ3", "NLGN3", "KLB", "SLC26A9",
                             "COL6A6", "MC5R", "THRSP")
highlight_rcv_control  <- c("RPL34", "CHGA", "FAM212B", "PCDH20", "KLK6", "C13H1orf189")
highlight_rcv_rhdv     <- c("CX3CL1", "SLC11A1", "CD163", "GPNMB", "GPX3", "S100A12",
                             "CCL2", "S100A9", "CCL5", "THBS4", "CNTN2", "VSIG4",
                             "MSR1", "CCR1", "ISG15", "RSAD2", "PTX3", "SERPINE1")

plot_volcano(
  res_rhdv_vs_control,
  "RHDV vs Control",
  file.path(dir_rhdv_ctrl, "volcano_RHDV_vs_CONTROL.png"),
  right_label = "RHDV",
  left_label = "Control",
  remove_loc = TRUE,
  point_size = 3,
  lfc_min = -3,
  lfc_max = 3,
  highlight_genes = highlight_rhdv_control,
  fig_num = "Figure10",
  fig_label = "Global_RHDV_vs_Control"
)

plot_volcano(
  res_rcv_vs_control,
  "RCV vs Control",
  file.path(dir_rcv_ctrl, "volcano_RCV_vs_CONTROL.png"),
  right_label = "RCV",
  left_label = "Control",
  remove_loc = TRUE,
  point_size = 3,
  lfc_min = -3,
  lfc_max = 3,
  highlight_genes = highlight_rcv_control,
  fig_num = "Figure15",
  fig_label = "Global_RCV_vs_Control"
)

plot_volcano(
  res_rcv_vs_rhdv,
  "RHDV vs RCV",
  file.path(dir_rhdv_rcv, "volcano_RCV_vs_RHDV.png"),
  right_label = "RCV",
  left_label = "RHDV",
  remove_loc = TRUE,
  point_size = 3,
  lfc_min = -3,
  lfc_max = 3,
  highlight_genes = highlight_rcv_rhdv,
  fig_num = "Figure18",
  fig_label = "Global_RCV_vs_RHDV"
)

plot_lollipop_sig <- function(res, title, path_plot,
                              padj_cutoff = 0.05,
                              lfc_cutoff = 0,
                              remove_loc = TRUE,
                              top_n = NULL) {
  
  df <- as.data.frame(res)
  df$Gene <- rownames(df)
  
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]
  
  if (remove_loc) {
    df <- df[!grepl("^LOC", df$Gene), ]
  }
  
  df$negLogPadj <- -log10(pmax(df$padj, .Machine$double.xmin))
  
  df_sig <- df %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  if (!is.null(top_n)) {
    df_sig <- df_sig %>%
      arrange(padj) %>%
      slice_head(n = top_n)
  }
  
  df_sig <- df_sig %>%
    mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
    arrange(log2FoldChange) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  p <- ggplot(df_sig, aes(x = log2FoldChange, y = Gene, color = Direction)) +
    geom_segment(
      aes(x = 0, xend = log2FoldChange, y = Gene, yend = Gene),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    geom_point(aes(size = negLogPadj)) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "Gene",
      color = "Direction",
      size = "-log10(padj)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  ggsave(
    path_plot,
    p,
    width = 10,
    height = max(6, 0.25 * nrow(df_sig)),
    dpi = 300,
    bg = "white"
  )

  return(p)
}

plot_lollipop_sig(
  res = res_rcv_vs_control,
  title = "RCV vs Control",
  path_plot = file.path(dir_rcv_ctrl, "lollipop_RCV_vs_CONTROL.png"),
  padj_cutoff = 0.05,
  lfc_cutoff = 0,
  remove_loc = TRUE,
  top_n = NULL
)

plot_lollipop_sig(
  res = res_rhdv_vs_control,
  title = "RHDV vs Control",
  path_plot = file.path(dir_rhdv_ctrl, "lollipop_RHDV_vs_CONTROL.png"),
  padj_cutoff = 0.05,
  lfc_cutoff = 0,
  remove_loc = TRUE,
  top_n = NULL
)

plot_lollipop_sig(
  res = res_rcv_vs_rhdv,
  title = "RCV vs RHDV",
  path_plot = file.path(dir_rhdv_rcv, "lollipop_RCV_vs_RHDV.png"),
  padj_cutoff = 0.05,
  lfc_cutoff = 0,
  remove_loc = TRUE,
  top_n = NULL
)

# ##resumo
# save_DEGs <- function(res, path_csv, n_top = 50, padj_cut = 0.05, lfc_cut = 1) {
#   df <- as.data.frame(res)
#   df$Gene <- rownames(df)
#   df <- df[!is.na(df$padj), ]
#   df <- df[!grepl("^LOC", df$Gene), ]
#   
#   sig <- subset(df, padj < padj_cut & abs(log2FoldChange) > lfc_cut)
#   
#   sig <- sig[order(sig$padj, -abs(sig$log2FoldChange), -sig$baseMean), ]
#   
#   top <- head(sig, n_top)
#   
#   write.csv(top, path_csv, row.names = TRUE)
#   
#   return(top)
# }



