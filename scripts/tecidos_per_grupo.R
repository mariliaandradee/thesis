library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

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

# remover duplicadas 181s e 183 no count_data
count_data <- count_data[, !duplicated(colnames(count_data))]

# alinhar a ordem do coldata ao count_data
coldata <- coldata[colnames(count_data), , drop = FALSE] 

#confirmar que estão iguais nos dois
all(rownames(coldata) == colnames(count_data))

base_plot <- "/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo"

dirs <- c("RHDV_CONTROL", "RCV_CONTROL", "RHDV_RCV")

for (d in dirs) {
  dir.create(file.path(base_plot, d), recursive = TRUE, showWarnings = FALSE)
}
#Função dos tecidos, no design mantive apenas a condição pois ja chamei os tecidos no col_subset
results_list <- list()
tissues <- unique(coldata$Tissue)

for (tissue in tissues) {
  
  col_subset <- coldata[coldata$Tissue == tissue, ]
  counts_subset <- count_data[, rownames(col_subset)]
  
  col_subset$Condition <- factor(col_subset$Condition,
                                 levels = c("CONTROL", "RCV", "RHDV"))

  # sequencing_run como covariavel de batch (ver overall_analyzes.R para a
  # justificacao - Lane fica confundida com Condition dentro de cada tecido,
  # sequencing_run (2 niveis) nao)
  col_subset$sequencing_run <- factor(trimws(col_subset$sequencing_run))

  dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                colData = col_subset,
                                design = ~ sequencing_run + Condition)
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res1 <- results(dds, contrast = c("Condition", "RHDV", "CONTROL"))
  res2 <- results(dds, contrast = c("Condition", "RCV", "CONTROL"))
  res3 <- results(dds, contrast = c("Condition", "RCV", "RHDV"))
  
  results_list[[tissue]] <- list(
    RHDV_vs_CONTROL = res1,
    RCV_vs_CONTROL  = res2,
    RCV_vs_RHDV     = res3
  )
}

base_plot <- "/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo"

dirs <- c("RHDV_CONTROL/sem_locs", "RCV_CONTROL/sem_locs", "RHDV_RCV/sem_locs",
          "RHDV_CONTROL/com_locs_current", "RCV_CONTROL/com_locs_current", "RHDV_RCV/com_locs_current")

for (d in dirs) {
  dir.create(file.path(base_plot, d), recursive = TRUE, showWarnings = FALSE)
}
# acumulador da Tabela 1 (nivel tecido, sem_locs) - preenchido no loop abaixo
deg_summary_tecidos <- data.frame(
  Comparison = character(0), Tissue = character(0),
  Total = integer(0), Up = integer(0), Down = integer(0)
)

count_updown <- function(sig_df) {
  data.frame(
    Total = nrow(sig_df),
    Up = sum(sig_df$log2FoldChange > 1, na.rm = TRUE),
    Down = sum(sig_df$log2FoldChange < -1, na.rm = TRUE)
  )
}

for (tissue in names(results_list)) {

  res1 <- results_list[[tissue]]$RHDV_vs_CONTROL
  res2 <- results_list[[tissue]]$RCV_vs_CONTROL
  res3 <- results_list[[tissue]]$RCV_vs_RHDV

  sig_res1 <- subset(as.data.frame(res1),
                     padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1)

  sig_res2 <- subset(as.data.frame(res2),
                     padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1)

  sig_res3 <- subset(as.data.frame(res3),
                     padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1)

  # COM LOC (versao completa, incluindo loci LOC-prefixados) - fonte para o
  # Material Suplementar; a pasta antiga "com_locs" ficou desatualizada
  # (nao era regravada por este script), por isso grava-se aqui numa pasta
  # nova para nao gerar confusao entre versao antiga e atual.
  write.csv(sig_res1,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RHDV_CONTROL/com_locs_current/",
                   tissue, "_RHDV_vs_CONTROL_sig.csv"))

  write.csv(sig_res2,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RCV_CONTROL/com_locs_current/",
                   tissue, "_RCV_vs_CONTROL_sig.csv"))

  write.csv(sig_res3,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RHDV_RCV/com_locs_current/",
                   tissue, "_RCV_vs_RHDV_sig.csv"))

  # SEM LOC
  sig_res1_sem_loc <- sig_res1[!grepl("^LOC", rownames(sig_res1)), ]
  sig_res2_sem_loc <- sig_res2[!grepl("^LOC", rownames(sig_res2)), ]
  sig_res3_sem_loc <- sig_res3[!grepl("^LOC", rownames(sig_res3)), ]

  write.csv(sig_res1_sem_loc,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RHDV_CONTROL/sem_locs/",
                   tissue, "_RHDV_vs_CONTROL_sig.csv"))

  write.csv(sig_res2_sem_loc,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RCV_CONTROL/sem_locs/",
                   tissue, "_RCV_vs_CONTROL_sig.csv"))

  write.csv(sig_res3_sem_loc,
            paste0("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/RHDV_RCV/sem_locs/",
                   tissue, "_RCV_vs_RHDV_sig.csv"))

  deg_summary_tecidos <- rbind(
    deg_summary_tecidos,
    cbind(Comparison = "RHDV vs CONTROL", Tissue = tissue, count_updown(sig_res1_sem_loc)),
    cbind(Comparison = "RCV vs CONTROL",  Tissue = tissue, count_updown(sig_res2_sem_loc)),
    cbind(Comparison = "RCV vs RHDV",     Tissue = tissue, count_updown(sig_res3_sem_loc))
  )
}

# Tabela 1 (nivel tecido) - fonte unica para o corpo da tese, gerada direto
# do codigo a partir dos ficheiros sem_locs
write.csv(
  deg_summary_tecidos,
  file.path(base_plot, "DEG_summary_tecidos_sem_locs.csv"),
  row.names = FALSE
)


##Volcano plot per tecidos e grupos
base_plot <- "/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo"

dirs <- c("RHDV_CONTROL/sem_locs", "RCV_CONTROL/sem_locs", "RHDV_RCV/sem_locs")

for (d in dirs) {
  dir.create(file.path(base_plot, d), recursive = TRUE, showWarnings = FALSE)
}

plot_volcano <- function(res, title, path_plot,
                         right_label, left_label,
                         lfc_min = -3, lfc_max = 3,
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

  # contagem total de significativos (a lista completa continua disponivel
  # no ficheiro *_sig.csv correspondente, so os genes destacados no texto
  # sao rotulados individualmente no grafico)
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

  # bg = "white" força fundo opaco no PNG (theme_minimal deixa o fundo
  # transparente, o que fazia os graficos aparecerem pretos/ilegiveis em
  # visualizadores com modo escuro - era a causa principal do texto
  # "impossivel de ler")
  ggsave(path_plot, plot = p, width = 11, height = 8.5, dpi = 320, bg = "white")

  if (!is.null(fig_num)) {
    ggsave(
      file.path(figuras_tese_dir, paste0(fig_num, "_", fig_label, ".png")),
      plot = p, width = 11, height = 8.5, dpi = 320, bg = "white"
    )
  }
}

# genes citados no texto dos Resultados, por tecido e comparacao (RHDV vs
# Control e RCV vs RHDV); usados para destacar os mesmos genes no grafico
highlight_lookup <- list(
  "Duodenum.RHDV_vs_CONTROL" = c("RSAD2", "MSMB", "IFIT2", "CMPK2", "ENPP3",
                                  "SERPINI2", "SPP1", "CFB", "IFIT3"),
  "Liver.RHDV_vs_CONTROL"    = c("ORM1", "SLC18A1", "SLC51B", "MAL", "S100A8",
                                  "SLC51A", "S100A9", "C17H15orf48", "FAM169B"),
  "Spleen.RHDV_vs_CONTROL"   = c("TRABD2B", "MC5R", "ADIPOQ", "VWF", "CIDEC",
                                  "CYP2B4", "LEP", "AUNIP", "NOXRED1"),
  "Thymus.RHDV_vs_CONTROL"   = character(0),  # 0 DEGs
  "Duodenum.RCV_vs_RHDV"     = c("MSMB", "GKN2", "CCL20", "USP25", "C17H14orf37",
                                 "PGC", "THBS4", "WSCD2", "CASQ1", "ACTN2", "BARX1",
                                 "CPA4", "VTCN1", "SERPINI2", "WNT9B", "CPA1",
                                 "PLA2G1B", "CPB1", "PNLIP"),
  "Liver.RCV_vs_RHDV"        = c("ORM1", "SLC51B", "SLC11A1", "SLC18A1", "CXCL11",
                                  "CCL5", "CXCL8", "ANKRD1", "S100A9", "SAA3",
                                  "ACVR1C", "MPP4"),
  "Spleen.RCV_vs_RHDV"       = c("ANPEP", "ADGRE2", "VWA2", "LURAP1", "KCNK10",
                                  "PAQR9", "CLEC4A"),
  "Thymus.RCV_vs_RHDV"       = c("IFGGA1", "TIMD4", "GRB7", "IL1R2", "PDK4",
                                  "CYP2A11_1", "TRIM63", "ISG15", "ANGPTL4"),
  "Duodenum.RCV_vs_CONTROL"  = c("USP25", "IFIT2", "RSAD2", "KRT20", "CMPK2",
                                  "CBLN3", "CFB", "IFIT3", "CCL20", "USP18",
                                  "HERC5", "CR1", "SPP1"),
  "Liver.RCV_vs_CONTROL"     = character(0),  # 0 DEGs
  "Spleen.RCV_vs_CONTROL"    = c("MC5R", "ADIPOQ", "LEP", "KLB", "THRSP", "PFKFB1",
                                  "COL6A6", "VAT1L", "ACSBG2"),
  "Thymus.RCV_vs_CONTROL"    = c("IFGGA1", "COL21A1")  # 2 DEGs
)

# numero da figura no capitulo de Resultados, por tecido/comparacao (tecido
# alvo listado primeiro dentro de cada comparacao, depois os nao-alvo)
fig_lookup <- list(
  "Liver.RHDV_vs_CONTROL"    = "Figure11",  # tecido-alvo do RHDV
  "Spleen.RHDV_vs_CONTROL"   = "Figure12",
  "Duodenum.RHDV_vs_CONTROL" = "Figure13",
  "Thymus.RHDV_vs_CONTROL"   = "Figure14",
  "Duodenum.RCV_vs_CONTROL"  = "Figure16",  # tecido-alvo do RCV
  "Spleen.RCV_vs_CONTROL"    = "Figure17",
  "Duodenum.RCV_vs_RHDV"     = "Figure19",  # tecido-alvo do RCV
  "Liver.RCV_vs_RHDV"        = "Figure20",  # tecido-alvo do RHDV
  "Spleen.RCV_vs_RHDV"       = "Figure21",
  "Thymus.RCV_vs_RHDV"       = "Figure22"
)

# IMPORTANTE: estas 3 chamadas tem de estar dentro de um loop por tecido.
# Antes nao estavam - corriam so UMA vez, usando o "tissue"/"res1"/"res3"
# deixados para tras pelo loop anterior (o ultimo tecido processado ali),
# por isso so saia figura de um tecido (Duodenum). Agora corre para os 4.
for (tissue in names(results_list)) {

  res1 <- results_list[[tissue]]$RHDV_vs_CONTROL
  res2 <- results_list[[tissue]]$RCV_vs_CONTROL
  res3 <- results_list[[tissue]]$RCV_vs_RHDV

  plot_volcano(
    res1,
    paste(tissue, "- RHDV vs Control"),
    file.path(base_plot, "RHDV_CONTROL", "sem_locs", "Volcano",
              paste0(tissue, "_volcano_RHDV_vs_CONTROL.png")),
    right_label = "RHDV",
    left_label = "Control",
    point_size = 3,
    lfc_min = -3,
    lfc_max = 3,
    max.overlaps = 50,
    highlight_genes = highlight_lookup[[paste0(tissue, ".RHDV_vs_CONTROL")]],
    fig_num = fig_lookup[[paste0(tissue, ".RHDV_vs_CONTROL")]],
    fig_label = paste0(tissue, "_RHDV_vs_Control")
  )

  # nota: esta versao de RCV vs Control e substituida pela versao com
  # regras especiais por tecido, mais abaixo no script (zoom no Duodeno,
  # etc.) - mantida aqui so para gerar o ficheiro base, sem fig_num
  plot_volcano(
    res2,
    paste(tissue, "- RCV vs Control"),
    file.path(base_plot, "RCV_CONTROL", "sem_locs", "Volcano",
              paste0(tissue, "_volcano_RHDV_vs_CONTROL.png")),
    right_label = "RCV",
    left_label = "Control",
    point_size = 3,
    lfc_min = -3,
    lfc_max = 3,
    max.overlaps = 50,
    highlight_genes = highlight_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]]
  )

  plot_volcano(
    res3,
    paste(tissue, "- RHDV vs RCV"),
    file.path(base_plot, "RHDV_RCV", "sem_locs", "Volcano",
              paste0(tissue, "_volcano_RCV_vs_RHDV.png")),
    right_label = "RCV",
    left_label = "RHDV",
    point_size = 3,
    lfc_min = -3,
    lfc_max = 3,
    max.overlaps = 50,
    highlight_genes = highlight_lookup[[paste0(tissue, ".RCV_vs_RHDV")]],
    fig_num = fig_lookup[[paste0(tissue, ".RCV_vs_RHDV")]],
    fig_label = paste0(tissue, "_RCV_vs_RHDV")
  )
}




##Condição especial para o RCV vs Control

dirs <- c(
  "RCV_CONTROL/sem_locs",
  "RCV_CONTROL/sem_locs/Condition"
)

for (d in dirs) {
  dir.create(file.path(base_plot, d), recursive = TRUE, showWarnings = FALSE)
}


plot_volcano <- function(res, title, path_plot,
                         right_label, left_label,
                         remove_loc = TRUE,
                         lfc_min = -3, lfc_max = 3,
                         padj_min = 0.05,
                         point_size = 2.5,
                         max.overlaps = 20,
                         label_all_in_window = FALSE,
                         x_zoom = NULL,
                         y_zoom = NULL,
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
  df$Highlight <- df$Gene %in% highlight_genes
  df_highlight <- df[df$Highlight, ]

  # rotular APENAS os genes destacados no texto; se nao houver nenhum,
  # rotular só os top_n_label mais significativos - nunca todos os pontos
  # da janela (label_all_in_window deixou de controlar isto, mantido so
  # por compatibilidade com as chamadas existentes)
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

  if (!is.null(x_zoom) || !is.null(y_zoom)) {
    p <- p + coord_cartesian(xlim = x_zoom, ylim = y_zoom)
  }

  ggsave(path_plot, plot = p, width = 11, height = 8.5, dpi = 320, bg = "white")

  if (!is.null(fig_num)) {
    ggsave(
      file.path(figuras_tese_dir, paste0(fig_num, "_", fig_label, ".png")),
      plot = p, width = 11, height = 8.5, dpi = 320, bg = "white"
    )
  }
}

for (tissue in names(results_list)) {
  
  res2 <- results_list[[tissue]]$RCV_vs_CONTROL
  
  # RCV vs CONTROL com regras especiais por tecido
  if (tissue == "Duodenum") {
    # NOTA: x_zoom/y_zoom removidos - a lista de genes destacados atual
    # (highlight_lookup[["Duodenum.RCV_vs_CONTROL"]]) tem log2FC entre
    # ~1.7 e ~9.3, e o zoom antigo (x_zoom = c(-2, 2)) cortava a maioria
    # dos pontos fora da area visivel. Como coord_cartesian() so recorta
    # a "janela" sem remover os pontos dos dados, o geom_text_repel
    # empurrava os rotulos desses genes cortados para a margem direita,
    # numa coluna vertical sem ligacao ao ponto real. Sem zoom, o grafico
    # ajusta-se automaticamente ao intervalo real dos dados, como nos
    # restantes tecidos.
    plot_volcano(
      res2,
      paste(tissue, "- RCV vs Control"),
      file.path(base_plot, "RCV_CONTROL", "sem_locs", "Volcano",
                paste0(tissue, "_volcano_RCV_vs_CONTROL.png")),
      right_label = "RCV",
      left_label = "Control",
      remove_loc = TRUE,
      point_size = 3,
      lfc_min = -3,
      lfc_max = 3,
      padj_min = 0.05,
      max.overlaps = 100,
      label_all_in_window = FALSE,
      highlight_genes = highlight_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_num = fig_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_label = paste0(tissue, "_RCV_vs_Control")
    )
  } else if (tissue %in% c("Liver", "Thymus", "Spleen")) {
    plot_volcano(
      res2,
      paste(tissue, "- RCV vs Control"),
      file.path(base_plot, "RCV_CONTROL", "sem_locs", "Volcano",
                paste0(tissue, "_volcano_RCV_vs_CONTROL.png")),
      right_label = "RCV",
      left_label = "Control",
      remove_loc = TRUE,
      point_size = 3,
      lfc_min = -3,
      lfc_max = 3,
      max.overlaps = Inf,
      label_all_in_window = TRUE,
      highlight_genes = highlight_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_num = fig_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_label = paste0(tissue, "_RCV_vs_Control")
    )
  } else {
    plot_volcano(
      res2,
      paste(tissue, "- RCV vs Control"),
      file.path(base_plot, "RCV_CONTROL", "sem_locs", "Volcano",
                paste0(tissue, "_volcano_RCV_vs_CONTROL.png")),
      right_label = "RCV",
      left_label = "Control",
      remove_loc = TRUE,
      point_size = 3,
      lfc_min = -3,
      lfc_max = 3,
      max.overlaps = 50,
      highlight_genes = highlight_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_num = fig_lookup[[paste0(tissue, ".RCV_vs_CONTROL")]],
      fig_label = paste0(tissue, "_RCV_vs_Control")
    )
  }
  
  
}


