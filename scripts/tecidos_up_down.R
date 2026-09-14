library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Se o objeto importado se chama tecidos_per_grupo e ele é a lista:

### counts genes by group 
count_deg <- function(res_obj, alpha = 0.05, lfc = 1){
  df <- as.data.frame(res_obj)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  sig <- df[df$padj < alpha & abs(df$log2FoldChange) > lfc, ]
  
  c(
    tested = nrow(df),
    sig_total = nrow(sig),
    up = sum(sig$log2FoldChange >  lfc),
    down = sum(sig$log2FoldChange < -lfc)
  )
}

deg_by_tissue <- function(contrast_name, alpha = 0.05, lfc = 1){
  out <- t(sapply(results_list, function(x){
    obj <- x[[contrast_name]]
    if (is.null(obj)) return(c(tested = NA, sig_total = NA, up = NA, down = NA))
    count_deg(obj, alpha = alpha, lfc = lfc)
  }))
  as.data.frame(out)
}

# Ordem desejada
tissue_order <- c("Duodenum", "Liver", "Spleen", "Thymus")
contrast_order <- c("RHDV_vs_CONTROL", "RCV_vs_CONTROL", "RCV_vs_RHDV")
contrast_labels <- c(
  RHDV_vs_CONTROL = "RHDV vs CONTROL",
  RCV_vs_CONTROL  = "RCV vs CONTROL",
  RCV_vs_RHDV     = "RCV vs RHDV"
)

### 1) Tabela resumo total
summary_total <- bind_rows(
  deg_by_tissue("RHDV_vs_CONTROL") %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RHDV_vs_CONTROL"),
  deg_by_tissue("RCV_vs_CONTROL")  %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RCV_vs_CONTROL"),
  deg_by_tissue("RCV_vs_RHDV")     %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RCV_vs_RHDV")
) %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Contrast = factor(Contrast, levels = contrast_order, labels = contrast_labels[contrast_order]),
    fill_value = log10(sig_total + 1)
  )

### Plot 1: total DEGs por tecido
p_total <- ggplot(summary_total, aes(x = Contrast, y = Tissue, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sig_total), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "grey95", high = "midnightblue",
                      name = expression(log[10]~"(Count+1)")) +
  labs(
    title = "Differentially Expressed Genes by Tissue",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 16),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p_total)

### 2) Tabela resumo up/down
summary_reg <- bind_rows(
  deg_by_tissue("RHDV_vs_CONTROL") %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RHDV_vs_CONTROL"),
  deg_by_tissue("RCV_vs_CONTROL")  %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RCV_vs_CONTROL"),
  deg_by_tissue("RCV_vs_RHDV")     %>% rownames_to_column("Tissue") %>% mutate(Contrast = "RCV_vs_RHDV")
) %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Contrast = factor(Contrast, levels = contrast_order, labels = contrast_labels[contrast_order])
  ) %>%
  select(Tissue, Contrast, up, down) %>%
  pivot_longer(cols = c(up, down), names_to = "Regulation", values_to = "Count") %>%
  mutate(
    Regulation = recode(Regulation,
                        up = "Upregulated",
                        down = "Downregulated"),
    Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated")),
    fill_value = log10(Count + 1)
  )

### Plot 2: Upregulated / Downregulated
p_reg <- ggplot(summary_reg, aes(x = Contrast, y = Tissue, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "grey95", high = "midnightblue",
                      name = expression(log[10]~"(Count+1)")) +
  facet_wrap(~ Regulation, ncol = 2) +
  labs(
    title = "Differentially Expressed Genes (DEGs) by Tissue",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 16),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = "grey95"),
    strip.text = element_text(size = 16)
  )

print(p_reg)

ggsave("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/tecidos_total.png",
       p_total, width = 12, height = 6, dpi = 300, bg = "white")

ggsave("/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/tecidos_up_down.png",
       p_reg, width = 14, height = 7, dpi = 300, bg = "white")

### 3) Tecidos-alvo vs nao-alvo (target vs non-target tissues) -----------
# Mesmo estilo de heatmap dos plots acima, mas usado especificamente para
# a seccao "target tissue vs non-target tissues" do texto: o tecido-alvo
# de cada comparacao fica marcado com um * junto ao numero.
#
# Os valores vem diretamente do DEG_summary_tecidos_sem_locs.csv gerado
# pelo tecidos_per_grupo.R - a mesma fonte oficial usada na Tabela 1 - e
# nao sao recalculados aqui, para garantir que os numeros batem sempre
# certo com a tabela e com o texto (count_deg() acima, usada nos plots 1
# e 2, nao remove os LOC e por isso da valores diferentes dos oficiais).

target_lookup <- list(
  "RHDV vs CONTROL" = c("Liver"),
  "RCV vs CONTROL"  = c("Duodenum"),
  "RCV vs RHDV"     = c("Duodenum", "Liver")
)

summary_target <- read.csv(
  "/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/DEG_summary_tecidos_sem_locs.csv",
  stringsAsFactors = FALSE
)

summary_target <- summary_target %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Comparison = factor(Comparison, levels = c("RHDV vs CONTROL", "RCV vs CONTROL", "RCV vs RHDV")),
    is_target = mapply(function(comp, tis) tis %in% target_lookup[[as.character(comp)]],
                        Comparison, Tissue),
    fill_value = log10(Total + 1),
    Label = ifelse(is_target, paste0(Total, "*"), as.character(Total))
  )

p_target <- ggplot(summary_target, aes(x = Comparison, y = Tissue, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Label), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "grey95", high = "midnightblue",
                      name = expression(log[10]~"(Total+1)")) +
  labs(
    title = "Total DEGs by Tissue and Comparison",
    subtitle = "* = target tissue for that comparison (based on known tropism)",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 13, face = "italic"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p_target)

ggsave(
  "/Users/mariliaandrade/Desktop/Tese/resultados_final/tecidos_grupo/tecidos_target_vs_nontarget.png",
  p_target, width = 9, height = 6, dpi = 300, bg = "white"
)

# copia para a pasta figuras_tese - Figure09, logo a seguir ao heatmap
# global (Figure08), na seccao "Overview of Differentially Expressed
# Genes" do texto (antes de entrar comparacao a comparacao)
figuras_tese_dir <- "/Users/mariliaandrade/Desktop/Tese/resultados_final2/figuras_tese"
dir.create(figuras_tese_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(figuras_tese_dir, "Figure09_Target_vs_NonTarget_Tissues.png"),
  p_target, width = 9, height = 6, dpi = 300, bg = "white"
)
