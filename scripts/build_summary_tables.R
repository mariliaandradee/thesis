# Consolida as tabelas-resumo (Tabela 1 e Tabela 2 da tese) diretamente a
# partir dos ficheiros gerados pelos 3 scripts de analise (sem_locs):
#   - geral/overall_analyzes.R          -> geral/DEG_counts_sem_locs.csv
#   - tecidos_grupo/tecidos_per_grupo.R -> tecidos_grupo/DEG_summary_tecidos_sem_locs.csv
#   - tempo/with_time.R                 -> tempo/*/sem_locs/*_DEG_counts.csv
#
# Corre este script DEPOIS de correr os tres scripts acima (produz
# Table1_DEG_summary.csv e Table2_DEG_summary_tempo.csv, usados diretamente
# no capitulo de Resultados, sem numeros escritos a mao).

base_final <- "/Users/mariliaandrade/Desktop/Tese/resultados_final"
out_dir    <- "/Users/mariliaandrade/Desktop/Tese/resultados_final2"

# ---------------------------------------------------------------
# Tabela 1: Global + por tecido (RHDV vs CONTROL, RCV vs CONTROL, RCV vs RHDV)
# ---------------------------------------------------------------

geral <- read.csv(file.path(base_final, "geral", "DEG_counts_sem_locs.csv"),
                   stringsAsFactors = FALSE)
geral$Level <- "Global"
geral <- geral[, c("Comparison", "Level", "Total", "Up", "Down")]

tecidos <- read.csv(file.path(base_final, "tecidos_grupo", "DEG_summary_tecidos_sem_locs.csv"),
                     stringsAsFactors = FALSE)
tecidos$Level <- tecidos$Tissue
tecidos <- tecidos[, c("Comparison", "Level", "Total", "Up", "Down")]

table1 <- rbind(geral, tecidos)

# ordenar: RHDV vs CONTROL, RCV vs CONTROL, RCV vs RHDV; dentro de cada,
# Global primeiro e depois tecidos por ordem alfabetica
comparison_order <- c("RHDV vs CONTROL", "RCV vs CONTROL", "RCV vs RHDV")
level_order <- c("Global", "Duodenum", "Liver", "Spleen", "Thymus")
table1$Comparison <- factor(table1$Comparison, levels = comparison_order)
table1$Level <- factor(table1$Level, levels = level_order)
table1 <- table1[order(table1$Comparison, table1$Level), ]

write.csv(table1, file.path(out_dir, "Table1_DEG_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# Tabela 2: Time-course (sem_locs)
# ---------------------------------------------------------------

tempo_files <- c(
  "RHDV_CONTROL_4DPI/sem_locs/RHDV_4DPI_vs_CONTROL_14DPI_DEG_counts.csv",
  "RHDV_CONTROL_14DPI/sem_locs/RHDV_14DPI_vs_CONTROL_14DPI_DEG_counts.csv",
  "RCV_CONTROL_4DPI/sem_locs/RCV_4DPI_vs_CONTROL_14DPI_DEG_counts.csv",
  "RCV_CONTROL_14DPI/sem_locs/RCV_14DPI_vs_CONTROL_14DPI_DEG_counts.csv",
  "RHDV_RCV/sem_locs/RHDV_vs_RCV_4DPI_DEG_counts.csv",
  "RHDV_RCV/sem_locs/RHDV_vs_RCV_14DPI_DEG_counts.csv",
  "RHDV_14dpi/sem_locs/RHDV_14DPI_vs_4DPI_DEG_counts.csv",
  "RCV_4dpi/sem_locs/RCV_14DPI_vs_4DPI_DEG_counts.csv"
)

# rotulos legiveis para a Tabela 2, na mesma ordem de tempo_files acima
labels <- data.frame(
  Comparison = c("RHDV vs Control", "RHDV vs Control",
                 "RCV vs Control", "RCV vs Control",
                 "RHDV vs RCV", "RHDV vs RCV",
                 "RHDV 14 dpi vs 4 dpi", "RCV 14 dpi vs 4 dpi"),
  Timepoint  = c("4 dpi", "14 dpi",
                 "4 dpi", "14 dpi",
                 "4 dpi", "14 dpi",
                 "within RHDV", "within RCV"),
  stringsAsFactors = FALSE
)

table2_rows <- lapply(seq_along(tempo_files), function(i) {
  f <- file.path(base_final, "tempo", tempo_files[i])
  d <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(
    Comparison = labels$Comparison[i],
    Timepoint  = labels$Timepoint[i],
    Total = d$Up[1] + d$Down[1],
    Up = d$Up[1],
    Down = d$Down[1]
  )
})

table2 <- do.call(rbind, table2_rows)
write.csv(table2, file.path(out_dir, "Table2_DEG_summary_tempo.csv"), row.names = FALSE)

cat("Tabela 1 (", nrow(table1), "linhas ) e Tabela 2 (", nrow(table2),
    "linhas ) escritas em", out_dir, "\n")
