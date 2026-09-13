# Host transcriptomic response to RHDV and RCV in the European rabbit

MSc dissertation — Faculty of Sciences of the University of Porto (FCUP), 2026
Marília Andrade

RNA-seq of liver, duodenum, spleen and thymus from *Oryctolagus cuniculus* infected with a pathogenic lagovirus (RHDV, genotype GI.1c) and a non-pathogenic one (RCV, genotype GI.4), sampled at 4 and 14 days post-infection.

This repository collects the material behind the dissertation: study design, the analysis pipeline and its parameters, the differential expression results, and the supplementary tables and scripts.

---

## Contents

- [Abstract](#abstract)
- [Study design](#study-design)
- [Analysis pipeline](#analysis-pipeline)
- [Results](#results)
- [Supplementary material](#supplementary-material)
- [Repository layout](#repository-layout)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Data availability](#data-availability)

---

## Abstract

RHDV (genotype GI.1c) and RCV (genotype GI.4) are lagoviruses of the European rabbit (*Oryctolagus cuniculus*) with opposite outcomes: RHDV causes an acute, frequently fatal hepatitis, while RCV infection remains subclinical. This thesis asks whether that contrast is visible in the host transcriptome, and in which organs. RNA-seq data came from an archived challenge experiment (RHDV, n = 5; RCV, n = 5; Control, n = 3) in which liver, duodenum, spleen and thymus were sampled at 4 and 14 days post-infection. Reads were processed by a pipeline written for this project, covering trimming, HISAT2 alignment, featureCounts quantification and DESeq2 differential expression in one standardized workflow.

Principal component analysis of variance-stabilized counts showed samples grouping by tissue rather than by infection status: PC1 (65% of the variance) separated liver from spleen and thymus, PC2 (26%) set the duodenum apart, and Control, RCV and RHDV samples overlapped within each tissue. Against control animals, RHDV altered 135 genes in the global model (127 up- and 8 downregulated), with the largest per-tissue response in the spleen (83 genes) rather than in the liver (32 genes), the organ where the virus replicates. RCV altered 6 genes in the same model, concentrated in the duodenum. Genes raised after RHDV infection included interferon-stimulated effectors (*ISG15*, *RSAD2*, *IRF7*), macrophage and scavenger-receptor markers (*CD163*, *VSIG4*, *S100A9*) and chemokines (*CCL5*, *CX3CL1*). The downregulated genes carried no clear immune annotation and fell mostly in adipocyte and lipid metabolism.

In the liver, the RHDV response was present at 4 dpi and had largely resolved by 14 dpi, combining inflammatory recruitment, macrophage scavenging and altered bile-acid transport. RCV reproduced part of that signature at 4 dpi but never the bile-acid transport component. What separated the two infections was therefore not the presence of an antiviral response, which both viruses triggered, but whether that response came with a short-lived disruption of liver function; and the organ with the largest response was not the organ the virus infects first.

**Keywords:** rabbit hemorrhagic disease virus, rabbit calicivirus, *Oryctolagus cuniculus*, lagovirus, RNA-seq, differential gene expression, host response, liver, spleen

<details>
<summary><strong>Resumo (português)</strong></summary>

O RHDV (genótipo GI.1c) e o RCV (genótipo GI.4) são lagovírus do coelho-europeu (*Oryctolagus cuniculus*) com desfechos opostos: o RHDV causa uma hepatite aguda e frequentemente fatal, enquanto a infeção por RCV permanece subclínica. Esta tese pergunta se esse contraste é visível no transcriptoma do hospedeiro e em que órgãos. Os dados de RNA-seq provêm de uma experiência de infeção arquivada (RHDV, n = 5; RCV, n = 5; Controlo, n = 3), com amostragem de fígado, duodeno, baço e timo aos 4 e aos 14 dias pós-infeção. As leituras foram processadas por um pipeline desenvolvido para este trabalho, que reúne numa única sequência normalizada o trimming, o alinhamento com HISAT2, a quantificação com featureCounts e a análise de expressão diferencial com DESeq2.

A análise de componentes principais mostrou que as amostras se agrupam por tecido e não por condição de infeção. Face aos controlos, o RHDV alterou 135 genes no modelo global, com a maior resposta por tecido no baço (83 genes) e não no fígado (32 genes), o órgão onde o vírus replica; o RCV alterou 6 genes, concentrados no duodeno. No fígado, a resposta ao RHDV estava presente aos 4 dpi e tinha desaparecido aos 14 dpi, combinando recrutamento inflamatório, atividade de macrófagos e alteração do transporte de ácidos biliares. O RCV reproduziu parte desta assinatura aos 4 dpi, mas nunca a componente de transporte de ácidos biliares. O que separa as duas infeções não é, portanto, a existência de uma resposta antiviral, que ambos os vírus desencadeiam, mas o facto de essa resposta vir acompanhada de uma perturbação transitória da função hepática.

**Palavras-chave:** vírus da doença hemorrágica do coelho, calicivírus do coelho, *Oryctolagus cuniculus*, lagovírus, RNA-seq, expressão diferencial, resposta do hospedeiro, fígado, baço

</details>

---

## Study design

The tissue samples were collected in 2014 at CSIRO (Black Mountain Laboratories, Canberra, Australia), under approval ESAEC 13-11 from the CSIRO Animal Ethics Committee and in accordance with the Australian Animal Research Acts and the Australian Code of Practice for the Care and Use of Animals for Scientific Purposes. Humane endpoints were applied throughout.

| Group | n | Inoculum |
|---|---|---|
| RHDV (GI.1c) | 5 | 1 mL oral, ≈ 1.7 × 10⁷ viral copies, 500 ID<sub>50</sub> |
| RCV (GI.4) | 5 | 1 mL oral, ≈ 1.7 × 10⁷ viral copies (no ID<sub>50</sub>: the virus is non-pathogenic, so there is no clinical endpoint to calculate one from) |
| Control | 3 | — |

Animals were 4 weeks old. Liver, duodenum, spleen and thymus were harvested for RNA extraction. Two animals from the RHDV group were euthanized at 4 dpi to capture the early phase of infection; the remaining animals were sampled at 14 dpi.

Libraries were prepared with the TruSeq RNA Sample Preparation Kit v2 (Illumina): oligo-dT enrichment of polyadenylated RNA, chemical fragmentation, first- and second-strand cDNA synthesis, end repair, 3′ adenylation, ligation of indexed adapters, and PCR amplification before pooling and paired-end sequencing.

---

## Analysis pipeline

Every sample passed through the same steps with the same parameters, which is what makes the three groups directly comparable.

| Stage | Tool | Version | Key settings |
|---|---|---|---|
| Trimming and adapter removal | Trimmomatic | 0.39 | Paired-end; `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True`, `LEADING:3`, `TRAILING:3`, `SLIDINGWINDOW:4:15`, `MINLEN:36` |
| Alignment | HISAT2 | — | Reference *O. cuniculus* mOryCun1.1 (NCBI `GCF_964237555.1`); index built with `hisat2-build`, known exons and splice sites supplied; `--rna-strandness` unstranded |
| Post-alignment | SAMtools | — | SAM → BAM conversion, coordinate sorting, indexing |
| Quantification | featureCounts | 1.4.2 | Gene-level counting, aggregation by `gene_id` from the reference GTF |
| Differential expression | DESeq2 | 1.34 (R 4.1.3) | Negative binomial GLM, Wald test, Benjamini–Hochberg correction |

Overall alignment rates ran from 87.4% to 96.3% (mean 91.5%) across all libraries.

**Design formulas**

| Level | Design |
|---|---|
| Global | `~ lane + tissue + condition` |
| Per tissue (matrix subset by organ) | `~ lane + condition` |
| Time course | `~ lane + tissue + group`, where `group` combines condition and time point (`RHDV_4dpi`, `RHDV_14dpi`, `RCV_4dpi`, `RCV_14dpi`, `Control_14dpi`) |

A gene was called differentially expressed at adjusted *p* < 0.05 and |log<sub>2</sub> fold change| > 1.

---

## Results

Counts below exclude loci that the mOryCun1.1 annotation reports only under a `LOC` identifier. The versions that keep those loci are listed under [Supplementary material](#supplementary-material).

### Table 1 — DEGs per comparison and level of analysis

"Global" means the model was fitted to all tissues together with tissue as a factor.

| Comparison | Level | Total | Up | Down |
|---|---|---:|---:|---:|
| RHDV vs Control | Global | 135 | 127 | 8 |
| | Duodenum | 9 | 7 | 2 |
| | Liver *(RHDV target)* | 32 | 30 | 2 |
| | Spleen | 83 | 20 | 63 |
| | Thymus | 0 | 0 | 0 |
| RCV vs Control | Global | 6 | 4 | 2 |
| | Duodenum *(RCV target)* | 19 | 17 | 2 |
| | Liver | 0 | 0 | 0 |
| | Spleen | 16 | 2 | 14 |
| | Thymus | 2 | 1 | 1 |
| RCV vs RHDV | Global | 289 | 157 | 132 |
| | Duodenum | 691 | 637 | 54 |
| | Liver | 15 | 2 | 13 |
| | Spleen | 7 | 6 | 1 |
| | Thymus | 36 | 0 | 36 |

### Table 2 — DEGs per comparison and time point

| Comparison | Time point | Total | Up | Down |
|---|---|---:|---:|---:|
| RHDV vs Control | 4 dpi | 1363 | 906 | 457 |
| RHDV vs Control | 14 dpi | 4 | 2 | 2 |
| RCV vs Control | 4 dpi | 124 | 112 | 12 |
| RCV vs Control | 14 dpi | 5 | 1 | 4 |
| RHDV vs RCV | 4 dpi | 787 | 342 | 445 |
| RHDV vs RCV | 14 dpi | 10 | 1 | 9 |
| RHDV, 14 dpi vs 4 dpi | within RHDV | 1489 | 498 | 991 |
| RCV, 14 dpi vs 4 dpi | within RCV | 58 | 17 | 41 |

Counts are not directly comparable between tissues: detection depends on within-group dispersion and on library size, both of which differ by organ and were estimated here from a small number of animals per group. They are read as descriptive rather than as a ranking of the intensity of the response.

---

## Supplementary material

| ID | Content | Files |
|---|---|---|
| Table S1 | DESeq2 design formulas and the significance thresholds | — |
| Note S1 | Trimmomatic call | — |
| Table S2 | DEGs in the global models, uncharacterised loci retained | `RHDV_vs_CONTROL_sig.csv`, `RCV_vs_CONTROL_sig.csv`, `RCV_vs_RHDV_sig.csv` |
| Table S3 | DEGs per tissue, uncharacterised loci retained | `<Tissue>_<comparison>_sig.csv`, one per organ and comparison (12 files) |
| Table S4 | DEGs per time point, uncharacterised loci retained | `RHDV_4DPI_vs_CONTROL_14DPI_sig.csv` and the other seven time contrasts |
| Table S5 | `LOC` identifiers that have since received an approved symbol, and the comparisons in which each appears | `LOC_genes_renamed_NCBI.csv`, `LOC_genes_onde_aparecem.csv` |
| Table S6 | Complete DESeq2 output, every gene that passed the count filter | `<comparison>_all_genes.csv`, `<Tissue>_<comparison>.csv` |
| Scripts S1–S3 | `preprocessing_fasta.py`, `alignment.py`, `feature_count.py` | — |

Each result table has one row per gene and the columns gene symbol, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj`. Tables S2 to S4 are restricted to genes meeting *p*<sub>adj</sub> < 0.05 and |log<sub>2</sub>FC| > 1, ordered by *p*<sub>adj</sub>.

---

## Repository layout

```
.
├── README.md          this page
└── _config.yml        GitHub Pages configuration
```

Still to add: the pipeline scripts (`preprocessing_fasta.py`, `alignment.py`, `feature_count.py`), the R analysis scripts (`overall_analyzes.R`, `tecidos_per_grupo.R`, `with_time.R`, `build_summary_tables.R`), the supplementary CSV tables, and the figures.

---

## Reproducing the analysis

1. **Pre-processing** — `preprocessing_fasta.py` runs Trimmomatic over each paired-end library and writes the trimmed and unpaired FASTQ files.
2. **Alignment** — `alignment.py` builds the HISAT2 index from mOryCun1.1, maps the trimmed pairs, and converts, sorts and indexes the BAM files.
3. **Quantification** — `feature_count.py` validates the `gene_id` attribute across the GTF, runs featureCounts and assembles the gene-level count matrix.
4. **Differential expression** — `overall_analyzes.R` (global models), `tecidos_per_grupo.R` (per-tissue models) and `with_time.R` (time course) take the count matrix and the sample table and write the result tables, both with and without the uncharacterised loci.
5. **Summary tables** — `build_summary_tables.R` assembles Tables 1 and 2 from the outputs of step 4, so the numbers in the dissertation are not written by hand.

---

## Data availability

The RNA-seq data come from an archived challenge experiment carried out at CSIRO in 2014. Requests for the raw sequencing data should be directed through the dissertation's supervision.
