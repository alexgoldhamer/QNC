library(ggplot2)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(DESeq2)
library(tidyverse)
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(data.table)
library(pheatmap)
library(stringr)
library(gridExtra)
library(EnhancedVolcano)
library(ggpubr)
library(edgeR)
library(readxl)

# Set working directory
setwd("~/Zhou_Lab/YY1_Stress_RNA_seq_10_16_25")

###############################################################
# Differential Expression Analysis (DESeq2)
# Separate runs for Males and Females
# Genotype: WT vs KO
# Stress: CTR vs aCUS
###############################################################

# Load data
countData <- read.csv("FeatureCounts_M39_Stress.txt", sep = "\t", header = TRUE)
metaData <- read.csv("Metadata_RNA_Seq_YY1_Stress_10_16_25.txt", sep = "\t", header = TRUE)

###############################################################
# ---- Male Analysis ----
###############################################################

# Subset metadata and count matrix for males
meta_male <- subset(metaData, Sex == "M")
count_male <- countData[, c("Geneid", paste0("X", meta_male$Sample))]

# Set gene IDs as rownames
rownames(count_male) <- count_male$Geneid
count_male <- count_male[, -1]

# Create DESeq2 object
dds_male <- DESeqDataSetFromMatrix(
  countData = count_male,
  colData = meta_male,
  design = ~ Stress + Genotype
)

# Set reference levels
dds_male$Stress <- relevel(dds_male$Stress, ref = "CTR")
dds_male$Genotype <- relevel(dds_male$Genotype, ref = "WT")

# Run DESeq2
dds_male <- DESeq(dds_male)

# View available result names
resultsNames(dds_male)

# Extract results for stress and genotype contrasts
res_male_stress <- results(dds_male, name = "Stress_aCUS_vs_CTR")
res_male_genotype <- results(dds_male, name = "Genotype_KO_vs_WT")

# Save normalized counts and rlog
norm_counts_male <- counts(dds_male, normalized = TRUE)
rlog_male <- rlog(dds_male, blind = TRUE)

# Write out results
write.csv(as.data.frame(res_male_stress), "DEG_Male_Stress_aCUS_vs_CTR.csv")
write.csv(as.data.frame(res_male_genotype), "DEG_Male_Genotype_KO_vs_WT.csv")

library(EnhancedVolcano)

res_male_stress <- results(dds_male, name = "Stress_aCUS_vs_CTR")
res_male_genotype <- results(dds_male, name = "Genotype_KO_vs_WT")

# ---- Volcano plots for Males ----
pdf("Volcano_Male_Stress_aCUS_vs_CTR.pdf", width = 7, height = 6)
EnhancedVolcano(res_male_stress,
                lab = rownames(res_male_stress),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Male: Stress (aCUS vs CTR)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3
)
dev.off()

pdf("Volcano_Male_Genotype_KO_vs_WT.pdf", width = 7, height = 6)
EnhancedVolcano(res_male_genotype,
                lab = rownames(res_male_genotype),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Male: Genotype (KO vs WT)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3
)
dev.off()

###############################################################
# ---- Female Analysis ----
###############################################################

# Subset metadata and count matrix for females
meta_female <- subset(metaData, Sex == "F")
count_female <- countData[, c("Geneid", paste0("X", meta_female$Sample))]

# Set gene IDs as rownames
rownames(count_female) <- count_female$Geneid
count_female <- count_female[, -1]

# Create DESeq2 object
dds_female <- DESeqDataSetFromMatrix(
  countData = count_female,
  colData = meta_female,
  design = ~ Stress + Genotype
)

# Set reference levels
dds_female$Stress <- relevel(dds_female$Stress, ref = "CTR")
dds_female$Genotype <- relevel(dds_female$Genotype, ref = "WT")

# Run DESeq2
dds_female <- DESeq(dds_female)

# View available result names
resultsNames(dds_female)

# Extract results for stress and genotype contrasts
res_female_stress <- results(dds_female, name = "Stress_aCUS_vs_CTR")
res_female_genotype <- results(dds_female, name = "Genotype_KO_vs_WT")

# Save normalized counts and rlog
norm_counts_female <- counts(dds_female, normalized = TRUE)
rlog_female <- rlog(dds_female, blind = TRUE)

# Write out results
write.csv(as.data.frame(res_female_stress), "DEG_Female_Stress_aCUS_vs_CTR.csv")
write.csv(as.data.frame(res_female_genotype), "DEG_Female_Genotype_KO_vs_WT.csv")

res_female_stress <- results(dds_female, name = "Stress_aCUS_vs_CTR")
res_female_genotype <- results(dds_female, name = "Genotype_KO_vs_WT")

# ---- Volcano plots for Females ----
pdf("Volcano_Female_Stress_aCUS_vs_CTR.pdf", width = 7, height = 6)
EnhancedVolcano(res_female_stress,
                lab = rownames(res_female_stress),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Female: Stress (aCUS vs CTR)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3
)
dev.off()

pdf("Volcano_Female_Genotype_KO_vs_WT.pdf", width = 7, height = 6)
EnhancedVolcano(res_female_genotype,
                lab = rownames(res_female_genotype),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Female: Genotype (KO vs WT)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3
)
dev.off()

###############################################################
# ---- Optional: GSEA Ranking Files ----
###############################################################

# Create ranked vectors for GSEA (using DESeq2 test statistic)

rank_male_stress <- res_male_stress$stat
names(rank_male_stress) <- rownames(res_male_stress)
write.table(rank_male_stress, "rank_male_stress_aCUS_vs_CTR.txt", sep = "\t", quote = FALSE, col.names = NA)

rank_male_genotype <- res_male_genotype$stat
names(rank_male_genotype) <- rownames(res_male_genotype)
write.table(rank_male_genotype, "rank_male_genotype_KO_vs_WT.txt", sep = "\t", quote = FALSE, col.names = NA)

rank_female_stress <- res_female_stress$stat
names(rank_female_stress) <- rownames(res_female_stress)
write.table(rank_female_stress, "rank_female_stress_aCUS_vs_CTR.txt", sep = "\t", quote = FALSE, col.names = NA)

rank_female_genotype <- res_female_genotype$stat
names(rank_female_genotype) <- rownames(res_female_genotype)
write.table(rank_female_genotype, "rank_female_genotype_KO_vs_WT.txt", sep = "\t", quote = FALSE, col.names = NA)
