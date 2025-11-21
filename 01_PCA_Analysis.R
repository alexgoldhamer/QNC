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

setwd("~/Zhou_Lab/YY1_Stress_RNA_seq_10_16_25")

countData <- read.csv("FeatureCounts_M39_Stress.txt", sep="\t", header = TRUE)
metaData <- read.csv("Metadata_RNA_Seq_YY1_Stress_10_16_25.txt", sep="\t", header = TRUE)
dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design = ~ Sex + Stress + Genotype, tidy = TRUE)
dds$Sex <- relevel(dds$Sex, ref = "M")
dds$Stress <- relevel(dds$Stress, ref = "aCUS")
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

colData(dds)
dds <- DESeq(dds)
resultsNames(dds)
normalized_counts <- counts(dds, normalized = TRUE)
rlog_expression <- rlog(dds, blind = TRUE)

#PCA
plotPCA(rlog_expression, intgroup = "Sex", pcsToUse = c(1,2), ntop = 300000) + theme_bw()
plotPCA(rlog_expression, intgroup = "Stress", pcsToUse = c(1,2), ntop = 300000) + theme_bw()
plotPCA(rlog_expression, intgroup = "Genotype", pcsToUse = c(1,2), ntop = 300000) + theme_bw()
plotPCA(rlog_expression, intgroup = c("Sex", "Stress", "Genotype"), pcsToUse = c(1,2), ntop = 300000) + theme_bw()

pcaData = plotPCA(rlog_expression, intgroup = "Sex", pcsToUse = c(1,3), ntop = 300000, returnData = T) 
ggplot(pcaData, aes(x = Sex, y = PC3, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  labs(
    title = "PC3 by Sex",
    y = "PC3 value",
    x = "Genotype"
  )

#Identiy significant principal components 
mat <- assay(rlog_expression)  # or vst(dds) if you used vst instead
meta <- as.data.frame(colData(rlog_expression))

# Perform PCA manually on transposed expression matrix
pca <- prcomp(t(mat))

# Create a dataframe of PC coordinates + metadata
pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, meta)

# Loop through all PCs, test for sex differences
pc_names <- colnames(pca$x)
pc_results <- data.frame(PC = pc_names, p_value = NA)

for (i in seq_along(pc_names)) {
  pc <- pc_names[i]
  test <- wilcox.test(pca_df[[pc]] ~ pca_df$Genotype)
  pc_results$p_value[i] <- test$p.value
}

# Adjust for multiple testing (optional but recommended)
pc_results$padj <- p.adjust(pc_results$p_value, method = "BH")

# Show significant PCs
significant_pcs <- pc_results %>% filter(padj < 0.1)
significant_pcs

