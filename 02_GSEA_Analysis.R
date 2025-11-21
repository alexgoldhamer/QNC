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

####Subset the data for males
metaData = subset(metaData[metaData$Sex == "M", ])
countData = countData[, c("Geneid", paste0("X", metaData$Sample))]

####Subset the data for females
metaData = subset(metaData[metaData$Sex == "F", ])
countData = countData[, c("Geneid", paste0("X", metaData$Sample))]

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design = ~ Stress + Genotype, tidy = TRUE)
dds$Sex <- relevel(dds$Sex, ref = "M")
dds$Stress <- relevel(dds$Stress, ref = "CTR")
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

colData(dds)
dds <- DESeq(dds)
resultsNames(dds)
normalized_counts <- counts(dds, normalized = TRUE)
rlog_expression <- rlog(dds, blind = TRUE)

#GSEA
#Gene Set Enrichment Analysis 
#run to see naming options
resultsNames(dds)

#DEG analysis looking genotype comparison
res <- results(dds, name="Genotype_KO_vs_WT", tidy = TRUE)

res2 <- res %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(row) %>%
  summarize(stat=mean(stat))
res2
ranks <- deframe(res2)

pathways.kegg <- msigdbr(species = "mouse", category = "C2", subcategory = "KEGG_LEGACY") #CHECK NAME
pathways.reactome <- msigdbr(species = "mouse", category = "C2", subcategory = "REACTOME")
pathways.hallmark <- msigdbr(species = "mouse", category = "H")
pathways.tf <- msigdbr(species = "mouse", category = "C3", subcategory = "TFT:GTRD")
pathway_list_kegg = split(x = pathways.kegg$gene_symbol, f = pathways.kegg$gs_name)
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list_reactome = split(x = pathways.reactome$gene_symbol, f = pathways.reactome$gs_name)
pathway_list_tf <- split(x = pathways.tf$gene_symbol, f = pathways.tf$gs_name)

#Set the pathway list variable (choose one from here)
pathway_list <- pathway_list_reactome
pathway_list <- c(pathway_list_kegg, pathway_list_hallmark) #Start with this, combination of KEGG and Hallmark pathways 
pathway_list <- pathway_list_tf

GSEA_result <-  fgsea(pathways = pathway_list, stats = ranks)

sig_pathways <- GSEA_result[GSEA_result$padj<0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways = na.omit(sig_pathways)

top_50_pathways <- head(sig_pathways, 25)
bottom_50_pathways <- tail(sig_pathways, 25)
sig_pathways_filter <- rbind(top_50_pathways, bottom_50_pathways)
sig_pathways_filter <- sig_pathways_filter[!duplicated(sig_pathways_filter[["pathway"]]), ]
sig_pathways_filter$pathway <- factor(
  sig_pathways_filter$pathway,
  levels = sig_pathways_filter$pathway[order(sig_pathways_filter$NES)]
)

p3 <- ggplot(sig_pathways_filter, aes(x = pathway, y = NES, fill = padj)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient(low = "blue",high = "yellow") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR (p-adjusted)") +  # Update axis and legend labels
  coord_flip()+ ggtitle("") + 
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, family = "sans", size = 8, color = "black"),
    axis.text.y = element_text(family = "sans", size = 12, color = "black"),
    axis.line = element_blank(), 
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt")
  ) + theme_bw()
p3

