library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(apeglm)

cts <- read.csv("input/miRNA_counts_matrix.csv")
sample_df <- read.csv("input/sample_info.csv")


names(cts)[1] <- "miRNA"

miRNA_counts_all <- cts |>
  rename_with(~ str_remove(., "X")) |> #remove prefix "X" from all column names
  column_to_rownames(var = names(cts)[1]) #convert column of miRNAs to row names

miRNA_all_counts_mat <- as.matrix(miRNA_counts_all)


# Formating Column Data
sample_df$condition <- factor(sample_df$condition)
sample_data <- sample_df |>
  column_to_rownames(var = names(sample_df)[1])
table(sample_df$condition)

#Check if column names and row names line up
all(rownames(sample_data) == colnames(miRNA_all_counts_mat))


############## PD vs Healthy ##################################
PD_v_Ctrl_df <- sample_df |>
  filter(condition %in% c("Healthy Control", "Parkinson's Disease")) |>
  column_to_rownames(var = names(sample_df)[1])

PD_v_Ctrl_intersect <- intersect(
  rownames(PD_v_Ctrl_df),
  colnames(miRNA_all_counts_mat)
)

miRNA_PD_v_Ctrl_counts_mat <- miRNA_all_counts_mat[, PD_v_Ctrl_intersect]

#Verify that column names and row names line up
#Check if column names and row names line up
all(rownames(PD_v_Ctrl_df) == colnames(miRNA_PD_v_Ctrl_counts_mat))


#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = miRNA_PD_v_Ctrl_counts_mat,
  colData = PD_v_Ctrl_df,
  design = ~condition
)

#Perform pre-filtering
keep <- rowSums(counts(dds) >= 10)
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds)

res_df <- as.data.frame(res)
#with shrinkage
#resLFC <- lfcShrink(dds, coef="condition_Parkinson.s.Disease_vs_Healthy.Control", type="apeglm")

res05_0 <- res_df |>
  filter(padj < 0.5, abs(log2FoldChange) > 0)

write.csv(res05_0, "results_PD_v_healthy_ctrl.csv")
############## PD w Dementia vs Healthy ##################################
PDD_v_Ctrl_df <- sample_df |>
  filter(
    condition %in%
      c("Healthy Control", "Dementia associated with Parkinson's Disease")
  ) |>
  column_to_rownames(var = names(sample_df)[1])

PDD_v_Ctrl_intersect <- intersect(
  rownames(PDD_v_Ctrl_df),
  colnames(miRNA_all_counts_mat)
)

miRNA_PDD_v_Ctrl_counts_mat <- miRNA_all_counts_mat[, PDD_v_Ctrl_intersect]

#Verify that column names and row names line up
#Check if column names and row names line up
all(rownames(PDD_v_Ctrl_df) == colnames(miRNA_PDD_v_Ctrl_counts_mat))


#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = miRNA_PDD_v_Ctrl_counts_mat,
  colData = PDD_v_Ctrl_df,
  design = ~condition
)

#Perform pre-filtering
keep <- rowSums(counts(dds) >= 10)
dds <- dds[keep, ]

dds <- DESeq(dds)

res_PDD <- results(dds)

res_PDD_df <- as.data.frame(res_PDD)

res05_0_PDD <- res_PDD_df |>
  filter(padj < 0.5, abs(log2FoldChange) > 0)

write.csv(res05_0_PDD, "results_PD_w_dementia_v_healthy_ctrl.csv")
############## Lewy Body Dementia vs Healthy ##################################
LBD_v_Ctrl_df <- sample_df |>
  filter(condition %in% c("Healthy Control", "Lewy Body Dementia")) |>
  column_to_rownames(var = names(sample_df)[1])

LBD_v_Ctrl_intersect <- intersect(
  rownames(LBD_v_Ctrl_df),
  colnames(miRNA_all_counts_mat)
)

miRNA_LBD_v_Ctrl_counts_mat <- miRNA_all_counts_mat[, LBD_v_Ctrl_intersect]

#Verify that column names and row names line up
#Check if column names and row names line up
all(rownames(LBD_v_Ctrl_df) == colnames(miRNA_LBD_v_Ctrl_counts_mat))


#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = miRNA_LBD_v_Ctrl_counts_mat,
  colData = LBD_v_Ctrl_df,
  design = ~condition
)

#Perform pre-filtering
keep <- rowSums(counts(dds) >= 10)
dds <- dds[keep, ]

dds <- DESeq(dds)

res_LBD <- results(dds)

res_LBD_df <- as.data.frame(res_LBD)

res05_0_LBD <- res_LBD_df |>
  filter(padj < 0.5, abs(log2FoldChange) > 0)

write.csv(res05_0_LBD, "results_lewy_body_dementia_v_healthy_ctrl.csv")
