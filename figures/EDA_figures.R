library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


cts <- read.csv("input/miRNA_counts_matrix.csv")
sample_df <- read.csv("input/sample_info.csv")

names(cts)[1] <- "miRNA"

miRNA_counts_all <- cts |>
  rename_with(~ str_remove(., "X")) |> #remove prefix "X" from all column names
  column_to_rownames(var = names(cts)[1]) #convert column of miRNAs to row names

miRNA_counts_df <- cts |>
  rename_with(~ str_remove(., "X"))
# Formating Column Data
sample_df$condition <- factor(sample_df$condition)
sample_data <- sample_df |>
  column_to_rownames(var = names(sample_df)[1])

sample_counts_mat <- t(as.matrix(miRNA_counts_all))
sample_counts_df <- rownames_to_column(
  as.data.frame(sample_counts_mat),
  var = "sample"
)


sample_disease_df <- sample_counts_df |>
  inner_join(sample_df, by = "sample") |>
  mutate(total_mirna_counts = rowSums(across(2:1977))) |>
  select(c(sample, condition, total_mirna_counts))

#Figure of miRNA counts by Disease
ggplot(
  sample_disease_df,
  aes(total_mirna_counts, condition, fill = condition)
) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  labs(title = "miRNA Counts by Disease", x = "Total miRNA Counts")

# Figure of Sample Counts
sample_disease_df |>
  dplyr::count(condition) |>
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "# of Samples")
