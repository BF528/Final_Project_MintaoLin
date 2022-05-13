suppressPackageStartupMessages(library(tidyverse))
setwd("/projectnb/bf528/students/minty/final/count_files/")
library(DESeq2)

#read in sample metadata, sample matrix, and counts matrix from problem 3
rna_info <- read_csv("/project/bf528/project_3/groups/group_5_rna_info.csv")
sample_matrix <- read_csv("/project/bf528/project_3/samples/control_counts.csv")
counts_df <- read_csv("counts_matrix.csv")

#extract control sample names
control_samples <- rna_info %>% 
  filter(mode_of_action == "Control") %>%
  pull(Run)

#keep control columns from sample matrix
sample_matrix <- sample_matrix %>% select(Geneid, control_samples)

#join counts_df and subsetted sample_matrix by gene id - this creates 4.1 matrix
combined_counts <- full_join(counts_df, sample_matrix, by = "Geneid") %>% 
  column_to_rownames("Geneid") %>%
  rowwise() %>%
  filter(!(sum(SRR1177963:SRR1178048) == 0))

inf_1 <- rna_info %>% 
  filter(mode_of_action %in% c("DNA_Damage", "Control") & vehicle == "CMC_.5_%")

inf_2 <- rna_info %>% 
  filter(mode_of_action %in% c("CAR/PXR", "Control") & vehicle == "CORN_OIL_100_%")

inf_3 <- rna_info %>% 
  filter(mode_of_action %in% c("PPARA", "Control") & vehicle == "CMC_.5_%")

run_deseq <- function(count_df, col_data, name, norm_name) {
  runs <- col_data %>% 
    pull(Run)
  
  counts <- count_df %>% 
    select(runs)
  
  # create the DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design= ~ mode_of_action
  )
  
  # relevel mode_of_action as factor
  dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')
  
  dds <- DESeq(dds)
  res <- results(dds)
  res <- lfcShrink(dds, coef=2)
  
  write.csv(res,name)
  write.csv(counts(dds,normalized=TRUE),norm_name) #4.7 normalized matrix
  
  return (res)
}

res1 <- run_deseq(combined_counts, inf_1, "AFLATOXIN_B1.csv", "norm_counts_AFLATOXIN_B1.csv")

res2 <- run_deseq(combined_counts, inf_2, "MICONAZOLE.csv", "norm_counts_MICONAZOLE.csv")

res3 <- run_deseq(combined_counts, inf_3, "PIRINIXIC_ACID.csv", "norm_counts_PIRINIXIC_ACID.csv")

#4.4 sort by p value and write to files 
deseq_res1 <- read_csv("AFLATOXIN_B1.csv") %>% 
  arrange(padj) %>%
  rename("...1" = "Geneid")
write_csv(deseq_res1, "AFLATOXIN_B1.csv")
deseq_res2 <- read_csv("MICONAZOLE.csv") %>% 
  arrange(padj) %>%
  rename("...1" = "Geneid")
write_csv(deseq_res2, "MICONAZOLE.csv")
deseq_res3 <- read_csv("PIRINIXIC_ACID.csv") %>% 
  arrange(padj) %>%
  rename("...1" = "Geneid")
write_csv(deseq_res3, "PIRINIXIC_ACID.csv")

#number of significant genes at 0.05 level: 10,839 genes total, 
#1242 from AFLATOXIN_B1 (DNA_Damage MOA), 
#6513 from MICONAZOLE (CAR/PXR MOA),
#3152 from PIRINIXIC_ACID (PPARA MOA)
sig_genes <- sum(deseq_res1$padj < 0.05, na.rm = TRUE) + 
  sum(deseq_res2$padj < 0.05, na.rm = TRUE) + 
  sum(deseq_res3$padj < 0.05, na.rm = TRUE)

#4.5 top 10 DE genes for each 
write_csv(deseq_res1[1:10,], "../table_data/DNA_Damage_table.csv")
write_csv(deseq_res2[1:10,], "../table_data/CAR-PXR_table.csv")
write_csv(deseq_res3[1:10,], "../table_data/PPARA_table.csv")


norm_counts1 <- read_csv("norm_counts_AFLATOXIN_B1.csv")
norm_counts2 <- read_csv("norm_counts_MICONAZOLE.csv")
norm_counts3 <- read_csv("norm_counts_PIRINIXIC_ACID.csv")


## 4.6
#Create histograms of fold change values from the significant DE genes 
#for each analysis. Also create scatter plots of fold change vs nominal p-value.
graph_df <- deseq_res1 %>% mutate(Group = "DNA_Damage") %>% 
  bind_rows(deseq_res2 %>% mutate(Group = "CAR/PXR"), 
            deseq_res3 %>% mutate(Group = "PPARA")) %>%
  filter(padj < 0.05)


graph_df %>%
  ggplot() + 
  geom_histogram(aes(x = log2FoldChange, fill = Group), bins = 50) + 
  facet_wrap(~ Group) + 
  labs(title = "Log2 Fold Change of Significant DE Genes", 
       x = "Log2 Fold Change", y = "Count") + 
  theme(legend.position="none")

ggsave("../plots/DEhist.png")

graph_df %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = padj, color = Group), alpha = 0.1) + 
  facet_wrap(~ Group) + 
  labs(title = "Log2 Fold Change vs. Adjusted p-value", x = "Log2 Fold Change", 
       y = "Adjusted p-value") + 
  theme(legend.position="none")

ggsave("../plots/log2foldchange.png")  
