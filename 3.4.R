suppressPackageStartupMessages(library(tidyverse))
setwd("/projectnb/bf528/students/minty/final/count_files/")

#read files
counts_df <- read_tsv("SRR1177963Aligned.txt", skip = 1) %>%
  pivot_longer(cols = contains("projectnb"), names_to = "BAM Path", values_to = "counts")

#other files to be added to the dataframe counts_df
files <- list.files(getwd(), pattern = "SRR117[0-9]{3}[0-2, 4-9]{1}Aligned.txt$")

#combine dataframes generated from each file
for (file in files) {
  counts_df <- bind_rows(counts_df, read_tsv(file, skip = 1) %>%
                           pivot_longer(cols = contains("projectnb"), 
                                        names_to = "BAM Path", 
                                        values_to = "counts"))
}

#get gene id, sample name, and counts
counts_df <- counts_df %>% 
  rowwise() %>%
  mutate(sample = str_extract(`BAM Path`, "SRR117[0-9]{4}")) %>%
  select(Geneid, sample, counts)

#color by sample type - join mode of action info from rna_info file
rna_info <- read_csv("/project/bf528/project_3/groups/group_5_rna_info.csv")

#create boxplot: 3.5 (for samples showing the distribution of counts.)
counts_df %>%
  inner_join(rna_info, by = c("sample" = "Run")) %>%
  ggplot() + 
  geom_boxplot(aes(x = log(counts), y = sample, color = mode_of_action)) + 
  labs(x = "log(Counts)", y = "Sample", title = "Distribution of Sample Counts",
       color = "Mode of Action")

ggsave("counts_boxplot.png")

counts_df <- counts_df %>%
  pivot_wider(names_from = sample, values_from = counts)

#create csv file
write_csv(counts_df, "counts_matrix.csv")

