library(tidyverse)

library(ggvenn)

tpm_data <- read.csv("tpm.csv", row.names = 1)

group_data <- read.csv("groupsize.csv")

colnames(tpm_data) <- gsub("\\.", "_", colnames(tpm_data))

tpm_long <- tpm_data %>%
  rownames_to_column(var = "MAG") %>%
  pivot_longer(cols = -MAG, names_to = "sample", values_to = "TPM")

plot_data_summary <- tpm_long %>%
  left_join(group_data, by = "sample") %>%
  filter(TPM > 0)

mags_small <- plot_data_summary %>%
  filter(group == "small") %>%
  pull(MAG) %>%
  unique()

mags_large <- plot_data_summary %>%
  filter(group == "large") %>%
  pull(MAG) %>%
  unique()

venn_list <- list(USF = mags_small, ConF = mags_large)

p_venn_fill <- ggvenn(venn_list, fill_color = c("#1B9E77", "#D95F02"), fill_alpha = 0.8, stroke_color = "black",
  stroke_size = 0.5, set_name_size = 10, text_size = 10, show_percentage = TRUE) + labs(title = "Number of detected candidate-pathogen ASVs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20))

ggsave("Figure_Venn_Sorted.pdf", p_venn_fill, width = 8, height = 8)

ggsave("Figure_Venn_Sorted.png", p_venn_fill, width = 8, height = 8, dpi = 300, bg = "white")

print("韦恩图已保存。")
