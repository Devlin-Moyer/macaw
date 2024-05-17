# fig_agora.R

library(tidyverse)
theme_set(theme_bw())

# table of information about the different microbes there are AGORA2 models for
agora_tax <- read_csv("input_data/AGORA2_taxonomy.csv", show_col_types = F)

# numbers of reactions flagged by each test in the relatively few AGORA2 models
# that we could read in successfully
agora_data <- read_csv("figure_data/fig_agora.csv", show_col_types = F)
plot_data <- agora_data %>%
  # strip the .xml off of the filenames to get AGORA2's "MicrobeIDs"
  mutate(organism = gsub(".xml", "", model)) %>%
  # replace counts with proportions, then drop the counts
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  select(-model, -duplicates, -`dead-ends`, -loops, -all_rxns) %>%
  # add in the other info about each organism
  merge(agora_tax, by.x = "organism", by.y = "MicrobeID") %>%
  pivot_longer(
    c(dupe_prop, dead_prop, loop_prop), names_to = "test", values_to = "prop"
  ) %>%
  mutate(test = case_when(
    test == "dupe_prop" ~ "Duplicate Test",
    test == "dead_prop" ~ "Dead-End Test",
    test == "loop_prop" ~ "Loop Test"
  )) %>%
  filter(Phylum %in% c(
    # only plot results from models in phyla we had at least 100 models for
    "Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"
  ))

ggplot(plot_data, aes(x = Phylum, y = prop, fill = test)) +
  geom_violin(show.legend = F) +
  facet_grid(rows = vars(test), scales = "free_y") +
  scale_fill_manual(values = c("#FFFF33", "#FB8072", "#80B1D3")) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  ) +
  labs(x = "", y = "Proportion of Reactions Flagged by Test")

ggsave("figures/fig_agora.png", width = 3, height = 6, units = "in")