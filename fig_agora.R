# fig_agora.R

library(patchwork)
suppressMessages(library(tidyverse))
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
  mutate(any_prop = flagged/all_rxns) %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(dil_prop = `dilution-blocked`/all_rxns) %>%
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  select(
    -model, -flagged, -`dead-ends`, `dilution-blocked`, duplicates, -loops,
    -all_rxns
  ) %>%
  # add in the other info about each organism
  merge(agora_tax, by.x = "organism", by.y = "MicrobeID") %>%
  pivot_longer(
    c(dupe_prop, dead_prop, loop_prop), names_to = "test", values_to = "prop"
  ) %>%
  mutate(test = case_when(
    test == "any_prop" ~ "Any Test",
    test == "dead_prop" ~ "Dead-End Test",
    test == "dil_prop" ~ "Dilution Test",
    test == "dupe_prop" ~ "Duplicate Test",
    test == "loop_prop" ~ "Loop Test"
  )) %>%
  # clean up the oxygen requirement column
  mutate(o2_req = case_when(
    `Oxygen Requirement` == "Facultative anaerobic" ~ "Facultative anaerobe",
    `Oxygen Requirement` == "Facultative" ~ "Facultative anaerobe",
    `Oxygen Requirement` == "Obligate  anaerobe" ~ "Obligate anaerobe",
    `Oxygen Requirement` == "Aerotolerant" ~ "Other",
    `Oxygen Requirement` == "0" ~ "Other",
    is.na(`Oxygen Requirement`) ~ "Other",
  )) %>%
  mutate(o2_req = factor(o2_req, levels = c(
    "Obligate aerobe",
    "Aerobe",
    "Microaerophile",
    "Facultative anaerobe",
    "Anaerobe",
    "Obligate anaerobe",
    "Other"
  ))

# do one violin plot with the phylum on the x-axis and another with the oxygen
# requirement (obligate anaerobe/obligate aerobe/facultative)
tax_panel <- ggplot(plot_data, aes(x = Phylum, y = prop, fill = test)) +
  geom_violin(show.legend = F) +
  facet_grid(rows = vars(test), scales = "free_y") +
  scale_fill_manual(values = c("#FFFF33", "#FB8072", "#80B1D3")) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  ) +
  labs(x = "", y = "Proportion of Reactions Flagged by Test")

o2_panel <- ggplot(plot_data, aes(x = o2_req, y = prop, fill = test)) +
  geom_violin(show.legend = F) +
  facet_grid(rows = vars(test), scales = "free_y") +
  scale_fill_manual(values = c("#FFFF33", "#FB8072", "#80B1D3")) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  ) +
  labs(x = "", y = "Proportion of Reactions Flagged by Test")

fig <- tax_panel/o2_panel
ggsave("figures/fig_agora.png", fig, width = 7.5, height = 10, units = "in")
