# fig_agora.R

library(ggplot2)
library(patchwork)
library(ggbeeswarm)
suppressMessages(library(tidyverse))
theme_set(theme_bw())

# table of information about the different microbes there are AGORA2 models for
agora_tax <- read_csv("figure_data/fig_agora_info.csv", show_col_types = F)
# summarized test results from the AGORA2 models
agora_data <- bind_rows(lapply(
  list.files("figure_data/", "fig_agora_data*", full.names = T),
  function(f) read_csv(f, show_col_types = F)
))
cat("Have data for ", nrow(agora_data), " models (", round(100*nrow(agora_data)/7302), "%)\n", sep = "")

plot_data <- agora_data %>%
  # strip the extensions off of the filenames to get AGORA2's "MicrobeIDs"
  mutate(organism = gsub(".mat", "", model)) %>%
  # replace counts with proportions, then drop the counts
  mutate(any_prop = flagged/all_rxns) %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(dil_prop = `dilution-blocked`/all_rxns) %>%
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  select(
    -model, -all_rxns, -flagged, -`dead-ends`, -`dilution-blocked`, -duplicates,
    -loops, -redoxes
  ) %>%
  # add in the other info about each organism
  merge(agora_tax, by.x = "organism", by.y = "MicrobeID") %>%
  pivot_longer(
    c(any_prop, dead_prop, dil_prop, dupe_prop, loop_prop),
    names_to = "test", values_to = "prop"
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
    TRUE ~ `Oxygen Requirement`
  )) %>%
  select(-`Oxygen Requirement`) %>%
  mutate(o2_req = factor(o2_req, levels = c(
    "Obligate aerobe",
    "Aerobe",
    "Microaerophile",
    "Facultative anaerobe",
    "Anaerobe",
    "Obligate anaerobe",
    "Other"
  ))) %>%
  # lump together the phyla with very few models
  group_by(Phylum) %>%
  mutate(Phylum = ifelse(n() <= 10, "Other", Phylum)) %>%
  ungroup() %>%
  mutate(Phylum = fct_relevel(as.factor(Phylum), "Other", after = Inf))

# do one violin plot with the phylum on the x-axis and another with the oxygen
# requirement (obligate anaerobe/obligate aerobe/facultative)
tax_panel <- ggplot(plot_data, aes(x = Phylum, y = prop)) +
  geom_quasirandom(aes(col = test), alpha = 0.4, show.legend = F) +
  stat_summary(fun = median, geom = "crossbar", linewidth = 0.25) +
  stat_summary(
    fun = quantile,
    fun.args = list(probs = 0.25),
    geom = "crossbar",
    linetype = "dashed",
    linewidth = 0.1
  ) +
  stat_summary(
    fun = quantile,
    fun.args = list(probs = 0.75),
    geom = "crossbar",
    linetype = "dashed",
    linewidth = 0.1
  ) +
  facet_grid(rows = vars(test), scales = "free_y") +
  scale_color_manual(
    values = c("#FDB462", "#FB8072", "#B3DE69", "#FCCDE5", "#80B1D3")
  ) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  ) +
  labs(x = "", y = "Proportion of Reactions Flagged by Test")

o2_panel <- ggplot(plot_data, aes(x = o2_req, y = prop)) +
  geom_quasirandom(aes(col = test), alpha = 0.4, show.legend = F) +
  stat_summary(fun = median, geom = "crossbar", linewidth = 0.25) +
  stat_summary(
    fun = quantile,
    fun.args = list(probs = 0.25),
    geom = "crossbar",
    linetype = "dashed",
    linewidth = 0.1
  ) +
  stat_summary(
    fun = quantile,
    fun.args = list(probs = 0.75),
    geom = "crossbar",
    linetype = "dashed",
    linewidth = 0.1
  ) +
  facet_grid(rows = vars(test), scales = "free_y") +
  scale_color_manual(
    values = c("#FDB462", "#FB8072", "#B3DE69", "#FCCDE5", "#80B1D3")
  ) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  ) +
  labs(x = "", y = "Proportion of Reactions Flagged by Test")

fig <- tax_panel | o2_panel
ggsave(
  "figures/fig_agora.tif", fig, width = 7.5, height = 7.5, units = "in",
  dpi = 300, compression = "lzw"
)

# quantify how different these distributions are
stats_data <- plot_data %>%
  filter((Phylum != "Other") & (o2_req != "Other"))

mean_diffs <- bind_rows(
  stats_data %>%
    group_by(test, Phylum) %>%
    summarize(mean_prop = mean(prop), .groups = "drop_last") %>%
    summarize(
      max_diff = max(mean_prop) - min(mean_prop),
      groups = paste(
        Phylum[which.max(mean_prop)], Phylum[which.min(mean_prop)], sep = "-"
      )
    ) %>%
    mutate(grouped_by = "phylum"),
  stats_data %>%
    group_by(test, o2_req) %>%
    summarize(mean_prop = mean(prop), .groups = "drop_last") %>%
    summarize(
      max_diff = max(mean_prop) - min(mean_prop),
      groups = paste(
        o2_req[which.max(mean_prop)], o2_req[which.min(mean_prop)], sep = "-"
      )
    ) %>%
    mutate(grouped_by = "aerobicity")
) %>% mutate(max_diff = round(max_diff, 2))

p_vals <- bind_rows(
  stats_data %>%
    group_by(test) %>%
    do(anova(lm(prop ~ Phylum, data = .))) %>%
    ungroup() %>%
    mutate(grouped_by = "phylum"),
  stats_data %>%
    group_by(test) %>%
    do(anova(lm(prop ~ o2_req, data = .))) %>%
    ungroup() %>%
    mutate(grouped_by = "aerobicity")
) %>%
  mutate(anova_q_val = signif(p.adjust(`Pr(>F)`), 2)) %>%
  filter(!is.na(anova_q_val)) %>%
  select(grouped_by, test, anova_q_val)

mean_diffs
#merge(mean_diffs, p_vals) %>%
#  select(grouped_by, test, groups, max_diff, anova_q_val) %>%
#  arrange(grouped_by, test)

# correlations
agora_data %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(dil_prop = `dilution-blocked`/all_rxns) %>%
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  select(
    -model, -all_rxns, -flagged, -`dead-ends`, -`dilution-blocked`, -duplicates,
    -loops, -redoxes
  ) %>%
  cor()
