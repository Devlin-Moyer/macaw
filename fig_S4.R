# fig_S4.R
# plot the proportions of reactions flagged by each test against the version
# number for each version of Human-GEM to see how that changed over time

suppressMessages(library(tidyverse))
theme_set(theme_bw())

all_human_data <- read_csv(
  "figure_data/fig_S4_data.csv", show_col_types = FALSE, col_types = "fiiiii"
)
plot_data <- all_human_data %>%
  # replace counts with proportions, then drop the counts
  mutate(any_prop = flagged/all_rxns) %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(dil_prop = `dilution-blocked`/all_rxns) %>%
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  select(
    -flagged, -`dead-ends`, -`dilution-blocked`, duplicates, -loops, -all_rxns
  ) %>%
  # pivot cuz it'll be easier to plot this way
  pivot_longer(-model_version, names_to = "test", values_to = "prop") %>%
  # make more human-readable for figure
  mutate(test = case_when(
    test == "any_prop" ~ "Any Test",
    test == "dead_prop" ~ "Dead-End Test",
    test == "dil_prop" ~ "Dilution Test",
    test == "dupe_prop" ~ "Duplicate Test",
    test == "loop_prop" ~ "Loop Test"
  ))

fig <- ggplot(
  plot_data, aes(x = model_version, y = prop, col = test, group = test)
) +
  geom_line() +
  scale_x_discrete(
    # tick mark every other version
    labels = c(
      "1.0", "", "1.2", "", "1.4", "", "1.6", "", "1.8", "", "1.10", "", "1.12",
      "", "1.14", "", "1.16", "", "1.18"
    )
  ) +
  labs(
    x = "Version of Human-GEM",
    y = "Proportion of All Reactions",
    col = "Flagged by"
  ) +
  theme(
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.title.align = 0.5,
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.box.spacing = unit(0, "in")
  )

ggsave(
  "figures/fig_S4.png", height = 2, width = 4, units = "in", dpi = 600
)
