# fig_4.R

suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
library(ggplot2)
library(patchwork)
library(ggbeeswarm)
suppressMessages(library(tidyverse))
theme_set(theme_bw())

fig_4a_data_raw <- read_csv("figure_data/fig_4a_data.csv", show_col_types = F)
fig_4a_data <- fig_4a_data_raw %>%
  # unpack information in model filenames into separate columns
  separate(model, c("organism", "method", NA)) %>%
  mutate(method = case_when(
    grepl("AU", method) ~ "AuReMe",
    grepl("CA", method) ~ "CarveMe",
    grepl("ME", method) ~ "Merlin",
    grepl("MD", method) ~ "MetaDraft",
    grepl("MS", method) ~ "ModelSEED",
    grepl("PT", method) ~ "Pathway Tools",
    grepl("RA", method) ~ "RAVEN"
  )) %>%
  # turn numbers of reactions into proportions of all reactions
  mutate(any_prop = flagged/all_rxns) %>%
  mutate(dead_prop = `dead-ends`/all_rxns) %>%
  mutate(dil_prop = `dilution-blocked`/all_rxns) %>%
  mutate(dupe_prop = duplicates/all_rxns) %>%
  mutate(loop_prop = loops/all_rxns) %>%
  # drop count columns now that we have proportions
  select(
    -flagged, -`dead-ends`, -`dilution-blocked`, -duplicates, -loops, -all_rxns
  ) %>%
  # pivot twice to get into appropriate format for heatmap
  pivot_longer(
    c(any_prop, dead_prop, dil_prop, dupe_prop, loop_prop),
    names_to = "test",
    values_to = "prop"
  ) %>%
  pivot_wider(
    id_cols = c(organism, test),
    names_from = "method",
    values_from = "prop",
    # average together proportions from the different parameter combinations for
    # each reconstruction method
    values_fn = mean
  ) %>%
  # get rows in a helpful order and give them unique rownames
  arrange(test, organism) %>%
  unite(org_test, c(organism, test)) %>%
  column_to_rownames("org_test") %>%
  # put columns in alphabetical order
  select(AuReMe, CarveMe, MetaDraft, Merlin, ModelSEED, `Pathway Tools`, RAVEN)

# make two named vectors to make annotation bars
fig_4a_row_names <- c(
  "BPE_any_prop", "LPL_any_prop", "PPU_any_prop",
  "BPE_dead_prop", "LPL_dead_prop", "PPU_dead_prop",
  "BPE_dil_prop", "LPL_dil_prop", "PPU_dil_prop",
  "BPE_dupe_prop", "LPL_dupe_prop", "PPU_dupe_prop",
  "BPE_loop_prop", "LPL_loop_prop", "PPU_loop_prop"
)
fig_4a_orgs <- c("B. pertussis", "L. plantarum", "P. putida")
fig_4a_org_ann_bar <- structure(
  rep(fig_4a_orgs, 5), names = fig_4a_row_names
)
fig_4a_test_ann_bar <- structure(
  c(
    rep("Any", 3),
    rep("Dead-end", 3),
    rep("Dilution", 3),
    rep("Duplicate", 3),
    rep("Loop", 3)
  ), names = fig_4a_row_names
)
# make two more named vectors to control colors in annotation bars
fig_4a_org_colors <- structure(
  c("#BEBADA", "#FFFFB3", "#8DD3C7"), names = fig_4a_orgs
)
fig_4a_test_colors <- structure(
  c("#FDB462", "#FD8072", "#B3DE69", "#FCCDE5", "#80B1D3"),
  names = c("Any", "Dead-end", "Dilution", "Duplicate", "Loop")
)
# make colormap for the main heatmap
fig_4a_col_fun = colorRamp2(c(0, max(fig_4a_data)), c("white", "red"))

# finally make the heatmap
fig_4a <- Heatmap(
  fig_4a_data,
  col = fig_4a_col_fun,
  rect_gp = gpar(col = "black", lwd = 1.5),
  cluster_rows = F,
  cluster_columns = F,
  left_annotation = rowAnnotation(
    Test = fig_4a_test_ann_bar,
    Organism = fig_4a_org_ann_bar,
    col = list(Organism = fig_4a_org_colors, Test = fig_4a_test_colors),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
      Organism = list(
        title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)
      ),
      Test = list(
        title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)
      )
    )
  ),
  show_row_names = F,
  column_names_gp = gpar(fontsize = 8),
  column_title = "Reconstruction Method",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 8),
  name = "Avg. Prop. of\nReactions\nFlagged by\nTest",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)
  ),
  height = unit(2.9, "in"),
  width = unit(1.4, "in")
) 

# summarized test results from the AGORA2 models
fig_4b_data_raw <- bind_rows(lapply(
  list.files("figure_data/", "fig_4b_data*", full.names = T),
  function(f) read_csv(f, show_col_types = F)
))
cat("Have data for ", nrow(fig_4b_data_raw), " models (", round(100*nrow(fig_4b_data_raw)/7302), "%)\n", sep = "")

fig_4b_data <- fig_4b_data_raw %>%
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
  ))

fig_4b <- ggplot(fig_4b_data, aes(x = prop, fill = test)) +
  geom_histogram(bins = 50, show.legend = F) +
  facet_wrap(. ~ test, scales = "free", ncol = 1) +
  scale_fill_manual(
    values = c("#8DD3C7", "#FB8072", "#B3DE69", "#FCCDE5", "#80B1D3")
  ) +
  theme(
    strip.background = element_blank(),
    strip.clip = "off",
    panel.grid = element_blank(),
    text = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    panel.spacing = unit(0, "in")
  ) +
  labs(x = "Prop. of Reactions Flagged by Test", y = "Number of GSMMs")

# turn ComplexHeatmap object into a grob so wrap_plots recognizes it
fig_4a <- fig_4a %>%
  draw(padding = unit(c(0,0,0,0), "in")) %>%
  grid.grabExpr() %>%
  wrap_plots()

fig_4 <- (fig_4a | fig_4b) +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.03,0.98)
  )

ggsave(
  "figures/fig_4.tif",
  fig_4,
  width = 5.5,
  height = 4.75,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)
