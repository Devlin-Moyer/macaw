# fig_mendoza.R

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(tidyverse))
theme_set(theme_bw())

# read in summarized test results from the models
mendoza_data <- read_csv("figure_data/fig_mendoza_data.csv", show_col_types = F)

# reformat data to prepare for making a heatmap
heatmap_data <- mendoza_data %>%
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
row_names <- c(
  "BPE_any_prop", "LPL_any_prop", "PPU_any_prop",
  "BPE_dead_prop", "LPL_dead_prop", "PPU_dead_prop",
  "BPE_dil_prop", "LPL_dil_prop", "PPU_dil_prop",
  "BPE_dupe_prop", "LPL_dupe_prop", "PPU_dupe_prop",
  "BPE_loop_prop", "LPL_loop_prop", "PPU_loop_prop"
)
org_names <- c("B. pertussis", "L. plantarum", "P. putida")
org_ann_bar <- structure(
  rep(org_names, 5), names = row_names
)
test_ann_bar <- structure(
  c(
    rep("Any", 3),
    rep("Dead-end", 3),
    rep("Dilution", 3),
    rep("Duplicate", 3),
    rep("Loop", 3)
  ), names = row_names
)

# make two more named vectors to control colors in annotation bars
org_colors <- structure(c("#BEBADA", "#FFFFB3", "#8DD3C7"), names = org_names)
test_colors <- structure(
  c("#FDB462", "#FD8072", "#B3DE69", "#FCCDE5", "#80B1D3"),
  names = c("Any", "Dead-end", "Dilution", "Duplicate", "Loop")
)
# make colormap for the main heatmap
col_fun = colorRamp2(c(0, max(heatmap_data)), c("white", "red"))

# finally make the heatmap
tiff(
  "figures/fig_mendoza.tif", height = 6.5, width = 4.25, unit = "in",
  res = 300, compression = "lzw"
)
Heatmap(
  heatmap_data, name = "Mean Proportion\nof Reactions\nFlagged by Test",
  cluster_rows = F, cluster_columns = F,
  show_row_names = F,
  column_title = "Reconstruction Method", column_title_side = "bottom",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  left_annotation = rowAnnotation(
    Test = test_ann_bar,
    Organism = org_ann_bar,
    col = list(Organism = org_colors, Test = test_colors)
  ),
  height = unit(5, "in"), width = unit(2.25, "in")
)
dev.off()
