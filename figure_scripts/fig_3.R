# fig_3.R

lib <- "/usr3/graduate/dcmoyer/R/x86_64-pc-linux-gnu-library/4.2"
library(tools, lib.loc = lib)
library(png, lib.loc = lib)
library(ggplot2, lib.loc = lib)
library(ggpubr, lib.loc = lib)
library(ggupset, lib.loc = lib)
library(scales, lib.loc = lib)
library(patchwork, lib.loc = lib)
suppressMessages(library(tidyverse, lib.loc = lib))
theme_set(theme_bw())

source("figure_scripts/shared_funcs.R")

fig_3_width <- 6
fig_3bc_height <- 2.25
human_results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
) %>% select(-diphosphate_test)
fig_3_stuff <- load_image_as_panel("figures/fig_3a.png")
fig_3a_height <- (fig_3_width * (fig_3_stuff[[2]]))
fig_3a <- fig_3_stuff[[1]] + theme(plot.tag.position = c(0.002,0.98))
fig_3b <- make_hists(human_results, 500, 1500) +
  theme(
    panel.spacing = unit(3, "mm"),
    plot.tag.position = c(-0.16, 1.1),
    plot.margin = unit(c(0,0,0.125,0), "in")
  )
fig_3c <- make_upset(human_results) +
  theme(
    plot.tag.position = c(-0.4, 0.96),
    plot.margin = unit(c(0, 0, 0.125, 0.3), "in")
  )
fig_3 <- (free(fig_3a) / (
  (free(fig_3b) | free(fig_3c)) + plot_layout(widths = c(1.7, 1))
)) +
  plot_layout(heights = unit(c(fig_3a_height, fig_3bc_height), "in")) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.location = "panel"
  )
ggsave(
  "figures/fig_3.tif",
  height = fig_3a_height + fig_3bc_height,
  width = fig_3_width,
  units = "in",
  dpi = 300,
  compression = "lzw"
)
