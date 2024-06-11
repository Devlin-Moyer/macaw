# fig_5.R

# load packages
lib <- "/usr3/graduate/dcmoyer/R/x86_64-pc-linux-gnu-library/4.2"
library(png, lib.loc = lib)
library(ggplot2, lib.loc = lib)
library(ggpubr, lib.loc = lib)
library(patchwork, lib.loc = lib)
suppressMessages(library(tidyverse, lib.loc = lib))

# read in a PNG file and turn it into a ggplot object so it can be patchworked
# together with actual ggplot plots into a single figure
load_image_as_panel <- function(path) {
  img <- readPNG(path)
  panel <- ggplot() + background_image(img) + theme_void() +
    coord_fixed(ratio = dim(img)[1] / dim(img)[2]) +
    theme(plot.margin = unit(c(0,0,0,0), "in"))
  return(panel)
}

# read in data
fig_5c_data <- read_csv("figure_data/fig_5c_data.csv", show_col_types = FALSE)

fig_5c <- fig_5c_data %>%
  # add a newline to "PDH blocked, no impact on GCS"
  mutate(impact = gsub(", no", ",\nno", impact)) %>%
  # make sure impact and condition have factor levels in a helpful order
  mutate(impact = factor(
    impact,
    c("No impact on either", "PDH blocked,\nno impact on GCS", "Both blocked")
  )) %>%
  mutate(condition = factor(
    condition,
    c(
      "Human-GEM 1.18+ with Dilution",
      "Human-GEM 1.18+ without Dilution",
      "Human-GEM 1.15 with Dilution",
      "Human-GEM 1.15 without Dilution",
      "Patients"
    )
  )) %>%
  # make a heatmap
  ggplot(aes(x = gene, y = condition, fill = impact)) +
    geom_tile(col = "black", linewidth = 0.5) +
    # Juan said to use a pastel color palette
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "", fill = "Measured or Predicted\nImpact of Knockout") +
    theme(
      # drop the background and axes
      panel.background = element_blank(),
      # remove extra space on the left and bottom of the plot from where the
      # axis titles would be
      plot.margin = unit(c(1,1,-3,-3), "mm"),
      # remove axis ticks and the extra space between the labels and the plot
      # area left behind without them
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "mm"),
      # make all text font size 8 and make sure there aren't any awkward gaps
      # between things
      axis.text = element_text(color = "black", size = 8),
      axis.text.y = element_text(margin = margin(r = -0.1)),
      axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 0)),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.box.spacing = unit(0, "mm"),
      axis.title.x = element_blank(),
      # add a bit more space between each item in the legend cuz we have two
      # lines of text for one of the things
      legend.spacing.y = unit(2, "mm")
    ) +
    # without this, legend.spacing.y only affects the space between the title
    # of the legend and the contents of the legend and not the space between
    # each row of the legend
    guides(fill = guide_legend(byrow = TRUE))

# now that we've made the plot, read in the Escher maps to use as the other two
# panels and patchwork them all together
fig_5a <- load_image_as_panel("figures/fig_5a.png")
fig_5b <- load_image_as_panel("figures/fig_5b.png")

fig_5 <- (fig_5a / fig_5b / free(fig_5c)) +
  plot_layout(heights = unit(c(1.8, 2.75, 1.5), "in")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold", hjust = 1, vjust = -5))

ggsave(
  "figures/fig_5.tif",
  fig_5,
  height = 6.6,
  width = 5,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)
