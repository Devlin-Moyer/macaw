# fig_4.R

library(png)
library(ggplot2)
library(ggpubr)
library(patchwork)

# read in a PNG file and turn it into a ggplot object so it can be patchworked
load_image_as_panel <- function(path) {
  img <- readPNG(path)
  asp_rat <- dim(img)[1] / dim(img)[2]
  panel <- ggplot() + background_image(img) + theme_void() +
    coord_fixed(ratio = dim(img)[1] / dim(img)[2]) +
    theme(plot.margin = unit(c(0,0,0,0), "in"))
  return(panel)
}

fig_4a <- load_image_as_panel("figures/fig_4a.png")
fig_4b <- load_image_as_panel("figures/fig_4b.png")
fig_4c <- load_image_as_panel("figures/fig_4c.png")
fig_4 <- (fig_4a | fig_4b | fig_4c) +
  plot_layout(widths = c(1.3, 2, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold", vjust = -5))
ggsave(
  "figures/fig_4.tif",
  fig_4,
  height = 3,
  width = 7,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)
