# fig_2c.R

library(ggupset)
library(tidyverse)
theme_set(theme_bw())

condense_results <- function(row) {
  dupe <- (row["duplicate_test_exact"] != "ok") |
    (row["duplicate_test_directions"] != "ok") |
    (row["duplicate_test_coefficients"] != "ok") |
    (row["duplicate_test_redox"] != "ok")
  dead <- !grepl("(ok|only)", row["dead_end_test"])
  diphos <- (row["diphosphate_test"] != "ok")
  loop <- (row["loop_test"] != "ok")
  dil <- !grepl("(ok|always)", row["dilution_test"])
  result <- character()
  if (dil) {result[length(result) + 1] <- "Dilution Test"}
  if (loop) {result[length(result) + 1] <- "Loop Test"}
  if (dupe) {result[length(result) + 1] <- "Duplicate Test"}
  if (dead) {result[length(result) + 1] <- "Dead-end Test"}
  if (diphos) {result[length(result) + 1] <- "Diphosphate Test"}
  return(result)
}

results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = F
)
results$flagged_by <- apply(results, 1, condense_results)
fig <- results %>%
  filter(sapply(results$flagged_by, length) > 0) %>%
  ggplot(aes(x = flagged_by)) +
    geom_bar(fill = "black", width = 0.5) +
    scale_x_upset() +
    #scale_y_log10() +
    labs(x = "", y = "# Reactions") +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(vjust = -20),
      axis.text.y = element_text(color = "black"),
      plot.margin = unit(c(0.125, 0.125, -0.125, -0.125), "in")
    ) +
    theme_combmatrix(combmatrix.label.text = element_text(color = "black"))

ggsave(
  "figures/fig_2c.png", height = 2.25, width = 3.25, units = "in", plot = fig,
  dpi = 600
)