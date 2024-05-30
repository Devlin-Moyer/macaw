# fig_3_S2_S3.R

lib <- "/usr3/graduate/dcmoyer/R/x86_64-pc-linux-gnu-library/4.2"
library(tools)
library(scales)
library(ggupset, lib.loc = lib)
library(png)
library(ggpubr, lib.loc = lib)
library(patchwork)
suppressMessages(library(tidyverse))
theme_set(theme_bw())

simplify_results <- function(df) {
  out <- df %>%
    mutate(
      # condense the duplicate test columns into one
      duplicate_test = ifelse(
        (
          (duplicate_test_exact != "ok") |
          (duplicate_test_directions != "ok") |
          (duplicate_test_coefficients != "ok") |
          (duplicate_test_redox != "ok")
        ), "bad", "ok"
      ),
      # make sure all the other columns just say "bad" or "ok"
      dead_end_test = ifelse(grepl("(ok|only)", dead_end_test), "ok", "bad"),
      diphosphate_test = ifelse(diphosphate_test != "ok", "bad", "ok"),
      loop_test = ifelse(loop_test != "ok", "bad", "ok"),
      dilution_test = ifelse(grepl("(ok|always)", dilution_test), "ok", "bad")
    ) %>%
    # drop all the specific duplicate test columns
    select(
      -duplicate_test_exact, -duplicate_test_directions,
      -duplicate_test_coefficients, -duplicate_test_redox
    )
  return(out)
}

get_pathway_sizes <- function(df) {
  df %>%
    # for all reactions that weren't in a pathway (have a 0 in the pathway
    # column) that were flagged by at least one test, change the pathway to the
    # total number of rows plus that row's row number so it gets a unique 
    # pathway label
    rownames_to_column(var = "row_num") %>%
    mutate(row_num = as.integer(row_num)) %>%
    mutate(pathway = ifelse(
      (pathway == 0) & (
        (dead_end_test != "ok") |
        (dilution_test != "ok") |
        (diphosphate_test != "ok") |
        (duplicate_test != "ok") |
        (loop_test != "ok")
      ), nrow(df) + row_num, pathway
    )) %>%
    select(-row_num, -reaction_equation) %>%
    # now get the number of reactions in each pathway
    filter(pathway != 0) %>%
    group_by(pathway) %>%
    mutate(pathway_size = n()) %>%
    ungroup()
}

filter_sizes <- function(df, test) {
  # if test was "any", keep all rows
  if (test != "any") {
    df <- filter(df, !!as.symbol(test) != "ok")
  }
  out <- df %>%
    # now that we're sure we only have reactions in pathways with at least one
    # reaction that was flagged by the given test, only keep one row for each
    # pathway with its size
    select(pathway, pathway_size) %>%
    unique() %>%
    mutate(flagged_by = toTitleCase(case_when(
      test == "any" ~ "Any Test",
      test == "dead_end_test" ~ "Dead-End Test",
      TRUE ~ gsub("_", " ", test)
    ))) %>%
    select(-pathway)
  return(out)
}

# for all reactions flagged by a particular test, get the list of the numbers
# of reactions in each pathway containing each reaction and find the median.
# i.e. if 4 of the 10 reactions in a particular pathway were flagged by the
# given test, you'd add 4 10s to the list of pathway sizes you're taking the
# median of. so it's kind of like a weighted median (Juan is calling it a C50
# because it apparently resembles a similarly named statistic used in the
# context of long read sequencing data)
compute_c50 <- function(df, test) {
  df %>%
    filter(!!as.symbol(test) != "ok") %>%
    pull(pathway_size) %>%
    median()
}

make_hists <- function(results, x_coord, y_coord) {
  better_results <- results %>%
    simplify_results() %>%
    get_pathway_sizes()
  tests <- c(
    "dilution_test", "loop_test", "dead_end_test", "duplicate_test",
    "diphosphate_test"
  )
  hist_data <- bind_rows(lapply(
    c("any", tests), function(t) filter_sizes(better_results, t)
  ))
  c50s <- bind_cols(
    c(
      "Dilution Test", "Loop Test", "Dead-End Test", "Duplicate Test",
      "Diphosphate Test"
    ), sapply(tests, function(t) compute_c50(better_results, t))
  )
  colnames(c50s) <- c("flagged_by", "label")
  c50s <- c50s %>%
    # add columns for the coordinates of each label
    mutate(x = x_coord, y = y_coord) %>%
    mutate(label = paste0("C50 = ", label))
  hist_fig <- ggplot() +
    geom_histogram(
      data = hist_data,
      mapping = aes(x = pathway_size, fill = flagged_by),
      bins = 30,
      show.legend = F
    ) +
    geom_text(
      data = c50s,
      mapping = aes(x = x, y = y, label = label),
      hjust = 1,
      size = 6 * 0.36 # approximate a font size of 6
    ) +
    facet_wrap(. ~ flagged_by) +
    scale_x_log10() +
    scale_y_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = c(1, 10, 100, 1000)
    ) +
    labs(
      x = "# Reactions in Pathway",
      y = "# Pathways"
    ) +
    theme(
      text = element_text(color = "black", size = 8),
      axis.text = element_text(color = "black", size = 8),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.clip = "off"
    )
  return(hist_fig)
}

# condense all the different test result columns into a single one that lists
# all tests each reaction was flagged by
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

# start with figure 3
human_results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
)
fig_3b <- make_hists(human_results, 500, 1500)

human_results$flagged_by <- apply(human_results, 1, condense_results)
fig_3c <- human_results %>%
  filter(sapply(human_results$flagged_by, length) > 0) %>%
  ggplot(aes(x = flagged_by)) +
    geom_bar(fill = "black", width = 0.5) +
    scale_x_upset() +
    labs(x = "", y = "# Reactions") +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(color = "black"),
      plot.margin = unit(c(0.125, 0.125, -0.125, 0.5), "in")
    ) +
    theme_combmatrix(combmatrix.label.text = element_text(color = "black"))

# now that we've created the ggplot objects for panels b and c, read in panel a
# (that we made by hand in Cytoscape) and patchwork all three together
fig_3a_raw <- readPNG("figures/fig_3a.png")
fig_3a <- ggplot() + background_image(fig_3a_raw) + theme_void() +
  # preserve aspect ratio so patchwork doesn't distort the image
  coord_fixed(ratio = dim(fig_3a_raw)[1] / dim(fig_3a_raw)[2])

(fig_3a / (free(fig_3b) | fig_3c + plot_layout(widths = unit(3, "in")))) +
  plot_layout(heights = c(4.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, hjust = 1, vjust = -5))
ggsave(
  "figures/fig_3.png", height = 10, width = 8, units = "in", dpi = 600
)

# then do figures S2 and S3
# TODO: patchwork together the histograms with the Cytoscape networks
fig_S2b <- make_hists(
  read_csv(
    "figure_data/yeast-GEMv9.0.0_test-results.csv", show_col_types = FALSE
  ), 1300, 300
)
ggsave("figures/fig_S2b.png")
fig_S3b <- make_hists(
  read_csv("figure_data/iML1515_test-results.csv", show_col_types = FALSE),
  1000, 125
)
ggsave("figures/fig_S3b.png")
