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

# read in a PNG file and turn it into a ggplot object so it can be patchworked
# together with actual ggplot plots into a single figure
load_image_as_panel <- function(path) {
  img <- readPNG(path)
  asp_rat <- dim(img)[1] / dim(img)[2]
  panel <- ggplot() + background_image(img) + theme_void() +
    coord_fixed(ratio = asp_rat)
  # return both the ggplot object and the aspect ratio
  return(list(panel, asp_rat))
}

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
      dilution_test = ifelse(grepl("(ok|always)", dilution_test), "ok", "bad"),
      loop_test = ifelse(loop_test != "ok", "bad", "ok"),
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
  tests <- c("dead_end_test", "dilution_test", "duplicate_test", "loop_test")
  hist_data <- bind_rows(lapply(
    c("any", tests), function(t) filter_sizes(better_results, t)
  ))
  c50s <- data.frame(
    flagged_by = c(
      "Dead-End Test", "Dilution Test", "Duplicate Test", "Loop Test"
    ), label = sapply(tests, function(t) compute_c50(better_results, t))
  )
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
    scale_fill_manual(
      values = c("#FDB462", "#FB8072", "#B3DE69", "#FCCDE5", "#80B1D3")
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
  loop <- (row["loop_test"] != "ok")
  dil <- !grepl("(ok|always)", row["dilution_test"])
  result <- character()
  if (dil) {result[length(result) + 1] <- "Dilution Test"}
  if (loop) {result[length(result) + 1] <- "Loop Test"}
  if (dupe) {result[length(result) + 1] <- "Duplicate Test"}
  if (dead) {result[length(result) + 1] <- "Dead-end Test"}
  return(result)
}

make_upset <- function(data) {
  data$flagged_by <- apply(data, 1, condense_results)
  fig <- data %>%
    filter(sapply(data$flagged_by, length) > 0) %>%
    ggplot(aes(x = flagged_by)) +
      geom_bar(fill = "black", width = 0.5) +
      scale_x_upset() +
      labs(x = "", y = "# Reactions") +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(color = "black"),
        plot.background = element_rect(fill = "transparent")
      ) +
      theme_combmatrix(combmatrix.label.text = element_text(color = "black"))
  return(fig)
}

fig_3_width <- 6
fig_3bc_height <- 2
human_results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
) %>% select(-diphosphate_test)
fig_3_stuff <- load_image_as_panel("figures/fig_3a.png")
fig_3a_height <- (fig_3_width * (fig_3_stuff[[2]]))
fig_3a <- fig_3_stuff[[1]] + theme(plot.tag.position = c(0.002,0.98))
fig_3b <- make_hists(human_results, 500, 1500) +
  theme(
    plot.tag.position = c(-0.22, 1.1),
    plot.margin = unit(c(0,0,0.125,0), "in")
  )
fig_3c <- make_upset(human_results) +
  theme(
    plot.tag.position = c(-0.4, 0.96),
    plot.margin = unit(c(0, 0, 0.125, 0.5), "in")
  )
fig_3 <- (free(fig_3a) / (free(fig_3b) | free(fig_3c))) +
  plot_layout(heights = unit(c(fig_3a_height, fig_3bc_height), "in")) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 12),
    plot.tag.location = "panel"
  )
ggsave(
  "figures/fig_3.tif",
  height = fig_3a_height + fig_3bc_height,
  width = fig_3_width,
  units = "in",
  dpi = 300
)
