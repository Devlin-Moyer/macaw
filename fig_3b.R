# fig_3b.R

library(tools)
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

results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
) %>%
  simplify_results() %>%
  get_pathway_sizes()

tests <- c(
  "any", "dilution_test", "loop_test", "dead_end_test", "duplicate_test",
  "diphosphate_test"
)
plot_data <- bind_rows(lapply(
  tests, function(t) filter_sizes(results, t)
))

fig <- plot_data %>%
  ggplot(aes(x = pathway_size, col = flagged_by)) +
    geom_density(aes(y = (after_stat(scaled) * after_stat(n)) + 1)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Size of Pathway",
      y = "# Pathways + 1",
      col = "Pathway Had\n\u22651 Reaction\nFlagged By"
    ) +
    theme(
      text = element_text(color = "black", size = 8),
      axis.text = element_text(color = "black", size = 6),
      panel.grid = element_blank(),
      legend.key.size = unit(1/8, "in"),
      legend.title.align = 0.5,
      legend.box.spacing = unit(0, "in")
    )

ggsave("figures/fig_3b.png", height = 1.5, width = 3.5, units = "in", dpi = 600)
