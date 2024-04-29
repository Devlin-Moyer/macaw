# fig_3b.R

suppressMessages(library(tidyverse))
theme_set(theme_bw())

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
        (dilution_test != "ok") | (loop_test != "ok") |
        (duplicate_test != "ok") | (dead_end_test != "ok") |
        (diphosphate_test != "ok")
      ), nrow(df) + row_num, pathway
    )) %>%
    select(-row_num, -reaction_equation) %>%
    # now get the number of reactions in each pathway
    filter(pathway != 0) %>%
    group_by(pathway) %>%
    mutate(pathway_size = n()) %>%
    ungroup()
}

get_cumul_rxns <- function(df, test) {
  # if test was "any", keep all rows
  if (test != "any") {
    df <- filter(df, !!as.symbol(test) == "bad")
  }
  df %>%
    group_by(pathway_size) %>%
    summarize(reactions = n(), .groups = "drop") %>%
    mutate(cumul_rxns = cumsum(reactions)) %>%
    mutate(flagged_by = toTitleCase(case_when(
      test == "any" ~ "Any Test",
      test == "dead_end_test" ~ "Dead-End Test",
      TRUE ~ gsub("_", " ", test)
    )))
}

results <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
)
tests <- c(
  "any", "dilution_test", "loop_test", "dead_end_test", "duplicate_test",
  "diphosphate_test"
)

fig <- bind_rows(lapply(
  tests, function(t) get_cumul_rxns(get_pathway_sizes(results), t)
)) %>%
  ggplot(aes(x = pathway_size, y = cumul_rxns, col = flagged_by)) +
    geom_line() +
    scale_x_log10() +
    labs(
      x = "# Reactions in Pathway",
      y = "Reactions in Pathways\nwith up to _ Reactions",
      col = "Reactions Flagged By"
    )

ggsave("figures/fig_3b.png", height = 1, width = 1.5, units = "in", dpi = 600)
