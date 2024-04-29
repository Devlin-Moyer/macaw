# figs_2_and_S1.R

# load packages

library(jsonlite)
library(tools)
suppressMessages(library(tidyverse))
theme_set(theme_bw())

simplify_results <- function(df) {
  out <- df %>%
    mutate(
      duplicate_test = ifelse(
        (
          (duplicate_test_exact != "ok") |
          (duplicate_test_directions != "ok") |
          (duplicate_test_coefficients != "ok") |
          (duplicate_test_redox != "ok")
        ), "bad", "ok"
      ),
      dead_end_test = ifelse(grepl("(ok|only)", dead_end_test), "ok", "bad"),
      diphosphate_test = ifelse(diphosphate_test != "ok", "bad", "ok"),
      loop_test = ifelse(loop_test != "ok", "bad", "ok"),
      dilution_test = ifelse(grepl("(ok|always)", dilution_test), "ok", "bad")
    ) %>%
    select(
      -duplicate_test_exact, -duplicate_test_directions,
      -duplicate_test_coefficients, -duplicate_test_redox
    )
  return(out)
}

# use a bunch of tables to add columns with the top-level KEGG ortholog groups
# a dataframe of test results
add_kegg_groups <- function(
  df, rxn_to_gene, gene_to_orth, orth_to_grp = ortholog_to_group
) {
  df %>%
    left_join(rxn_to_gene, by = "reaction_id") %>%
    left_join(gene_to_orth, by = "gene_id") %>%
    left_join(orth_to_grp, by = "ortholog_id") %>%
    # no longer need the gene or ortholog IDs
    select(-gene_id, -ortholog_id) %>%
    # one of the KEGG ortholog groups is called "not included in regular maps"
    # which is uh not what we need to know rn, so change those to "none"
    mutate(kegg_group = ifelse(
      kegg_group == "Not included in regular maps", "Not in KEGG", kegg_group
    )) %>%
    # KEGG reactions can be associated with multiple different KEGG orthologs
    # (presumably cuz isozymes are a thing), but sometimes these different
    # orthologs are all in the same higher-level grouping of orthologs that
    # we're interested in here, so there are a bunch of duplicate rows rn
    unique() %>%
    mutate(kegg_group = ifelse(is.na(kegg_group), "Not in KEGG", kegg_group)) %>%
    # some reactions might be associated with some genes that are in KEGG's
    # ortholog hierarchy in addition to genes that aren't, and we don't want a
    # reaction to simultaneously have a group of "none" and a real group, so
    # filter out all the rows where kegg_group is "none" for reactions that also
    # have rows where kegg_group isn't none
    group_by(reaction_id) %>%
    filter((n_distinct(kegg_group) == 1) | (kegg_group != "Not in KEGG")) %>%
    ungroup() %>%
    # abbreviate the KEGG group names cuz many of them are quite long
    mutate(kegg_group = case_when(
      kegg_group == "Amino acid metabolism" ~ "Amino acids",
      kegg_group == "Biosynthesis of other secondary metabolites" ~
        "Secondary metabolites",
      kegg_group == "Carbohydrate metabolism" ~ "Carbohydrates",
      kegg_group == "Energy metabolism" ~ "Energy",
      kegg_group == "Glycan biosynthesis and metabolism" ~ "Glycans",
      kegg_group == "Lipid metabolism" ~ "Lipids",
      kegg_group == "Metabolism of cofactors and vitamins" ~ "Cofactors",
      kegg_group == "Metabolism of other amino acids" ~ "Other amino acids",
      kegg_group == "Metabolism of terpenoids and polyketides" ~ "Terpenoids",
      kegg_group == "Nucleotide metabolism" ~ "Nucleotides",
      kegg_group == "Xenobiotics biodegradation and metabolism" ~ "Xenobiotics",
      TRUE ~ kegg_group
    ))
}

# read in the test results for both versions of Yeast-GEM, iML1515, and both
# versions of Human-GEM
human15_tests <- read_csv(
  "figure_data/Human-GEMv1.15_test-results.csv", show_col_types = FALSE
) %>% simplify_results()
human18_tests <- read_csv(
  "figure_data/Human-GEMv1.18_test-results.csv", show_col_types = FALSE
) %>% simplify_results()
yeast_tests <- read_csv(
  "figure_data/yeast-GEMv9.0.0_test-results.csv", show_col_types = FALSE
) %>% simplify_results()
ecoli_tests <- read_csv(
  "figure_data/iML1515_test-results.csv", show_col_types = FALSE
) %>% simplify_results()

# get tables of reaction IDs and NCBI gene IDs
human_to_gene <- read_csv(
  "figure_data/Human-GEMv1.15_reactions-to-genes.csv", show_col_types = FALSE
)
yeast_to_gene <- read_csv(
  "figure_data/yeast-GEMv9.0.0_reactions-to-genes.csv", show_col_types = FALSE
)
ecoli_to_gene <- read_csv(
  "figure_data/iML1515_reactions-to-genes.csv", show_col_types = FALSE
)

# get tables of gene IDs and KEGG ortholog IDs
human_to_ortholog <- read.delim(
  "input_data/Human-GEMv1.15_genes-to-KEGG.tsv", header = FALSE
) %>%
  # remove the prefixes KEGG added to all IDs just to inconvenience us
  separate(V1, c(NA, "ortholog_id"), sep = ":", extra = "drop") %>%
  separate(V2, c(NA, "gene_id"), sep = ":", extra = "drop")

yeast_to_ortholog <- read.delim(
  "input_data/yeast-GEMv9.0.0_genes-to-KEGG.tsv", header = FALSE
) %>%
  separate(V1, c(NA, "ortholog_id"), sep = ":", extra = "drop") %>%
  separate(V2, c(NA, "gene_id"), sep = ":", extra = "drop")

ecoli_to_ortholog <- read.delim(
  "input_data/iML1515_genes-to-KEGG.tsv", header = FALSE
) %>%
  separate(V1, c(NA, "ortholog_id"), sep = ":", extra = "drop") %>%
  separate(V2, c(NA, "gene_id"), sep = ":", extra = "drop")

# turn JSON file with hierarchy of all KEGG functional orthologs into a two-
# column dataframe we can use to map individual KEGG ortholog IDs to the top-
# level groupings of KEGG orthologs deemed to be "metabolic"
ortholog_to_group <- fromJSON("input_data/KEGG_ortholog_hierarchy.json") %>%
  # this hierarchy includes a bunch of genes that are not metabolic enzymes, so
  # start by filtering down to those
  pluck("children") %>%
  pull(children) %>%
  first() %>%
  # now simplify the structure here cuz we don't actually care about anything
  # other than the relationships between individual KEGG orthologs and these
  # categories/labels right below the "metabolism" branch
    mutate(ortholog_id = unlist(lapply(
    children, function(x) x %>%
      pull(children) %>%
      bind_rows() %>%
      separate("name", c("some_ids"), sep = " ", extra = "drop") %>%
      summarize(more_ids = toString(some_ids))
  ))) %>%
  select(-children) %>%
  mutate(ortholog_id = strsplit(ortholog_id, ", ")) %>%
  unnest(ortholog_id) %>%
  # strip those numerical prefixes off of the names of the ortholog groups
  mutate(kegg_group = gsub("\\d+\ ", "", name, perl = TRUE)) %>%
  select(-name)

# use all the tables to label each reaction with a KEGG ortholog group
human15_kegg <- add_kegg_groups(human15_tests, human_to_gene, human_to_ortholog)
human18_kegg <- add_kegg_groups(human18_tests, human_to_gene, human_to_ortholog)
yeast_kegg <- add_kegg_groups(yeast_tests, yeast_to_gene, yeast_to_ortholog)
ecoli_kegg <- add_kegg_groups(ecoli_tests, ecoli_to_gene, ecoli_to_ortholog)

fig_2_data <- bind_rows(
  human15_kegg %>% mutate(Model = "Human-GEM"),
  yeast_kegg %>% mutate(Model = "Yeast-GEM"),
  ecoli_kegg %>% mutate(Model = "iML1515")
) %>%
  group_by(Model, kegg_group) %>%
  summarize(
    dilution_pct = as.integer(100*sum(dilution_test == "bad")/n()),
    loop_pct = as.integer(100*sum(loop_test == "bad")/n()),
    duplicate_pct = as.integer(100*sum(duplicate_test == "bad")/n()),
    dead_end_pct = as.integer(100*sum(dead_end_test == "bad")/n()),
    diphosphate_pct = as.integer(100*sum(diphosphate_test == "bad")/n()),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(dilution_pct, loop_pct, duplicate_pct, dead_end_pct, diphosphate_pct),
    names_to = "test", values_to = "pct_rxns"
  ) %>%
  mutate(test = gsub("_pct", " test", test)) %>%
  mutate(test = gsub("dead_end", "dead-end", test)) %>%
  mutate(kegg_group = relevel(as.factor(kegg_group), "Not in KEGG")) %>%
  mutate(Model = factor(Model, c("Human-GEM", "Yeast-GEM", "iML1515"))) %>%
  mutate(test = toTitleCase(test)) %>%
  mutate(test = factor(test, c(
    "Dilution Test", "Loop Test", "Dead-End Test", "Duplicate Test",
    "Diphosphate Test"
  )))

# make a separate dataframe of coordinates to use for grey rectangles to put
# behind every other group of bars instead of using normal gridlines
rect_df <- data.frame(
  starts = seq(0.5, n_distinct(fig_2_data$kegg_group), by = 2),
  ends = seq(1.5, n_distinct(fig_2_data$kegg_group), by = 2)
)

fig_2 <- ggplot() +
  scale_x_discrete() +
  geom_rect(
    data = rect_df,
    mapping = aes(xmin = starts, xmax = ends),
    ymin = -Inf,
    ymax = Inf,
    fill = "gray",
    color = NA,
    alpha = 0.5
  ) +
  geom_bar(
    data = fig_5_data,
    mapping = aes(x = kegg_group, y = pct_rxns, fill = Model),
    stat = "identity",
    position = "dodge"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_grid(rows = vars(test), scales = "free_y") +
  labs(x = "Metabolism of ___", y = "% Reactions in Group") +
  theme(
    text = element_text(color = "black", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(color = "black", size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.clip = "off",
    panel.grid = element_blank(),
    legend.key.size = unit(1/8, "in"),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0, "in")
  )

ggsave(
  "figures/fig_2.png", fig_2, height = 5.25, width = 3.25, units = "in", dpi = 600
)

fig_S1_data <- bind_rows(
  human15_kegg %>% mutate(Version = "1.15"),
  human18_kegg %>% mutate(Version = "1.18")
) %>%
  group_by(Version, kegg_group) %>%
  summarize(
    dilution_pct = as.integer(100*sum(dilution_test == "bad")/n()),
    loop_pct = as.integer(100*sum(loop_test == "bad")/n()),
    duplicate_pct = as.integer(100*sum(duplicate_test == "bad")/n()),
    dead_end_pct = as.integer(100*sum(dead_end_test == "bad")/n()),
    diphosphate_pct = as.integer(100*sum(diphosphate_test == "bad")/n()),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(dilution_pct, loop_pct, duplicate_pct, dead_end_pct, diphosphate_pct),
    names_to = "test", values_to = "pct_rxns"
  ) %>%
  mutate(kegg_group = relevel(as.factor(kegg_group), "Not in KEGG")) %>%
  mutate(test = gsub("_pct", " test", test)) %>%
  mutate(test = gsub("dead_end", "dead-end", test)) %>%
  mutate(test = toTitleCase(test)) %>%
  mutate(test = factor(test, c(
    "Dilution Test", "Loop Test", "Dead-End Test", "Duplicate Test",
    "Diphosphate Test"
  )))

fig_S1 <- ggplot() +
  scale_x_discrete() +
  geom_rect(
    data = rect_df,
    mapping = aes(xmin = starts, xmax = ends),
    ymin = -Inf,
    ymax = Inf,
    fill = "gray",
    color = NA,
    alpha = 0.5
  ) +
  geom_bar(
    data = fig_S1_data,
    mapping = aes(x = kegg_group, y = pct_rxns, fill = Version),
    stat = "identity",
    position = "dodge",
    width = 0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_grid(rows = vars(test), scales = "free_y") +
  labs(x = "Metabolism of ___", y = "% Reactions in Group") +
  theme(
    text = element_text(color = "black", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.clip = "off",
    panel.grid = element_blank(),
    legend.key.size = unit(1/8, "in"),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0, "in")
  )

ggsave(
  "figures/fig_S1.png", height = 5.25, width = 3, units = "in", plot = fig_S1,
  dpi = 600
)
