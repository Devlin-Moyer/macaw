# shared_funcs.R

# get proportions of reactions in each KEGG group flagged by each or any test
summ_func <- function(df) {
  summarize(
    df,
    any_pct = as.integer(
      100*sum(
        (dead_end_test == "bad") |
          (dilution_test == "bad") |
          (duplicate_test == "bad") |
          (loop_test == "bad")
      )/n()
    ),
    dead_end_pct = as.integer(100*sum(dead_end_test == "bad")/n()),
    dilution_pct = as.integer(100*sum(dilution_test == "bad")/n()),
    duplicate_pct = as.integer(100*sum(duplicate_test == "bad")/n()),
    loop_pct = as.integer(100*sum(loop_test == "bad")/n()),
    .groups = "drop"
  )
}