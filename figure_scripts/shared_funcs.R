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

# read in a PNG file and turn it into a ggplot object so it can be patchworked
# together with actual ggplot plots into a single figure
load_image_as_panel <- function(path) {
  img <- png::readPNG(path)
  asp_rat <- dim(img)[1] / dim(img)[2]
  panel <- ggplot() + background_image(img) + theme_void() +
    coord_fixed(ratio = asp_rat)
  # return both the ggplot object and the aspect ratio
  return(list(panel, asp_rat))
}
