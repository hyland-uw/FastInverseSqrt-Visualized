#!/usr/bin/env Rscript
# Convert CSV artifacts into internal data/ objects for the visualfrsr package.

extdata_dir <- file.path("inst", "extdata")
data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

load_dataset <- function(filename) {
  read.csv(
    file.path(extdata_dir, filename),
    stringsAsFactors = FALSE
  )
}

equal05_20bins <- load_dataset("equal05_20bins.csv")
equal_bins <- load_dataset("equal_bins.csv")
manybins <- load_dataset("manybins.csv")
too_many_bins <- load_dataset("too_many_bins.csv")
varied_bins <- load_dataset("varied_bins.csv")

varied_bins$Bottom <- factor(
  varied_bins$Bottom,
  levels = c(
    "point00781", "point0156", "point0312", "point0625",
    "point0125", "point25", "point5", "one"
  ),
  labels = c(0.0078, 0.0156, 0.0312, 0.0623, 0.125, 0.25, 0.5, 1),
  ordered = TRUE
)

output_objects <- list(
  equal05_20bins = equal05_20bins,
  equal_bins = equal_bins,
  manybins = manybins,
  too_many_bins = too_many_bins,
  varied_bins = varied_bins
)

purrr::iwalk(output_objects, function(obj, name) {
  save(obj,
       file = file.path(data_dir, paste0(name, ".rda")),
       compress = "xz")
})
