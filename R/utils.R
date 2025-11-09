# utils.R

# Load required packages
library(tibble)
library(data.table)
library(Matrix)
library(Rglpk)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(scales)
library(stringr)
library(lpSolve)
library(purrr)
library(knitr)
library(frsrr)

#### Colors and such

## from https://stackoverflow.com/a/9568659/1188479
## useful for false categorical coloring
false_categorical_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

#### helper functions

# Define the logStratifiedSampler function
logStratifiedSampler <- function(min, max, n) {
  exp(runif(n, log(min), log(max)))
}

# Slicing Function
create_slices <- function(df, N, min_input = 0.5, max_input = 2.0) {
  slice_width <- (max_input - min_input) / N

  df %>%
    filter(input >= min_input, input < max_input) %>%
    mutate(slice = cut(input,
                       breaks = seq(min_input, max_input, by = slice_width),
                       labels = FALSE,
                       include.lowest = TRUE))
}

# Optimization Function
find_optimal_magic <- function(slice_data) {
  unique_magics <- unique(slice_data$magic)

  results <- sapply(unique_magics, function(m) {
    slice_data %>%
      filter(magic == m) %>%
      summarise(max_error = max(error)) %>%
      pull(max_error)
  })

  optimal_index <- which.min(results)
  list(minimum = unique_magics[optimal_index], objective = results[optimal_index])
}

## more compact function for annotation
mc_annotate <- function(magic_value, label,
                        color, x_start = -0.035, x_end = 0.036,
                        text_size = 8) {
  list(
    annotate("segment",
             x = x_start, xend = x_end,
             y = magic_value, yend = magic_value,
             color = color, linetype = 2, linewidth = 1.5),
    annotate("point", x = x_end, y = magic_value, color = color, size = 3),
    annotate("text", x = x_end + 0.002, y = magic_value, label = label,
             hjust = -0.05, vjust = 0.5, color = color, size = text_size)
  )
}

## force specific plotting order so we can plot low iterations "on top of"
## higher iterations so overplotting doesn't cover the optimal range
# required there to be a variable "iters"
# which we can walk over
create_geom_points <- function(data, iter_range, shape, size, alpha = 1) {
  lapply(iter_range, function(i) {
    geom_point(data = data[data$iters == i, ],
               shape = shape,
               size = size,
               alpha = alpha)
  })
}
