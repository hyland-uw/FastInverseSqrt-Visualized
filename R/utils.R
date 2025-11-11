#' Internal color palette for categorical plots.
#' @keywords internal
#' @noRd
.false_categorical_25 <- c(
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

# Legacy `%>%` pipelines were replaced with base pipes (`|>`)

#' Log-stratified sampler for floats.
#' 
#' @param min Minimum sample value.
#' @param max Maximum sample value.
#' @param n Number of points to draw.
#' @return Numeric vector of length `n`.
logStratifiedSampler <- function(min, max, n) {
  exp(runif(n, log(min), log(max)))
}

#' Slice frsr sample data into equal-width bins.
#'
#' @param df Data frame containing an `input` column.
#' @param N Number of slices.
#' @param min_input Minimum input to keep.
#' @param max_input Maximum input to keep.
#' @return Tidy data frame with `slice` column.
create_slices <- function(df, N, min_input = 0.5, max_input = 2.0) {
  slice_width <- (max_input - min_input) / N

  df |>
    filter(input >= min_input, input < max_input) |>
    mutate(slice = cut(input,
                       breaks = seq(min_input, max_input, by = slice_width),
                       labels = FALSE,
                       include.lowest = TRUE))
}

#' Find the optimal magic constant for a slice.
#'
#' @param slice_data Slice produced by [create_slices()].
#' @return List containing the minimum magic and its objective value.
find_optimal_magic <- function(slice_data) {
  unique_magics <- unique(slice_data$magic)

  results <- sapply(unique_magics, function(m) {
    slice_data |>
      filter(magic == m) |>
      summarise(max_error = max(error)) |>
      pull(max_error)
  })

  optimal_index <- which.min(results)
  list(minimum = unique_magics[optimal_index], objective = results[optimal_index])
}
