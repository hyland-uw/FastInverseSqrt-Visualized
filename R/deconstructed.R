# internal helper to ensure equal-sized blocks
sample_block_rows <- function(df, block_size = 2048) {
  limit <- nrow(df) - (nrow(df) %% block_size)
  if (limit <= 0) {
    return(df)
  }
  df[sample(limit), ]
}

#' Prepare a deconstruction dataframe for plotting.
#' 
#' @param df Output from [frsrr::frsr_sample()] with `keep_params = TRUE`.
#' @param NRmax Newton-Raphson iteration ceiling used during sampling.
#' @return Tidy data frame with derived ranking columns.
prep_deconstruction_df <- function(df, NRmax) {
  df |>
    sample_block_rows() |> # sample evenly divisible block
    filter(iters <= (NRmax - 1)) |>
    mutate(
      reference = 1 / sqrt(input),
      iters = factor(iters, levels = 1:max(iters), labels = 1:max(iters)),
      input_rank = ntile(input, 128),
      magic_rank = ntile(magic, 128),
      initial_rank = ntile(initial - reference, 128),
      after_rank = ntile(after_one - reference, 128),
      iter_rank = cut(as.numeric(iters),
                      breaks = c(0, 1, 2, 3, 4, 24, 60, 94),
                      labels = c("1", "2", "3", "4", "5-24", "25-60", "61-94"))
    )
}

#' Sample datasets used by the deconstruction visualizations.
#'
#' @param tolerance Convergence tolerance passed to `frsr_sample()`.
#' @param NRmax Iteration ceiling for the main and wide samples.
#' @param narrow_NRmax Iteration ceiling for the zoom sample.
#' @param magic_core Range of integer magics for the main sample.
#' @param magic_wide Range of integer magics for the wide sample.
#' @param magic_narrow Range of integer magics for the zoom sample.
#' @param n_core,n_wide,n_narrow Sample sizes for each dataset.
#' @return List containing `deconstructed`, `widened`, and `narrowed` data frames.
#' @export
generate_deconstruction_samples <- function(
    tolerance = 2^-15,
    NRmax = 95,
    narrow_NRmax = 2,
    magic_core = c(1.5935e9, 1601175552),
    magic_wide = c(1.3e9, 1.6025e9),
    magic_narrow = c(1597461647, 1597469647),
    n_core = 524288,
    n_wide = 262144,
    n_narrow = 131072) {
  deconstructed <- frsr_sample(
    n = n_core,
    NRmax = NRmax,
    tol = tolerance,
    magic_min = magic_core[1],
    magic_max = magic_core[2],
    keep_params = TRUE
  )

  widened <- frsr_sample(
    n = n_wide,
    NRmax = NRmax,
    tol = tolerance,
    magic_min = magic_wide[1],
    magic_max = magic_wide[2],
    keep_params = TRUE
  )

  narrowed <- frsr_sample(
    n = n_narrow,
    NRmax = narrow_NRmax,
    magic_min = magic_narrow[1],
    magic_max = magic_narrow[2],
    keep_params = TRUE
  ) |>
    mutate(reference = 1 / sqrt(input))

  list(
    deconstructed = prep_deconstruction_df(deconstructed, NRmax),
    widened = prep_deconstruction_df(widened, NRmax),
    narrowed = narrowed
  )
}
