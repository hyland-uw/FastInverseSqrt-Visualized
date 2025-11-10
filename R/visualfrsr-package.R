#' visualfrsr: Visual Explorations of Fast Reciprocal Square Root Approximations
#'
#' Provides data generators, binning utilities, clustering helpers, and
#' standardized plotting functions that were previously scattered across scripts
#' in the FastInverseSqrt-Visualized project. The package exposes those
#' capabilities as reusable functions so datasets and visualizations can be
#' reproduced programmatically or integrated into pipelines.
#' @import ggplot2 frsrr
#' @importFrom viridis turbo viridis_pal
#' @importFrom dplyr filter mutate group_by summarise summarize arrange select
#'   rename row_number ungroup slice slice_min n ntile bind_rows bind_cols
#'   distinct if_else transmute pull lag
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom data.table setDT rbindlist
#' @importFrom purrr map_dfr map2_dfr
#' @importFrom tibble tibble
#' @importFrom grDevices colorRampPalette
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom scales trans_breaks rescale
#' @importFrom stats runif sd kmeans na.omit setNames
#' @importFrom utils read.csv
#' @importFrom rlang sym
#' @importFrom Rglpk Rglpk_solve_LP
#' @keywords internal
"_PACKAGE"
