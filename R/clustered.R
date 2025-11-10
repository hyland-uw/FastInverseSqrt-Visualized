cluster_magic_constants <- function(data, n_clusters = 3) {
  grouped <- data |>
    group_by(magic) |>
    summarize(
      mean_error = mean(error),
      min_error = min(error),
      max_error = max(error),
      .groups = "drop"
    )

  kmeans_result <- kmeans(grouped[, c("mean_error", "min_error", "max_error")],
                          centers = n_clusters)
  grouped$cluster <- kmeans_result$cluster

  grouped |>
    group_by(cluster) |>
    summarize(
      average_magic = mean(magic),
      min_magic = min(magic),
      max_magic = max(magic),
      count = n(),
      .groups = "drop"
    ) |>
    pivot_longer(
      cols = c(average_magic, min_magic, max_magic),
      names_to = "statistic",
      values_to = "value"
    ) |>
    mutate(statistic = factor(
      statistic,
      levels = c("average_magic", "min_magic", "max_magic"),
      labels = c("Average", "Minimum", "Maximum")
    ))
}

#' Build representative cluster bands across the main error regimes.
#'
#' @param n_floats Number of floats sampled per cluster sweep.
#' @param n_ints Number of magic integers sampled in each sweep.
#' @return Tibble describing the min/max bands for three illustrative clusters.
#' @export
build_cluster_bands <- function(
    n_floats = 2048,
    n_ints = 16384) {
  low_err <- frsr_sample(
    n_floats = n_floats,
    NRmax = 1,
    keep_params = TRUE,
    x_min = 0.5,
    x_max = 2.0,
    n_ints = n_ints,
    magic_min = 1596980000L,
    magic_max = 1598050000L
  ) |>
    select(input, error, magic)

  wide_err <- frsr_sample(
    n_floats = n_floats,
    NRmax = 1,
    keep_params = TRUE,
    x_min = 0.5,
    x_max = 2.0,
    n_ints = n_ints,
    magic_min = 1596910000L,
    magic_max = 1598100000L
  ) |>
    select(input, error, magic)

  low_clusters <- cluster_magic_constants(low_err) |>
    filter(cluster %in% c(1, 2))

  high_cluster <- cluster_magic_constants(wide_err) |>
    filter(cluster == 1) |>
    mutate(cluster = 3)

  rbind(low_clusters, high_cluster) |>
    filter(statistic != "Average") |>
    select(!count)
}

#' Sample representative clusters for downstream visualizations.
#'
#' @param n Number of samples per cluster.
#' @param cluster_ranges List of length three containing magic ranges.
#' @return Tibble with input/error/magic/cluster columns.
#' @export
sample_clusters <- function(
    n = 16384,
    cluster_ranges = list(
      c(1597622357, 1597176038),
      c(1598049921, 1597918602),
      c(1597884263, 1596910052)
    )) {
  representative_clusters <- map2_dfr(
    cluster_ranges,
    seq_along(cluster_ranges),
    function(bounds, cluster_id) {
      frsr_sample(
        n = n,
        NRmax = 1,
        keep_params = TRUE,
        x_min = 0.5,
        x_max = 2.0,
        magic_max = bounds[1],
        magic_min = bounds[2]
      ) |>
        select(input, error, magic) |>
        mutate(cluster = factor(cluster_id))
    }
  )

  representative_clusters
}
