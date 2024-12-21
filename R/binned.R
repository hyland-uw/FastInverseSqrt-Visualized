source("utils.R")

generate_bins<- function(sizes = c(4, 8, 16, 32, 64, 128, 256), floats, ints, NRmax) {
  do.call(rbind, lapply(sizes, function(n) {
    frsr_bin(n_bins = n, NRmax = NRmax, x_max = 0.5,
             float_samples = floats, magic_samples = ints,
             magic_min = 0x5F290000, magic_max = 0x5F410000)
  }))
}

binnedNR <- generate_bins(floats = 4096, ints = 16384, NRmax = 1)

binnedGuess <- generate_bins(floats = 4096, ints = 16384, NRmax = 0)


## helper function for linear programming
find_optimal_buckets <- function(binned, M) {
  # Create objective vector (minimize total max error)
  obj <- binned$Max_Relative_Error

  # Create constraint matrix
  n_slices <- nrow(binned)

  # Constraint 1: Each input value must be covered exactly once
  input_points <- sort(unique(c(binned$Range_Min, binned$Range_Max)))
  coverage_matrix <- matrix(0, nrow = length(input_points) - 1, ncol = n_slices)

  for(i in 1:(length(input_points) - 1)) {
    mid_point <- (input_points[i] + input_points[i + 1]) / 2
    coverage_matrix[i,] <- as.numeric(
      binned$Range_Min <= mid_point & binned$Range_Max >= mid_point
    )
  }

  # Constraint 2: Use exactly M slices
  slice_constraint <- matrix(1, nrow = 1, ncol = n_slices)

  # Combine constraints
  const.mat <- rbind(coverage_matrix, slice_constraint)
  const.dir <- c(rep("==", nrow(coverage_matrix)), "==")
  const.rhs <- c(rep(1, nrow(coverage_matrix)), M)

  # Solve using lpSolve
  result <- lp("min", obj, const.mat, const.dir, const.rhs, all.bin = TRUE)

  # Extract selected slices
  selected <- binned[result$solution == 1, ]

  return(selected %>% arrange(Range_Min))
}


### plots

## pseudo-waterfall

compare_n_values <- function(binned, n_small, n_large) {
  # Create a combined dataset with error values from both N
  combined_data <- binned %>%
    filter(N %in% c(n_small, n_large)) %>%
    mutate(N_type = if_else(N == n_small, "small", "large")) %>%
    pivot_wider(
      names_from = N_type,
      values_from = c(Range_Min, Range_Max, Max_Relative_Error),
      id_cols = c(Magic)
    ) %>%
    # Add reference values by matching ranges
    filter(!is.na(Range_Min_small) & !is.na(Range_Min_large))

  ggplot() +
    # Horizontal segments for smaller N
    geom_segment(data = filter(binned, N == n_small),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Max_Relative_Error, yend = Max_Relative_Error),
                 color = "red", linewidth = 0.5) +
    # Horizontal segments for larger N
    geom_segment(data = filter(binned, N == n_large),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Max_Relative_Error, yend = Max_Relative_Error),
                 color = "blue", linewidth = 0.5) +
    # Rectangles showing difference
    geom_rect(data = filter(binned, N == n_large),
              aes(xmin = Range_Min, xmax = Range_Max,
                  ymin = Max_Relative_Error,
                  ymax = filter(binned, N == n_small)$Max_Relative_Error[
                    findInterval(Range_Min,
                                 filter(binned, N == n_small)$Range_Min)]),
              fill = "blue", alpha = 0.3) +
    labs(title = sprintf("Comparison of Max Relative Error N=%d vs N=%d",
                         n_small, n_large),
         x = "Input Value",
         y = "Max Relative Error")
}

compare_buckets <- function(bucket1, bucket2, error = "max") {
  error_col <- if(error == "max") "Max_Relative_Error" else "Avg_Relative_Error"

  # Get all unique x-points where rectangles should start or end
  x_points <- sort(unique(c(
    bucket1$Range_Min, bucket1$Range_Max,
    bucket2$Range_Min, bucket2$Range_Max
  )))

  # Create rectangles for each adjacent pair of x-points
  rectangles <- tibble(
    xmin = x_points[-length(x_points)],
    xmax = x_points[-1]
  ) %>%
    mutate(
      # Find error values for each range
      error1 = sapply(xmin, function(x) {
        bucket1[[error_col]][x >= bucket1$Range_Min & x < bucket1$Range_Max][1]
      }),
      error2 = sapply(xmin, function(x) {
        bucket2[[error_col]][x >= bucket2$Range_Min & x < bucket2$Range_Max][1]
      }),
      # Determine rectangle properties
      ymin = pmin(error1, error2),
      ymax = pmax(error1, error2),
      fill_color = if_else(error2 < error1, "blue", "red")
    )

  # Create horizontal segments for both buckets
  segments_h1 <- bucket1 %>%
    select(Range_Min, Range_Max, Error = !!sym(error_col))

  segments_h2 <- bucket2 %>%
    select(Range_Min, Range_Max, Error = !!sym(error_col))

  ggplot() +
    # Horizontal segments for both buckets
    geom_segment(data = segments_h1,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error),
                 color = "black", linewidth = 0.5) +
    geom_segment(data = segments_h2,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error),
                 color = "black", linewidth = 0.5) +
    # Rectangles
    geom_rect(data = rectangles,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill_color),
              alpha = 0.3) +
    scale_fill_identity() +
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25),
                       limits = c(0.25, 1.0)) +
    labs(title = "Comparison of Error Values Between Base and Extension",
         x = "Input Value",
         y = paste(error, "Relative Error")) +
    theme_minimal()
}


## working single slice visualization
plot_bucket <- function(df) {
  # Create horizontal segments dataset
  segments_h <- df %>%
    pivot_longer(
      cols = c(Avg_Relative_Error, Max_Relative_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) %>%
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Relative_Error", "Max_Relative_Error"),
                               labels = c("Avg Error", "Max Error")))

  # Create vertical segments dataset - now including both ends of each range
  segments_v <- bind_rows(
    # Segments for the ending ranges
    df %>%
      transmute(
        x = Range_Max,
        y_start = Avg_Relative_Error,
        y_end = Max_Relative_Error
      ),
    # Segments for the starting ranges
    df %>%
      transmute(
        x = Range_Min,
        y_start = Avg_Relative_Error,
        y_end = Max_Relative_Error
      )
  )
  ggplot() +
    # Horizontal segments
    geom_segment(data = segments_h,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error,
                     color = Error_Type),
                 linewidth = 0.5) +
    # Vertical segments at range breaks
    geom_segment(data = segments_v,
                 aes(x = x, xend = x,
                     y = y_start, yend = y_end),
                 linetype = "dotted",
                 color = "black",
                 linewidth = 0.25,
                 alpha = 0.75) +
    scale_color_manual(values = c("Avg Error" = "blue", "Max Error" = "red")) +
    labs(x = "Input Value",
         y = "Error") +
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25),
                       limits = c(0.25, 1.0)) +
    theme_minimal()
}

## working multiple via facet_wrap
plot_multiple_n <- function(binned, n_values = unique(binned$N)) {
  # Create horizontal segments dataset
  segments_h <- binned %>%
    filter(N %in% n_values) %>%
    pivot_longer(
      cols = c(Avg_Relative_Error, Max_Relative_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) %>%
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Relative_Error", "Max_Relative_Error"),
                               labels = c("Avg Error", "Max Error")))

  # Create vertical segments dataset
  segments_v <- bind_rows(
    # Segments for the ending ranges
    binned %>%
      filter(N %in% n_values) %>%
      group_by(N) %>%
      slice(1:(n()-1)) %>%
      transmute(
        N = N,
        x = Range_Max,
        y_start = Avg_Relative_Error,
        y_end = Max_Relative_Error
      ),
    # Segments for the starting ranges
    binned %>%
      filter(N %in% n_values) %>%
      group_by(N) %>%
      slice(2:n()) %>%
      transmute(
        N = N,
        x = Range_Min,
        y_start = Avg_Relative_Error,
        y_end = Max_Relative_Error
      )
  )

  ggplot() +
    # Horizontal segments
    geom_segment(data = segments_h,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error,
                     color = Error_Type),
                 linewidth = 0.5) +
    # Vertical segments at range breaks
    geom_segment(data = segments_v,
                 aes(x = x, xend = x,
                     y = y_start, y_end = y_end),
                 linetype = "dotted",
                 color = "black",
                 linewidth = 0.25,
                 alpha= 0.75) +
    facet_grid(cols = vars(N)) +
    scale_color_manual(values = c("Avg Error" = "blue", "Max Error" = "red")) +
    labs(title = "Error Trends Across Input Range",
         x = "Input Value",
         y = "Error") +
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25),
                       limits = c(0.25, 1.0)) +
    theme_minimal()
}


## "weight" error by the fraction of the domain
## it covers.
norm_errorN <- function(df, bins) {
  bucket <- find_optimal_buckets(df, bins)
  bucket$Width <- (bucket[, "Range_Max"] - bucket[, "Range_Min"]) / 0.75
  bucket <- bucket %>%
    # if you don't pick N, dplyr complains and does it anyway
    select(N, Max_Relative_Error, Width) %>%
    rename(Error = Max_Relative_Error)
  ## avoids privileging small bucket sizes
  sum_err <- with(bucket, sum(Error*Width))
  return(sum_err)
}

# Use purrr::map_dfr to apply the function to each bin value
# why this is better than a for loop is unclear
map_dfr(4:36, ~tibble(bins = .x, error = norm_errorN(binned, .x))) %>%
  ggplot(aes(x = bins, y = error)) +
  geom_line() +
  geom_point() +
  labs(x = "Bins",
       y = "Normalized Error",
       title = "Optimal bucket selection error reduction slows after 24 bins") +
  scale_x_continuous(breaks = seq(4, 36, by = 4))
