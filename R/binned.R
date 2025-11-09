source("utils.R")



### helper functions
frsr_bin_multi <- function(sizes = c(4, 6, 8, 12, 18, 24, 32, 48, 64, 96, 128), x_min = 1.0, x_max = 2.0) {
  # Define constants inside the function
  floats <- 16384
  ints <- 524288
  NRmax <- 1
  magic_min <- 1596980000L
  magic_max <- 1598050000L

  do.call(rbind, lapply(sizes, function(n) {
    frsr_bin(n_bins = n,
             NRmax = NRmax,
             x_min = x_min,
             x_max = x_max,
             float_samples = floats,
             magic_samples = ints,
             magic_min = magic_min,
             magic_max = magic_max)
  }))
}

## helper function for linear programming
find_optimal_buckets <- function(binned, M, chunk_size = 1e4) {
  process_chunk <- function(chunk, M) {
    obj <- chunk$Max_Error
    n_slices <- nrow(chunk)

    input_points <- sort(unique(c(chunk$Range_Min, chunk$Range_Max)))
    n_points <- length(input_points) - 1

    # Pre-allocate vectors for sparse matrix creation
    max_nonzero <- n_points * n_slices
    i <- integer(max_nonzero)
    j <- integer(max_nonzero)
    x <- numeric(max_nonzero)
    idx <- 1

    mid_points <- (input_points[-1] + input_points[-length(input_points)]) / 2

    # Vectorized creation of coverage matrix
    for (k in seq_along(mid_points)) {
      indices <- which(chunk$Range_Min <= mid_points[k] & chunk$Range_Max >= mid_points[k])
      n_indices <- length(indices)
      if (n_indices > 0) {
        i[idx:(idx + n_indices - 1)] <- rep(k, n_indices)
        j[idx:(idx + n_indices - 1)] <- indices
        x[idx:(idx + n_indices - 1)] <- rep(1, n_indices)
        idx <- idx + n_indices
      }
    }

    # Trim excess pre-allocated space
    i <- i[1:(idx-1)]
    j <- j[1:(idx-1)]
    x <- x[1:(idx-1)]

    coverage_matrix <- sparseMatrix(i = i, j = j, x = x, dims = c(n_points, n_slices))

    slice_constraint <- Matrix(1, nrow = 1, ncol = n_slices, sparse = TRUE)

    const.mat <- rbind(coverage_matrix, slice_constraint)
    const.dir <- c(rep("==", nrow(coverage_matrix)), "==")
    const.rhs <- c(rep(1, nrow(coverage_matrix)), M)

    result <- Rglpk_solve_LP(obj = obj,
                             mat = as.matrix(const.mat),
                             dir = const.dir,
                             rhs = const.rhs,
                             types = rep("B", n_slices),
                             max = FALSE)

    selected <- chunk[result$solution == 1]

    return( selected[order(Range_Min)] )
  }
  setDT(binned)
  n_total <- nrow(binned)

  if (n_total > chunk_size) {
    # Pre-allocate list for results
    n_chunks <- ceiling(n_total / chunk_size)
    chunk_results <- vector("list", n_chunks)

    # Process data in chunks
    for (i in seq_len(n_chunks)) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_total)
      chunk <- binned[start_idx:end_idx]
      chunk_results[[i]] <- process_chunk(chunk, M)
    }

    # Combine results from all chunks & Final optimization on the combined results
    process_chunk(rbindlist(chunk_results), M)
  } else {
    # if chunking isn't needed, process right away
    process_chunk(binned, M)
  }
}


bucket_selection <- function(df, range = 4:100) {
  compute_statistics <- function(bucket_df, M) {
    mean_N <- mean(bucket_df$N)
    sd_N <- sd(bucket_df$N)

    bucket_df %>%
      mutate(
        Depth = N / mean_N,
        Width = Range_Max - Range_Min,
        Midpoint = (Range_Max + Range_Min) / 2,
        Left = Range_Min,
        Variance = ((N - mean_N) / sd_N)^2,
        M = M,
        Error = Max_Error
      ) %>%
      dplyr::select(M, Location, N, Depth, Error, Variance, Width, Midpoint, Left)
  }

  deviance_list <- lapply(range, function(M) {
    bucket_df <- find_optimal_buckets(df, M = M, chunk_size = 1e4)
    compute_statistics(bucket_df, M)
  })

  output <- do.call(rbind, deviance_list)
  # Calculate rarity across all buckets
  location_n_table <- with(output, table(Location, N))
  rarity_matrix <- location_n_table / sum(location_n_table)

  # Add Rarity to the output dataframe
  output$Rarity <- sapply(1:nrow(output), function(i) {
    rarity_matrix[as.character(output$Location[i]), as.character(output$N[i])]
  })

  output %>%
    group_by(M) %>%
    arrange(M, Midpoint) %>%
    mutate(order = row_number()) %>%
    ungroup() %>%
    as.data.frame()
}

analyze_bin_changes <- function(df) {
  # Create bins
  bins <- seq(0.5, 2.0, length.out = 2049)  # 2049 points to create 2048 bins

  # Assign data to bins
  df$bin <- cut(df$input,
                breaks = bins,
                labels = FALSE,
                include.lowest = TRUE)

  # Calculate bin statistics
  bin_stats <- df %>%
    group_by(bin) %>%
    summarize(
      mean_error = mean(error),
      mean_magic = mean(magic),
      .groups = 'drop'
    ) %>%
    # Calculate changes between consecutive bins
    mutate(
      delta_error = mean_error - lag(mean_error),
      delta_magic = mean_magic - lag(mean_magic)
    ) %>%
    # Remove first row which has NA in deltas
    filter(!is.na(delta_error))

  return(bin_stats)
}


### performance of magic constant
### outside their bins
compute_oob_performance <- function(data) {
  frsr_oob <- function(magic, bin_min, bin_max, total_min, total_max) {
    # Sample the lower out-of-bin range
    lower_samples <- frsr_sample(
      n = 1000,  # Adjust sample size as needed
      magic_min = magic,
      magic_max = NULL,
      x_min = total_min,
      x_max = bin_min
    )

    # Sample the upper out-of-bin range
    upper_samples <- frsr_sample(
      n = 1000,  # Adjust sample size as needed
      magic_min = magic,
      magic_max = NULL,
      x_min = bin_max,
      x_max = total_max
    )

    # Combine samples
    all_samples <- rbind(lower_samples, upper_samples)

    # Compute average error
    avg_error <- mean(all_samples$error)

    return(avg_error)
  }
  data %>%
    group_by(N) %>%
    mutate(
      OOB_Error = sapply(row_number(), function(i) {
        frsr_oob(Magic[i], Range_Min[i], Range_Max[i], 0.5, 2.0)
      })
    ) %>%
    ungroup()
}


## deepers

# this took hours to compute, so saving the result
# bin sizes: c(4, 6, 8, 12, 18, 24, 32, 48, 64, 96, 128)
varied_bins <- read.csv("../data/varied_bins.csv")

factor(varied_bins$Bottom,
       levels = c("point00781", "point0156",
                  "point0312", "point0625",
                  "point0125", "point25",
                  "point5", "one"),
       labels = c(0.0078, 0.0156, 0.0312, 0.0623, 0.125, 0.25, 0.5, 1),
       ordered = TRUE) -> varied_bins$Bottom


## 0.5 - 2.0 is a good range

too_many_bins <- frsr_bin_multi(sizes = c(8:24, 28,
                                          seq(32, 64, 2),
                                          seq(128, 256, 16),
                                          384),
                                x_min = 0.5)

too_many_bins <- read.csv("../data/too_many_bins.csv")


## Map magic and location/n combos

heatmap_data <- too_many_bins %>%
  select(Location, N, Magic) %>%
  distinct()

### this is fun!
ggplot(heatmap_data, aes(x = factor(Location), y = Magic, fill = factor(N))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_void() + guides(fill = "none")

### plots

## this is the stuff

ggplot(filter(full_bucket_plot) ) + geom_col(aes(x = factor(M), y = Width, group = factor(order), fill = factor(N)), color = "black", linewidth = 0.05) + guides(color = "none", fill = "none") + coord_flip() + theme_void() + scale_fill_discrete(type = sample(viridis::turbo(128)))

ggplot(filter(full_bucket_plot) ) + geom_col(aes(y = factor(M), x = Width, group = factor(order), fill = factor(N)), color = "black", linewidth = 0.05) + guides(color = "none", fill = "none") + scale_fill_discrete(type = sample(viridis::turbo(141)))

ggplot(full_bucket_plot, aes(ymin = as.numeric(factor(M)) - 0.5,
                             ymax = as.numeric(factor(M)) + 0.5,
                             xmin = Left,
                             xmax = Left + Width,
                             fill = factor(N))) +
  geom_rect(color = "black", linewidth = 0.1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 1.5)) +
  scale_y_continuous(breaks = 4:max(full_bucket_plot$M),
                     labels = unique(full_bucket_plot$M),
                     expand = c(0, 0), limits = c(12, NA)) +
  scale_fill_viridis_d() +
  theme_void() +
  labs(x = "Position", y = "M", fill = "N") + guides(fill = "none")

pair_plot <- function(bucket_df, spread_df) {
  # Scale the error values to match the M range
  y_range <- range(as.numeric(factor(bucket_df$M)))
  error_scaled <- scales::rescale(spread_df$error,
                                  to = y_range + c(-0.4, 0.4))


  p <- ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = factor(N))) +
    geom_point(data = spread_df,
              aes(x = input, y = error_scaled, color = cluster), shape = 16, size = 0.2, alpha = 0.4) +
    geom_smooth(data = spread_df, method = "lm",
                aes(x = input, y = error_scaled, color = cluster, group = cluster),
                formula =y ~ poly(x, 25)) +
    scale_fill_discrete(type = viridis::turbo(y_range[2])) +
    scale_x_continuous(limits = c(0.5, 2.0)) +
    scale_y_continuous(name = "Bucket (M)",
                       breaks = unique(as.numeric(factor(bucket_df$M))),
                       labels = unique(bucket_df$M)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + guides(fill = "none")

  return(p)
}

bp2_plot <- function(bucket_df) {
  # Scale the error values to match the M range
  y_range <- range(as.numeric(factor(bucket_df$M)))

  p <- ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = factor(N))) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + guides(fill = "none")
  return(p)
}

bp3_plot <- function(bucket_df) {
  # Scale the error values to match the M range
  y_range <- range(as.numeric(factor(bucket_df$M)))

  p <- ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = scale(Rarity) - scale(Depth) )) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + guides(fill = "none")
  return(p)
}

bucket_plot <- function(df, range = 4:64) {
  extent <- max(range) - min(range)

  df <- filter(bucket_selection(df, range = range))
  # Swap low with high values
  # there is probably a shorter way to do this
  d_sorted <- sort(df$Depth, index.return = TRUE)
  d_sorted$x <- rev(d_sorted$x)
  d_scale <- numeric()
  # d_scale should be low where Depth is high
  # and vice versa
  d_scale[d_sorted$ix] <- d_sorted$x

  plot_out <- ggplot(df, aes(ymin = as.numeric(factor(M)) - 0.5,
                 ymax = as.numeric(factor(M)) + 0.5,
                 xmin = Left, xmax = Left + Width,
                 fill = factor(N))) +
    geom_rect(linewidth = 0.2,
              color = "black") +
    scale_y_continuous(breaks = range,
                       labels = unique(df$M)) +
    scale_fill_discrete(type = sample(viridis::turbo(extent))) +
    theme_void() +
    labs(x = "Position", y = "M", fill = "N") +
    guides(fill = "none", color = "none") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = "white"))
  plot_out
}

bucket_plot_facets <- function(df, range = 4:64) {
  # Create a list of dataframes, each filtered by Bottom level
  bottom_levels <- levels(df$Bottom)
  df_list <- lapply(bottom_levels, function(level) {
    filtered_df <- filter(df, Bottom >= level)
    bucket_selection(filtered_df, range = range)  # Add this line
  })

  # Combine the dataframes and add a facet column
  combined_df <- bind_rows(df_list, .id = "facet")
  combined_df$facet <- factor(combined_df$facet,
                              labels = paste(">=", bottom_levels))

  # Create the plot
  plot_out <- ggplot(combined_df, aes(ymin = as.numeric(factor(M)) - 0.5,
                                      ymax = as.numeric(factor(M)) + 0.5,
                                      xmin = Left, xmax = Left + Width,
                                      fill = factor(N))) +
    geom_rect(color = "black", linewidth = 0.1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = range[1]:max(combined_df$M),
                       labels = unique(combined_df$M),
                       expand = c(0, 0)) +
    scale_fill_viridis_d() +
    theme_void() +
    labs(x = "Position", y = "M", fill = "N") +
    guides(fill = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          strip.text = element_text(size = 8)) +
    facet_wrap(~ facet, ncol = 1)

  plot_out
}



prepare_error_df <- function(df) {
  # Compute differences
  df$input_diff <- c(NA, diff(df$input))
  df$error_diff <- c(NA, diff(df$error))

  # Remove rows with NA values caused by diff()
  df_clean <- na.omit(df)

  # Return the cleaned data frame
  return(df_clean[, c("input", "error", "input_diff", "error_diff")])
}

## pseudo-waterfall

compare_n_values <- function(binned, n_small, n_large) {
  # Create a combined dataset with error values from both N
  combined_data <- binned %>%
    filter(N %in% c(n_small, n_large)) %>%
    mutate(N_type = if_else(N == n_small, "small", "large")) %>%
    pivot_wider(
      names_from = N_type,
      values_from = c(Range_Min, Range_Max, Avg_Error),
      id_cols = c(Magic)
    ) %>%
    # Add reference values by matching ranges
    filter(!is.na(Range_Min_small) & !is.na(Range_Min_large))

  ggplot() +
    # Horizontal segments for smaller N
    geom_segment(data = filter(binned, N == n_small),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Avg_Error, yend = Avg_Error),
                 color = "red", linewidth = 0.5) +
    # Horizontal segments for larger N
    geom_segment(data = filter(binned, N == n_large),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Avg_Error, yend = Avg_Error),
                 color = "blue", linewidth = 0.5) +
    # Rectangles showing difference
    geom_rect(data = filter(binned, N == n_large),
              aes(xmin = Range_Min, xmax = Range_Max,
                  ymin = Avg_Error,
                  ymax = filter(binned, N == n_small)$Avg_Error[
                    findInterval(Range_Min,
                                 filter(binned, N == n_small)$Range_Min)]),
              fill = "blue", alpha = 0.3) +
    labs(title = sprintf("Comparison of Avg Relative Error N=%d vs N=%d",
                         n_small, n_large),
         x = "Input Value",
         y = "Avg Relative Error")
}

compare_buckets <- function(bucket1, bucket2, error = "max") {
  error_col <- if(error == "max") "Max_Error" else "Avg_Error"

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
      cols = c(Avg_Error, Max_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) %>%
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Error", "Max_Error"),
                               labels = c("Avg Error", "Max Error")))

  # Create vertical segments dataset - now including both ends of each range
  segments_v <- bind_rows(
    # Segments for the ending ranges
    df %>%
      transmute(
        x = Range_Max,
        y_start = Avg_Error,
        y_end = Max_Error
      ),
    # Segments for the starting ranges
    df %>%
      transmute(
        x = Range_Min,
        y_start = Avg_Error,
        y_end = Max_Error
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
      cols = c(Avg_Error, Max_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) %>%
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Error", "Max_Error"),
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
        y_start = Avg_Error,
        y_end = Max_Error
      ),
    # Segments for the starting ranges
    binned %>%
      filter(N %in% n_values) %>%
      group_by(N) %>%
      slice(2:n()) %>%
      transmute(
        N = N,
        x = Range_Min,
        y_start = Avg_Error,
        y_end = Max_Error
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
    select(N, Max_Error, Width) %>%
    rename(Error = Max_Error)
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
