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

    bucket_df |>
      mutate(
        Depth = N / mean_N,
        Width = Range_Max - Range_Min,
        Midpoint = (Range_Max + Range_Min) / 2,
        Left = Range_Min,
        Variance = ((N - mean_N) / sd_N)^2,
        M = M,
        Error = Max_Error
      ) |>
      select(M, Location, N, Depth, Error, Variance, Width, Midpoint, Left)
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

  output |>
    group_by(M) |>
    arrange(M, Midpoint) |>
    mutate(order = row_number()) |>
    ungroup() |>
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
  bin_stats <- df |>
    group_by(bin) |>
    summarize(
      mean_error = mean(error),
      mean_magic = mean(magic),
      .groups = 'drop'
    ) |>
    # Calculate changes between consecutive bins
    mutate(
      delta_error = mean_error - lag(mean_error),
      delta_magic = mean_magic - lag(mean_magic)
    ) |>
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
  data |>
    group_by(N) |>
    mutate(
      OOB_Error = sapply(row_number(), function(i) {
        frsr_oob(Magic[i], Range_Min[i], Range_Max[i], 0.5, 2.0)
      })
    ) |>
    ungroup()
}

load_binned_csv <- function(filename) {
  path <- system.file("extdata", filename, package = "visualfrsr")
  if (path == "") {
    stop(sprintf("Cannot find %s in inst/extdata.", filename), call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

#' Load the varied bin dataset with ordered Bottom factor.
#'
#' @return Data frame.
load_varied_bins <- function() {
  df <- load_binned_csv("varied_bins.csv")
  factor(df$Bottom,
         levels = c(
           "point00781", "point0156", "point0312", "point0625",
           "point0125", "point25", "point5", "one"
         ),
         labels = c(0.0078, 0.0156, 0.0312, 0.0623, 0.125, 0.25, 0.5, 1),
         ordered = TRUE) -> df$Bottom
  df
}

#' Load precomputed bucket sweeps.
#'
#' @param filename Name of the CSV stored under `inst/extdata`.
#' @return Data frame.
load_bucket_sweep <- function(filename = "too_many_bins.csv") {
  load_binned_csv(filename)
}

#' Heatmap-ready data summarizing location, bin count, and magic.
#'
#' @param bucket_df Data frame returned by [load_bucket_sweep()].
#' @return Tidy data frame.
build_heatmap_data <- function(bucket_df = load_bucket_sweep()) {
  bucket_df |>
    select(Location, N, Magic) |>
    distinct()
}

## deepers

# this took hours to compute, so saving the result
# bin sizes: c(4, 6, 8, 12, 18, 24, 32, 48, 64, 96, 128)
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

## "weight" error by the fraction of the domain
## it covers.
norm_errorN <- function(df, bins) {
  bucket <- find_optimal_buckets(df, bins)
  bucket$Width <- (bucket[, "Range_Max"] - bucket[, "Range_Min"]) / 0.75
  bucket <- bucket |>
    # if you don't pick N, dplyr complains and does it anyway
    select(N, Max_Error, Width) |>
    rename(Error = Max_Error)
  ## avoids privileging small bucket sizes
  sum_err <- with(bucket, sum(Error*Width))
  return(sum_err)
}
