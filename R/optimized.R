source("utils.R")


# Compute results for each input and parameter combination
compute_result_block <- function(slices = 100000, GRID_SIZE = 50) {
  half_extent <- 0.15625
  inputs <- frsr_sample(slices, x_min = 2^-5, x_max = 1.0)$input
  # Initialize "grids" for A and B parameters
  Ain <- runif(slices, min = 1.5 - half_extent, max = 1.5 + half_extent)
  Bin <- runif(slices, min = 0.5 - half_extent, max = 0.5 + half_extent)
  frsr(x = inputs, magic = 0x5F400000,
       NR = 1, A = Ain, B = Bin,
       detail = TRUE, keep_params = TRUE) -> results
  results <- results[order(results[, "input"]), ]
  results$A_rank <- ntile(results$A, GRID_SIZE)
  results$B_rank <- ntile(results$B, GRID_SIZE)
  results$err_rank <- ntile(results$error, GRID_SIZE)
  results$pair <- paste0("(", round(results$A, 4), ", ", round(results$B, 4), ")")
  return(results)
}

# Compute results
optimized <- compute_result_block()


## tiling plot
ggplot(optimized, aes(x = A_rank, y = B_rank, fill = err_rank)) + geom_tile()


## good artistic plot of errors by grid location
ggplot(optimized, aes(x = input,
                      y = error,
                      color = pair)) +
  geom_path() +
  guides(color = "none") +
  coord_polar(theta = "x") +
  theme_void()

# gets the point across
ggplot(optimized, aes(x = A, y = B, color = error)) + geom_point()
