source("utils.R")


## helper functions
cluster_magic_constants <- function(data, n_clusters = 3) {
  # Group by 'magic' and calculate statistics for 'error'
  grouped <- data %>%
    group_by(magic) %>%
    summarize(
      mean_error = mean(error),
      min_error = min(error),
      max_error = max(error)
    )
  
  # Perform k-means clustering
  kmeans_result <- kmeans(grouped[, c("mean_error", "min_error", "max_error")], centers = n_clusters)
  # Add cluster assignments to the grouped data
  grouped$cluster <- kmeans_result$cluster
  
  # Calculate cluster statistics
  cluster_stats <- grouped %>%
    group_by(cluster) %>%
    summarize(
      average_magic = mean(magic),
      min_magic = min(magic),
      max_magic = max(magic),
      count = n()
    )
  
  # Prepare the result as a tidy data frame
  result <- cluster_stats %>%
    pivot_longer(cols = c(average_magic, min_magic, max_magic),
                 names_to = "statistic",
                 values_to = "value") %>%
    mutate(statistic = factor(statistic, 
                              levels = c("average_magic", "min_magic", "max_magic"),
                              labels = c("Average", "Minimum", "Maximum")))
  
  return(result)
}



### some scratch

low_err_05_2 <- frsr_sample(n_floats = 2048, NRmax = 1, keep_params = TRUE,
                           x_min = 0.5, x_max = 2.0,
                           n_ints = 16384,
                           magic_min = 1596980000L,
                           magic_max = 1598050000L) %>%
  select(input, error, magic)


filter(cluster_magic_constants(low_err_05_2), cluster %in% c(1, 2)) -> low_error_clusters



wide_05_2 <- frsr_sample(n_floats = 2048, NRmax = 1, keep_params = TRUE,
                            x_min = 0.5, x_max = 2.0,
                            n_ints = 16384,
                            magic_min = 1596910000L,
                            magic_max = 1598100000L) %>%
  select(input, error, magic)

filter(cluster_magic_constants(wide_05_2), cluster == 1) -> high_err_cluster
high_err_cluster$cluster <- 3

total_clust <- rbind(low_error_clusters, high_err_cluster) %>%
  filter(statistic != "Average") %>% select(!count)



# these can be used to make "bands" of probably representative
# error structures

# clusters
# cluster   count  statistic         value
# --------  ------  ----------  -----------
#   1         6901  Average      1597398039
#   1         6901  Minimum      1597176038
#   1         6901  Maximum      1597622357
#   2         2007  Average      1597984798
#   2         2007  Minimum      1597918602
#   2         2007  Maximum      1598049921
#   3         5144  Average      1597435684
#   3         5144  Minimum      1596910052
#   3         5144  Maximum      1597884263

sample_clusters <- function(n = 16384) {
  c(1597622357, 1597176038) -> cl1
  c(1598049921, 1597918602) -> cl2
  c(1597884263, 1596910052) -> cl3
  cln <- rbind(cl1, cl2, cl3)
  
  representative_clusters <- data.frame()
  for (i in 1:3) {
    fast_sample <- frsr_sample(n = n, NRmax = 1, keep_params = TRUE,
                               x_min = 0.5, x_max = 2.0,
                               magic_max = cln[i,1], magic_min = cln[i,2])
    fast_sample %>% select(input, error, magic) %>% 
      mutate(cluster = i) -> fast_sample
      
    representative_clusters <- rbind(representative_clusters, fast_sample)
  }
  
  return(representative_clusters %>% mutate(cluster = factor(cluster)))
}

ggplot(sample_clusters(),
       aes(x = input, y = error, color = cluster, group = cluster)) +
  geom_point(shape = ".") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 25), se = FALSE)


