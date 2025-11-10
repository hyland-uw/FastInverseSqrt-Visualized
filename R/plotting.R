# Internal plotting helpers ------------------------------------------------

.deconstruction_palettes <- function(df) {
  iter_levels <- as.integer(levels(df$iters))
  iter_colors <- colorRampPalette(c("dodgerblue2", "red"))(length(iter_levels))
  iter_rank_hue <- c("lightblue", colorRampPalette(c("white", "orange1", "red"))(7))

  list(
    iter_levels = rev(iter_levels),
    iter_colors = setNames(iter_colors, seq_along(iter_colors)),
    iter_rank_hue = iter_rank_hue
  )
}

.mc_annotate <- function(magic_value, label,
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

.create_geom_points <- function(data, iter_range, shape, size, alpha = 1) {
  lapply(iter_range, function(i) {
    geom_point(data = data[data$iters == i, ],
               shape = shape,
               size = size,
               alpha = alpha)
  })
}

.blinncomp_labels <- c(
  "Blinn" = "Blinn\n(1997)",
  "QuakeIII" = "Quake III\n(1999)",
  "Moroz" = "Moroz\n(2016)",
  "Kahan" = "Kahan\n(1999)"
)

# Deconstruction plots -----------------------------------------------------

deconstruction_zoom_plot <- function(narrowed) {
  ggplot(narrowed, aes(
    x = (initial - reference) / reference,
    y = magic
  )) +
    geom_point(shape = ".", alpha = 0.9) +
    xlab("Relative Error") + ylab("Magic Constant") +
    labs(color = "Iterations\nto converge",
         title = "Zooming in on three similar constants") +
    scale_y_continuous(
      labels = function(x) sprintf("0x%X", as.integer(x)),
      limits = c(0x5f37642f - 3200, 0x5F376D60)
    ) +
    .mc_annotate(0x5f3759df, "0x5f3759df", "blue", x_end = 0.04) +
    .mc_annotate(0x5f37642f, "0x5f37642f", "red", x_end = 0.043) +
    .mc_annotate(0x5f375a86, "0x5f375a86", "orange") +
    .mc_annotate(0x5F376908, "0x5F376908", "purple") +
    xlim(-0.035, 0.08)
}

deconstruction_rate_plot <- function(deconstructed) {
  pals <- .deconstruction_palettes(deconstructed)
  ggplot(
    deconstructed,
    aes(
      x = input,
      y = initial - reference,
      color = iters
    )
  ) +
    .create_geom_points(deconstructed, pals$iter_levels, 16, 2, 0.95) +
    guides(alpha = "none") +
    scale_color_manual(values = pals$iter_colors, breaks = 1:6) +
    labs(
      color = "Iteration\nCount",
      title = "Rate of convergence is not symmetric about first guess errors"
    ) +
    ylab("Error before NR Iteration") +
    xlab("Input float")
}

deconstruction_magic_span_plot <- function(deconstructed) {
  ggplot(
    deconstructed |>
      filter(as.numeric(iters) <= 4) |>
      group_by(iters) |>
      summarize(min_magic = min(magic), max_magic = max(magic), .groups = "drop"),
    aes(x = iters, ymin = min_magic, ymax = max_magic)
  ) +
    scale_y_continuous(labels = function(x) sprintf("0x%X", as.integer(x))) +
    geom_errorbar(width = 0.5) +
    geom_point(aes(y = min_magic), color = "blue") +
    geom_point(aes(y = max_magic), color = "red") +
    labs(
      x = "Iterations until convergence",
      y = "Integer value",
      title = "Good constants exist only in a narrow range"
    ) +
    theme_minimal()
}

deconstruction_quadratic_plot <- function(deconstructed) {
  pals <- .deconstruction_palettes(deconstructed)
  ggplot(
    deconstructed,
    aes(
      x = (initial - reference) / reference,
      y = (after_one - reference) / reference,
      color = iter_rank
    )
  ) +
    geom_point(shape = 16, size = 0.8, alpha = 0.9) +
    scale_color_manual(
      values = pals$iter_rank_hue,
      guide = guide_legend(override.aes = list(size = 1.5))
    ) +
    labs(
      title = "NR converges quadratically",
      x = "Relative error before Newton-Raphson",
      y = "Relative error after one iteration",
      color = "Iterations\nuntil\neventual\nconvergence"
    )
}

deconstruction_improvement_plot <- function(deconstructed) {
  pals <- .deconstruction_palettes(deconstructed)
  ggplot(
    deconstructed,
    aes(
      x = (initial - reference) / reference,
      y = abs(initial - reference) / reference -
        abs(after_one - reference) / reference,
      color = iter_rank
    )
  ) +
    geom_point(shape = 16, size = 0.5) +
    scale_color_manual(
      values = pals$iter_rank_hue,
      guide = guide_legend(override.aes = list(size = 3))
    ) +
    labs(
      x = "Relative error of first guess",
      y = "Improvement from one Newton-Raphson step",
      color = "Iterations\nto convergence",
      title = "Plotted against relative improvement, optimal region is visible"
    )
}

deconstruction_painterly_plot <- function(deconstructed) {
  ggplot(
    deconstructed,
    aes(
      x = magic,
      y = input,
      color = iters
    )
  ) +
    .create_geom_points(deconstructed, 6:1, 16, 6) +
    guides(alpha = "none", color = "none", shape = "none", size = "none") +
    theme_void()
}

deconstruction_polar_plot <- function(deconstructed) {
  ggplot(
    deconstructed,
    aes(
      x = initial - reference,
      y = after_one,
      color = iters
    )
  ) +
    geom_point(shape = ".") +
    guides(color = "none") +
    coord_polar(theta = "x") +
    theme_void()
}

deconstruction_combined_plot <- function(deconstructed, widened) {
  pals <- .deconstruction_palettes(deconstructed)

  subset_plot <- deconstructed |>
    filter(iters %in% levels(iters)[1:5]) |>
    ggplot(aes(
      x = (initial - reference) / reference,
      y = magic, color = iters
    )) +
    geom_point(shape = 16, size = 0.65, alpha = 0.95) +
    labs(color = "Iterations\nto converge",
         title = "Shaded region is 0.024% of the 32 bit integers") +
    guides(colour = "none") +
    scale_color_manual(values = pals$iter_rank_hue[1:5]) +
    scale_x_continuous(breaks = c(0), limits = c(-0.25, 0.25)) +
    ylim(1.5935e9, 1601175552) +
    theme(
      plot.title = element_text(hjust = 0.45),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    annotate("rect", alpha = 0.2, fill = "blue",
             xmin = -Inf, xmax = -0.125,
             ymin = 0x5F300000, ymax = 0x5F400000) +
    annotate("rect", alpha = 0.2, fill = "blue",
             xmin = 0.125, xmax = Inf,
             ymin = 0x5F300000, ymax = 0x5F400000) +
    geom_hline(yintercept = 0x5F400000,
               color = "blue", alpha = 0.6, lty = 4, linewidth = 0.2) +
    geom_hline(yintercept = 0x5F300000,
               color = "blue", alpha = 0.6, lty = 4, linewidth = 0.2)

  ggplot(
    widened,
    aes(
      x = (initial - reference) / reference,
      y = magic,
      color = iter_rank
    )
  ) +
    geom_point(shape = 16, size = 0.65, alpha = 0.95, show.legend = TRUE) +
    xlab("Relative error") + ylab("Restoring constant (in billions)") +
    labs(
      color = "Iterations\nto converge",
      title = "The region where the approximation is optimal is tiny.",
      fill = "Range of\noptimal integers"
    ) +
    scale_y_continuous(
      labels = function(x) sprintf("%.1f", x / 1e9),
      limits = c(1.3e9 - 1, NA)
    ) +
    scale_color_manual(values = pals$iter_rank_hue, na.translate = FALSE, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    annotate("rect", alpha = 0.2, fill = "blue",
             xmin = -Inf, xmax = -0.1,
             ymin = 0x5F300000, ymax = 0x5F400000) +
    annotate("rect", alpha = 0.2, fill = "blue",
             xmin = 0.1, xmax = Inf,
             ymin = 0x5F300000, ymax = 0x5F400000) +
    annotate("rect",
             fill = NA, color = "black",
             xmin = -0.5, xmax = 0.5,
             ymin = 1.325e9, ymax = 1.525e9) +
    annotate("segment", x = -0.15, xend = -0.5, y = 1602500000, yend = 1.525e9,
             linetype = "dashed") +
    annotate("segment", x = 0.15, xend = 0.5, y = 1602500000, yend = 1.525e9,
             linetype = "dashed") +
    annotate("rect", color = "black", fill = NA, linetype = "dashed",
             xmin = -0.15, xmax = 0.15,
             ymin = 1592500000, ymax = 1602500000) +
    annotate("segment", x = 0.15, xend = 0.5, y = 1592500000, yend = 1.325e9,
             linetype = "dashed") +
    annotate("segment", x = -0.15, xend = -0.5, y = 1592500000, yend = 1.325e9,
             linetype = "dashed") +
    geom_point(aes(fill = "0x5F300000 to\n0x5F400000"), alpha = 0) +
    scale_fill_manual(values = c("0x5F300000 to\n0x5F400000" = "blue"),
                      name = "Optimal\ninteger values") +
    guides(fill = guide_legend(override.aes = list(
      color = "blue",
      alpha = 0.2,
      size = 5,
      shape = 15
    ))) +
    annotation_custom(
      grob = ggplot2::ggplotGrob(subset_plot),
      xmin = -0.5, xmax = 0.5,
      ymin = 1.325e9, ymax = 1.525e9
    )
}

# Historical example plots -------------------------------------------------

plot_example_comparison <- function(approximated = generate_example_samples()) {
  ggplot(
    approximated,
    aes(
      x = input,
      y = (after_one - reference) / reference,
      color = method,
      linetype = method
    )
  ) +
    geom_line(linewidth = 1.25) +
    ylab("Relative Error") +
    xlab("Input") +
    labs(
      color = "Algorithm",
      linetype = "Algorithm",
      title = "Performance of four Fast Inverse Square Root algorithms"
    ) +
    scale_color_manual(
      values = c(
        "Blinn" = "blue",
        "QuakeIII" = "green",
        "Moroz" = "red",
        "Kahan" = "orange"
      ),
      labels = .blinncomp_labels
    ) +
    scale_linetype_manual(
      values = c(
        "Blinn" = "dashed",
        "QuakeIII" = "solid",
        "Moroz" = "solid",
        "Kahan" = "solid"
      ),
      labels = .blinncomp_labels
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
    theme(
      legend.key.size = grid::unit(1.5, "cm"),
      legend.text = element_text(margin = margin(l = 10, unit = "pt"))
    ) +
    scale_x_continuous(
      trans = "log2",
      breaks = trans_breaks("log2", function(x) 2^x, n = 4),
      labels = function(x) round(x, 4),
      limits = c(0.25, 1)
    )
}

plot_iteration_histogram <- function(approximated = generate_example_samples()) {
  ggplot(approximated) +
    geom_bar(aes(x = iters, fill = method), position = "dodge") +
    labs(
      title = "All four approximations converge to a tight tolerance in two iterations",
      y = "Count of samples",
      x = "Iterations",
      fill = "Approximation\nmethod"
    )
}

nrplot <- function(df = generate_example_samples(), approx = "QuakeIII") {
  target_xs <- seq(0.25, 2, by = 0.125)
  temp <- df |>
    filter(method %in% approx)
  segment_data <- temp |>
    mutate(closest_target = target_xs[sapply(input, function(x)
      which.min(abs(target_xs - x)))]) |>
    group_by(closest_target) |>
    slice_min(abs(input - closest_target), n = 1) |>
    ungroup()
  temp |>
    mutate(relErrGuess = (initial - reference) / reference,
           relErrNR = (after_one - reference) / reference) |>
    ggplot(aes(x = input)) +
    geom_ribbon(aes(ymin = pmin(relErrGuess, relErrNR),
                    ymax = pmax(relErrGuess, relErrNR),
                    fill = method),
                alpha = 0.3) +
    geom_line(aes(y = relErrGuess,
                  linetype = "Initial guess"),
              linewidth = 1.1) +
    geom_line(aes(y = relErrNR,
                  linetype = "After one iteration")) +
    scale_linetype_manual(name = "Newton-Raphson",
                          breaks = c("Initial guess", "After one iteration"),
                          values = c("Initial guess" = "solid",
                                     "After one iteration" = "dotted")) +
    guides(fill = "none") +
    geom_segment(data = segment_data,
                 aes(x = input,
                     xend = input,
                     y = (initial - reference) / reference,
                     yend = (after_one - reference) / reference),
                 linewidth = 0.5,
                 arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed"),
                 show.legend = FALSE) +
    labs(
      y = "Relative error",
      title = "One iteration of Newton-Raphson markedly reduces error",
      x = "Input"
    )
}

# Binned pipeline plots ----------------------------------------------------

bucket_heatmap_plot <- function(bucket_df = load_bucket_sweep()) {
  ggplot(build_heatmap_data(bucket_df),
         aes(x = factor(Location), y = Magic, fill = factor(N))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_void() + guides(fill = "none")
}

plot_bucket_widths <- function(selection_df) {
  ggplot(selection_df,
         aes(x = factor(M),
             y = Width,
             group = factor(order),
             fill = factor(N))) +
    geom_col(color = "black", linewidth = 0.05) +
    guides(color = "none", fill = "none") +
    coord_flip() +
    theme_void() +
    scale_fill_discrete(type = sample(turbo(length(unique(selection_df$N)))))
}

plot_bucket_widths_horizontal <- function(selection_df) {
  ggplot(selection_df,
         aes(y = factor(M),
             x = Width,
             group = factor(order),
             fill = factor(N))) +
    geom_col(color = "black", linewidth = 0.05) +
    guides(color = "none", fill = "none") +
    scale_fill_discrete(type = sample(turbo(length(unique(selection_df$N)))))
}

plot_bucket_rectangles <- function(selection_df) {
  ggplot(selection_df,
         aes(ymin = as.numeric(factor(M)) - 0.5,
             ymax = as.numeric(factor(M)) + 0.5,
             xmin = Left,
             xmax = Left + Width,
             fill = factor(N))) +
    geom_rect(color = "black", linewidth = 0.1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, 1.5)) +
    scale_y_continuous(
      breaks = 4:max(selection_df$M),
      labels = unique(selection_df$M),
      expand = c(0, 0),
      limits = c(12, NA)
    ) +
    viridis::scale_fill_viridis_d() +
    theme_void() +
    labs(x = "Position", y = "M", fill = "N") +
    guides(fill = "none")
}

pair_plot <- function(bucket_df, spread_df) {
  y_range <- range(as.numeric(factor(bucket_df$M)))
  error_scaled <- rescale(spread_df$error,
                                  to = y_range + c(-0.4, 0.4))

  ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = factor(N))) +
    geom_point(data = spread_df,
               aes(x = input, y = error_scaled, color = cluster),
               shape = 16, size = 0.2, alpha = 0.4) +
    geom_smooth(data = spread_df, method = "lm",
                aes(x = input, y = error_scaled, color = cluster, group = cluster),
                formula = y ~ poly(x, 25)) +
    scale_fill_discrete(type = turbo(y_range[2])) +
    scale_x_continuous(limits = c(0.5, 2.0)) +
    scale_y_continuous(name = "Bucket (M)",
                       breaks = unique(as.numeric(factor(bucket_df$M))),
                       labels = unique(bucket_df$M)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + guides(fill = "none")
}

bp2_plot <- function(bucket_df) {
  ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = factor(N))) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + guides(fill = "none")
}

bp3_plot <- function(bucket_df) {
  ggplot() +
    geom_rect(data = bucket_df,
              aes(ymin = as.numeric(factor(M)) - 0.5,
                  ymax = as.numeric(factor(M)) + 0.5,
                  xmin = Left, xmax = Left + Width,
                  fill = scale(Rarity) - scale(Depth))) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + guides(fill = "none")
}

bucket_plot <- function(df, range = 4:64) {
  extent <- max(range) - min(range)

  df <- filter(bucket_selection(df, range = range))
  d_sorted <- sort(df$Depth, index.return = TRUE)
  d_sorted$x <- rev(d_sorted$x)
  d_scale <- numeric()
  d_scale[d_sorted$ix] <- d_sorted$x

  ggplot(df, aes(ymin = as.numeric(factor(M)) - 0.5,
                 ymax = as.numeric(factor(M)) + 0.5,
                 xmin = Left, xmax = Left + Width,
                 fill = factor(N))) +
    geom_rect(linewidth = 0.2,
              color = "black") +
    scale_y_continuous(breaks = range,
                       labels = unique(df$M)) +
    scale_fill_discrete(type = sample(turbo(extent))) +
    theme_void() +
    labs(x = "Position", y = "M", fill = "N") +
    guides(fill = "none", color = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"))
}

bucket_plot_facets <- function(df, range = 4:64) {
  bottom_levels <- levels(df$Bottom)
  df_list <- lapply(bottom_levels, function(level) {
    filtered_df <- filter(df, Bottom >= level)
    bucket_selection(filtered_df, range = range)
  })

  combined_df <- bind_rows(df_list, .id = "facet")
  combined_df$facet <- factor(combined_df$facet,
                              labels = paste(">=", bottom_levels))

  ggplot(combined_df, aes(ymin = as.numeric(factor(M)) - 0.5,
                          ymax = as.numeric(factor(M)) + 0.5,
                          xmin = Left, xmax = Left + Width,
                          fill = factor(N))) +
    geom_rect(color = "black", linewidth = 0.1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = range[1]:max(combined_df$M),
                       labels = unique(combined_df$M),
                       expand = c(0, 0)) +
    viridis::scale_fill_viridis_d() +
    theme_void() +
    labs(x = "Position", y = "M", fill = "N") +
    guides(fill = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          strip.text = element_text(size = 8)) +
    facet_wrap(~ facet, ncol = 1)
}

compare_n_values <- function(binned, n_small, n_large) {
  combined_data <- binned |>
    filter(N %in% c(n_small, n_large)) |>
    mutate(N_type = if_else(N == n_small, "small", "large")) |>
    pivot_wider(
      names_from = N_type,
      values_from = c(Range_Min, Range_Max, Avg_Error),
      id_cols = c(Magic)
    ) |>
    filter(!is.na(Range_Min_small) & !is.na(Range_Min_large))

  ggplot() +
    geom_segment(data = filter(binned, N == n_small),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Avg_Error, yend = Avg_Error),
                 color = "red", linewidth = 0.5) +
    geom_segment(data = filter(binned, N == n_large),
                 aes(x = Range_Min, xend = Range_Max,
                     y = Avg_Error, yend = Avg_Error),
                 color = "blue", linewidth = 0.5) +
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
  error_col <- if (error == "max") "Max_Error" else "Avg_Error"

  x_points <- sort(unique(c(
    bucket1$Range_Min, bucket1$Range_Max,
    bucket2$Range_Min, bucket2$Range_Max
  )))

  rectangles <- tibble(
    xmin = x_points[-length(x_points)],
    xmax = x_points[-1]
  ) |>
    mutate(
      error1 = sapply(xmin, function(x) {
        bucket1[[error_col]][x >= bucket1$Range_Min & x < bucket1$Range_Max][1]
      }),
      error2 = sapply(xmin, function(x) {
        bucket2[[error_col]][x >= bucket2$Range_Min & x < bucket2$Range_Max][1]
      }),
      ymin = pmin(error1, error2),
      ymax = pmax(error1, error2),
      fill_color = if_else(error2 < error1, "blue", "red")
    )

  segments_h1 <- bucket1 |>
    select(Range_Min, Range_Max, Error = !!sym(error_col))

  segments_h2 <- bucket2 |>
    select(Range_Min, Range_Max, Error = !!sym(error_col))

  ggplot() +
    geom_segment(data = segments_h1,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error),
                 color = "black", linewidth = 0.5) +
    geom_segment(data = segments_h2,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error),
                 color = "black", linewidth = 0.5) +
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

plot_bucket <- function(df) {
  segments_h <- df |>
    pivot_longer(
      cols = c(Avg_Error, Max_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) |>
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Error", "Max_Error"),
                               labels = c("Avg Error", "Max Error")))

  segments_v <- bind_rows(
    df |>
      transmute(
        x = Range_Max,
        y_start = Avg_Error,
        y_end = Max_Error
      ),
    df |>
      transmute(
        x = Range_Min,
        y_start = Avg_Error,
        y_end = Max_Error
      )
  )
  ggplot() +
    geom_segment(data = segments_h,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error,
                     color = Error_Type),
                 linewidth = 0.5) +
    geom_segment(data = segments_v,
                 aes(x = x, xend = x,
                     y = y_start, y_end = y_end),
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

plot_multiple_n <- function(binned, n_values = unique(binned$N)) {
  segments_h <- binned |>
    filter(N %in% n_values) |>
    pivot_longer(
      cols = c(Avg_Error, Max_Error),
      names_to = "Error_Type",
      values_to = "Error"
    ) |>
    mutate(Error_Type = factor(Error_Type,
                               levels = c("Avg_Error", "Max_Error"),
                               labels = c("Avg Error", "Max Error")))

  segments_v <- bind_rows(
    binned |>
      filter(N %in% n_values) |>
      group_by(N) |>
      slice(1:(n()-1)) |>
      transmute(
        N = N,
        x = Range_Max,
        y_start = Avg_Error,
        y_end = Max_Error
      ),
    binned |>
      filter(N %in% n_values) |>
      group_by(N) |>
      slice(2:n()) |>
      transmute(
        N = N,
        x = Range_Min,
        y_start = Avg_Error,
        y_end = Max_Error
      )
  )

  ggplot() +
    geom_segment(data = segments_h,
                 aes(x = Range_Min, xend = Range_Max,
                     y = Error, yend = Error,
                     color = Error_Type),
                 linewidth = 0.5) +
    geom_segment(data = segments_v,
                 aes(x = x, xend = x,
                     y = y_start, y_end = y_end),
                 linetype = "dotted",
                 color = "black",
                 linewidth = 0.25,
                 alpha = 0.75) +
    facet_grid(cols = vars(N)) +
    scale_color_manual(values = c("Avg Error" = "blue", "Max Error" = "red")) +
    labs(title = "Error Trends Across Input Range",
         x = "Input Value",
         y = "Error") +
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25),
                       limits = c(0.25, 1.0)) +
    theme_minimal()
}

plot_normalized_error_curve <- function(binned, bins = 4:36) {
  map_dfr(bins, ~tibble(bins = .x, error = norm_errorN(binned, .x))) |>
    ggplot(aes(x = bins, y = error)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Bins",
      y = "Normalized Error",
      title = "Optimal bucket selection error reduction slows after 24 bins"
    ) +
    scale_x_continuous(breaks = seq(min(bins), max(bins), by = 4))
}

# Clustering plots ---------------------------------------------------------

plot_sampled_clusters <- function(clusters = sample_clusters()) {
  ggplot(
    clusters,
    aes(x = input, y = error, color = cluster, group = cluster)
  ) +
    geom_point(shape = ".") +
    geom_smooth(method = "lm", formula = y ~ poly(x, 25), se = FALSE) +
    labs(
      title = "Representative cluster bands for FRSR magic constants",
      x = "Input value",
      y = "Relative error",
      color = "Cluster"
    )
}

# Optimization plots -------------------------------------------------------

plot_optimized_tile <- function(optimized = compute_result_block()) {
  ggplot(optimized, aes(x = A_rank, y = B_rank, fill = err_rank)) +
    geom_tile() +
    labs(
      title = "Heatmap of NR parameter errors",
      x = "A rank",
      y = "B rank",
      fill = "Error rank"
    )
}

plot_optimized_paths <- function(optimized = compute_result_block()) {
  ggplot(
    optimized,
    aes(
      x = input,
      y = error,
      color = pair
    )
  ) +
    geom_path() +
    guides(color = "none") +
    coord_polar(theta = "x") +
    theme_void()
}

plot_optimized_pairs <- function(optimized = compute_result_block()) {
  ggplot(optimized, aes(x = A, y = B, color = error)) +
    geom_point() +
    labs(
      title = "Parameter search landscape",
      x = "Half-three parameter",
      y = "Half-one parameter",
      color = "Error"
    )
}
