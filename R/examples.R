source("utils.R")

FLOAT_TOL <- 0.0004882812
SLICES <- 16384

# frsr_sample gives us log stratified samples
float_vector <- frsr_sample(n = SLICES * 2, x_min = 0.5, x_max = 2.0)$input

Blinn <- bind_cols(frsr(x = float_vector, NRmax = 10,
                     tol = FLOAT_TOL, magic = 0x5F400000,
                     A = 1.47, B = 0.47,
                     detail = TRUE, keep_params = TRUE),
                   method = "Blinn")

QuakeIII <- bind_cols(frsr(x = float_vector, NRmax = 10,
                        tol = FLOAT_TOL, magic = 0x5F3759DF,
                        detail = TRUE, keep_params = TRUE),
                      method = "QuakeIII")

Moroz <- bind_cols(frsr(x = float_vector, NRmax = 10,
                     tol = FLOAT_TOL, magic = 0x5F37ADD5,
                     detail = TRUE, keep_params = TRUE),
                    method = "Moroz")

Kahan <- bind_cols(frsr(x = float_vector, NRmax = 10,
                      tol = FLOAT_TOL, magic = 0x5f39d015,
                      detail = TRUE, keep_params = TRUE),
                    method = "Kahan")

approximated <- bind_rows(Blinn, QuakeIII, Moroz, Kahan)
rm(Blinn, Kahan, Moroz, QuakeIII)
## reference is generated here and not in C
approximated$reference <- with(approximated, 1/sqrt(input))


### Useful to compare Blinn with three tuned algorithms
### The shape of the error function of the other two are different becuase
### the impact of the constant is non-linear
### but the impact of a specifically good guess over Blinn's result is
### limited because of his specific NR parameter choice.
blinncomp_labels <- c(
  "Blinn" = "Blinn\n(1997)",
  "QuakeIII" = "Quake III\n(1999)",
  "Moroz" = "Moroz\n(2016)",
  "Kahan" = "Kahan\n(1999)")

approximated %>%
  ggplot(aes(x = input,
             y = (after_one - reference) / reference,
             color = method,
             linetype = method)) +
    geom_line(linewidth = 1.25) +
    ylab("Relative Error") +
    xlab("Input") +
    labs(color = "Algorithm", linetype = "Algorithm",
         title = "Performance of four Fast Inverse Square Root algorithms") +
    scale_color_manual(
      values = c(
        "Blinn" = "blue",
        "QuakeIII" = "green",
        "Moroz" = "red",
        "Kahan" = "orange"),
      labels = blinncomp_labels) +
    scale_linetype_manual(
      values = c("Blinn" = "dashed",
                 "QuakeIII" = "solid",
                 "Moroz" = "solid",
                 "Kahan" = "solid"),
      labels = blinncomp_labels
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
    theme(
      legend.key.size = unit(1.5, 'cm'),
      legend.text = element_text(margin = margin(l = 10, unit = "pt"))
    ) +
    scale_x_continuous(trans = 'log2',
                       breaks = scales::trans_breaks("log2", function(x) 2^x, n = 4),
                       labels = function(x) round(x, 4),
                       limits = c(0.25, 1))


## simple bar plot of iterations
approximated %>%
  ggplot() + geom_bar(aes(x = iters, fill = method), position = "dodge") +
  labs(title = "All four approximations converge to a tight tolerance in two iterations",
       y = "Count of samples",
       x = "Iterations",
       fill = "Approximation\nmethod")


### combines a more robust arrow system with the
### normal NR plot above
nrplot <- function(df = approximated, approx = "QuakeIII") {
  # First create the desired regularly spaced x values
  target_xs <- seq(0.25, 2, by = 0.125)
  temp <- df %>%
    filter(method %in% approx)
  # Then find the closest actual input values in the dataset
  segment_data <- temp %>%
    mutate(closest_target = target_xs[sapply(input, function(x)
      which.min(abs(target_xs - x)))]) %>%
    group_by(closest_target) %>%
    slice_min(abs(input - closest_target), n = 1) %>%
    ungroup()
  temp %>%
    mutate(relErrGuess = (initial - reference) / reference,
           relErrNR = (after_one - reference) / reference) %>%
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
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 show.legend = FALSE) +
    labs(y = "Relative error",
         title = "One iteration of Newton-Raphson markedly reduces error",
         x = "Input") -> output
  return(output)
}
nrplot()



