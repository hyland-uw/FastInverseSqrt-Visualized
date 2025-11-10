#' Generate deterministic comparison samples across classical FRSR variants.
#'
#' @param slices Base number of slices for the log-stratified samples.
#' @param float_tol Tolerance passed to [frsrr::frsr()].
#' @return Tidy data frame combining four approximation families.
generate_example_samples <- function(
    slices = 16384,
    float_tol = 0.0004882812) {
  float_vector <- frsr_sample(n = slices * 2, x_min = 0.5, x_max = 2.0)$input

  build_row <- function(magic, method, A = 1.5, B = 0.5) {
    bind_cols(
      frsr(
        x = float_vector,
        NRmax = 10,
        tol = float_tol,
        magic = magic,
        A = A,
        B = B,
        detail = TRUE,
        keep_params = TRUE
      ),
      method = method
    )
  }

  approximated <- bind_rows(
    build_row(0x5F400000, "Blinn", A = 1.47, B = 0.47),
    build_row(0x5F3759DF, "QuakeIII"),
    build_row(0x5F37ADD5, "Moroz"),
    build_row(0x5f39d015, "Kahan")
  )

  approximated$reference <- with(approximated, 1 / sqrt(input))
  approximated
}
