# Helper fixtures for deterministic bucket tests

tiny_bin_catalog <- data.frame(
  Range_Min = c(0, 1.5, 0, 2, 0, 3, 1.5),
  Range_Max = c(1.5, 4, 2, 4, 3, 4, 3),
  Max_Error = c(2, 4, 6, 6, 7, 2, 3),
  N = c(10, 5, 3, 4, 6, 2, 8),
  Location = c("A", "B", "C", "D", "E", "F", "G"),
  Magic = c(1.0, 1.2, 1.1, 1.3, 1.4, 1.5, 1.6),
  stringsAsFactors = FALSE
)
