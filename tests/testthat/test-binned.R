find_optimal_buckets <- visualfrsr:::find_optimal_buckets
bucket_selection <- visualfrsr:::bucket_selection
analyze_bin_changes <- visualfrsr:::analyze_bin_changes
compute_oob_performance <- visualfrsr:::compute_oob_performance

test_that("find_optimal_buckets picks the known optimal slices", {
  domain_points <- sort(unique(c(tiny_bin_catalog$Range_Min, tiny_bin_catalog$Range_Max)))
  domain_midpoints <- (domain_points[-1] + domain_points[-length(domain_points)]) / 2

  m2_result <- find_optimal_buckets(tiny_bin_catalog, M = 2)
  expect_equal(nrow(m2_result), 2)
  expect_equal(m2_result$Location, c("A", "B"))
  expect_equal(m2_result$Max_Error, c(2, 4))
  coverage_counts_m2 <- sapply(domain_midpoints, function(mid) sum(m2_result$Range_Min <= mid & m2_result$Range_Max >= mid))
  expect_true(all(coverage_counts_m2 == 1))

  m3_result <- find_optimal_buckets(tiny_bin_catalog, M = 3)
  expect_equal(nrow(m3_result), 3)
  expect_equal(m3_result$Location, c("A", "G", "F"))
  expect_equal(m3_result$Max_Error, c(2, 3, 2))
  coverage_counts_m3 <- sapply(domain_midpoints, function(mid) sum(m3_result$Range_Min <= mid & m3_result$Range_Max >= mid))
  expect_true(all(coverage_counts_m3 == 1))
})


test_that("bucket_selection reports deterministic statistics for the fixture", {
  selection <- bucket_selection(tiny_bin_catalog, range = 2:3)
  expect_equal(unique(selection$M), 2:3)
  expect_equal(nrow(selection), 5)

  m2_rows <- selection[selection$M == 2, , drop = FALSE]
  expect_equal(m2_rows$Location, c("A", "B"))
  expect_equal(m2_rows$Width, c(1.5, 2.5))
  expect_equal(m2_rows$Midpoint, c(0.75, 2.75))
  expect_equal(m2_rows$order, c(1, 2))
  expect_equal(m2_rows$Depth, c(10 / 7.5, 5 / 7.5), tolerance = 1e-8)
  expect_equal(m2_rows$Variance, rep(0.5, 2), tolerance = 1e-8)

  m3_rows <- selection[selection$M == 3, , drop = FALSE]
  expect_equal(m3_rows$order, 1:3)
  expect_equal(m3_rows$Width, c(1.5, 1.5, 1.0))
  expect_equal(m3_rows$Midpoint, c(0.75, 2.25, 3.5))

  expected_rarity <- c(0.4, 0.2, 0.4, 0.2, 0.2)
  expect_equal(selection$Rarity, expected_rarity, tolerance = 1e-8)
  expect_equal(sum(selection$Rarity), sum(expected_rarity), tolerance = 1e-8)
})


test_that("analyze_bin_changes reports correct deltas on a simple input", {
  simple_df <- data.frame(
    input = c(1, 1.5, 2),
    error = c(0.1, 0.2, 0.4),
    magic = c(10, 12, 13)
  )

  deltas <- analyze_bin_changes(simple_df)
  expect_equal(nrow(deltas), 2)
  expect_equal(deltas$delta_error, c(0.1, 0.2), tolerance = 1e-8)
  expect_equal(deltas$delta_magic, c(2, 1), tolerance = 1e-8)
  expect_equal(deltas$mean_error, c(0.2, 0.4), tolerance = 1e-8)
  expect_equal(deltas$mean_magic, c(12, 13), tolerance = 1e-8)
})

test_that("compute_oob_performance maintains structure for a tiny catalog", {
  tiny_oob_input <- data.frame(
    Magic = c(1.0, 1.3),
    Range_Min = c(0.6, 1.3),
    Range_Max = c(1.2, 1.8),
    N = c(10, 5),
    stringsAsFactors = FALSE
  )

  oob_result <- compute_oob_performance(tiny_oob_input)
  expect_s3_class(oob_result, "data.frame")
  expect_equal(nrow(oob_result), nrow(tiny_oob_input))
  expect_true(all(c("Magic", "Range_Min", "Range_Max", "N", "OOB_Error") %in% names(oob_result)))
  expect_type(oob_result$OOB_Error, "double")
})
