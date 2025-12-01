# tests/testthat/test-annotation.R

library(testthat)
library(dplyr)
library(readr)

test_that("Annotation candidate building works correctly", {
  
  # Load test data using stable project-relative paths
  vip_df     <- read_csv("tests/testdata/sample_vip.csv",     show_col_types = FALSE)
  limma_df   <- read_csv("tests/testdata/sample_limma.csv",   show_col_types = FALSE)
  featdef_df <- read_csv("tests/testdata/sample_featdef.csv", show_col_types = FALSE)
  
  # Simulate annotation join logic
  anno_candidates <- vip_df %>%
    left_join(limma_df,  by = "feature_id") %>%
    left_join(featdef_df, by = "feature_id")
  
  # 1) All VIP features should be present → 3 rows
  expect_equal(nrow(anno_candidates), 3)
  
  # 2) Expected number of columns after two joins
  expect_equal(ncol(anno_candidates), 7)
  
  # 3) F003 has no limma stats → logFC should be NA
  f003_data <- anno_candidates %>% filter(feature_id == "F003")
  expect_true(is.na(f003_data$logFC))
  
  # 4) F001: check correct mzmed join
  f001_data <- anno_candidates %>% filter(feature_id == "F001")
  expect_equal(f001_data$mzmed, 150.05)
})