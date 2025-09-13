  # We locate files relative to tests to parse them here
  this_test_dir <- dirname(test_path())
  r_dir <- file.path(this_test_dir, "..", "R")
  r_files <- list.files(r_dir, pattern = "\\.[rR]$", full.names = TRUE)

test_that("Package files are parsable", {
  expect_true(requireNamespace("plsRcox", quietly = TRUE))
  #  rdir <- system.file("..", package = "base")
  for (f in r_files) {
      expect_no_error(parse(file = f))
    }
  }
)
