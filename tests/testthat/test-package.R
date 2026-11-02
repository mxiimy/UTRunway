test_that("package namespace loads", {
  expect_true(requireNamespace("UTRunway", quietly = TRUE))
})

test_that("package metadata is as declared", {
  expect_identical(
    object = utils::packageDescription(pkg = "UTRunway")$Package,
    expected = "UTRunway"
  )
})
