test_that("the return value has the documented shape", {
  result <- scoreKozak(sequence = "GCCACCATGGCC", position = 7)

  expect_type(object = result, type = "list")
  expect_identical(
    object = names(x = result),
    expected = c("score", "score_scaled", "strength", "context")
  )
  expect_type(object = result$score, type = "double")
  expect_type(object = result$score_scaled, type = "double")
  expect_type(object = result$strength, type = "character")
  expect_type(object = result$context, type = "character")
})

test_that("the full consensus context scores the maximum", {
  # GCCACCATGG is the consensus at every scored position, so it must
  # sit at the top of the scale.
  result <- scoreKozak(sequence = "GCCACCATGGCC", position = 7)

  expect_identical(object = result$strength, expected = "strong")
  expect_equal(object = result$score_scaled, expected = 1)
  expect_gt(object = result$score, expected = 0)
  expect_identical(object = result$context, expected = "gccacc ATG g")
})

test_that("strength follows the -3 purine and +4 G rule", {
  # Only -3 and +4 differ between these four contexts.
  expect_identical(
    object = scoreKozak(
      sequence = "CCCACCATGGCC", position = 7
    )$strength,
    expected = "strong"
  )
  expect_identical(
    object = scoreKozak(
      sequence = "CCCTCCATGGCC", position = 7
    )$strength,
    expected = "adequate"
  )
  expect_identical(
    object = scoreKozak(
      sequence = "CCCACCATGACC", position = 7
    )$strength,
    expected = "adequate"
  )
  expect_identical(
    object = scoreKozak(
      sequence = "CCCTCCATGACC", position = 7
    )$strength,
    expected = "weak"
  )

  # A purine at -3 is either base, so G must behave as A does.
  expect_identical(
    object = scoreKozak(
      sequence = "CCCGCCATGGCC", position = 7
    )$strength,
    expected = "strong"
  )
})

test_that("a better context scores higher than a worse one", {
  strong <- scoreKozak(sequence = "GCCACCATGGCC", position = 7)
  adequate <- scoreKozak(sequence = "CCCTCCATGGCC", position = 7)
  weak <- scoreKozak(sequence = "TTTTTTATGCCC", position = 7)

  expect_gt(object = strong$score, expected = adequate$score)
  expect_gt(object = adequate$score, expected = weak$score)
  expect_true(
    object = all(
      c(strong$score_scaled, adequate$score_scaled, weak$score_scaled) >= 0 &
        c(strong$score_scaled, adequate$score_scaled, weak$score_scaled) <= 1
    )
  )
})

test_that("planted Kozak contexts score as planted", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- scoreKozak(
      sequence = paste0(row_i$utr5, row_i$cds),
      position = nchar(x = row_i$utr5) + 1L
    )

    if (row_i$has_strong_kozak) {
      expect_identical(
        object = result$strength,
        expected = "strong",
        info = row_i$construct_id
      )
      expect_equal(
        object = result$score_scaled,
        expected = 1,
        info = row_i$construct_id
      )
    }
    if (row_i$has_weak_kozak) {
      expect_identical(
        object = result$strength,
        expected = "weak",
        info = row_i$construct_id
      )
      expect_lt(
        object = result$score,
        expected = 0,
        label = row_i$construct_id
      )
    }
  }
})

test_that("lower case and RNA input score the same as upper case DNA", {
  expect_identical(
    object = scoreKozak(sequence = "gccaccauggcc", position = 7),
    expected = scoreKozak(sequence = "GCCACCATGGCC", position = 7)
  )
})

test_that("a context window that does not fit gives NA and a warning", {
  # Too few bases before the AUG.
  expect_warning(
    object = scoreKozak(sequence = "AATGGCCCCC", position = 2),
    regexp = "does not fit"
  )
  near_start <- suppressWarnings(
    scoreKozak(sequence = "AATGGCCCCC", position = 2)
  )
  expect_identical(object = near_start$score, expected = NA_real_)
  expect_identical(object = near_start$score_scaled, expected = NA_real_)
  expect_identical(object = near_start$strength, expected = NA_character_)
  expect_identical(object = near_start$context, expected = NA_character_)

  # Enough bases before the AUG, but no base at +4.
  expect_warning(
    object = scoreKozak(sequence = "GCCACCATG", position = 7),
    regexp = "does not fit"
  )
  near_end <- suppressWarnings(
    scoreKozak(sequence = "GCCACCATG", position = 7)
  )
  expect_identical(object = near_end$strength, expected = NA_character_)

  # One more base is enough to score it.
  expect_no_warning(
    object = scoreKozak(sequence = "GCCACCATGG", position = 7)
  )
})

test_that("invalid input is rejected with an informative message", {
  expect_error(
    object = scoreKozak(sequence = "GCCACCATGGCC", position = 8),
    regexp = "must point at the A of an ATG start codon"
  )
  expect_error(
    object = scoreKozak(sequence = "ACGTACGTACGTATG", position = 14),
    regexp = "leaves no room for a start codon"
  )
  expect_error(
    object = scoreKozak(sequence = "GCCACCATGGCC", position = 7.5),
    regexp = "must be a whole number"
  )
  expect_error(
    object = scoreKozak(sequence = "GCCACCATGGCC", position = 0),
    regexp = "must be at least 1"
  )
  expect_error(
    object = scoreKozak(sequence = "GCCACCATGGCC", position = "7"),
    regexp = "must be a single number"
  )
  expect_error(
    object = scoreKozak(sequence = 42, position = 7),
    regexp = "must be a character string"
  )
  expect_error(
    object = scoreKozak(sequence = NULL, position = 7),
    regexp = "not NULL"
  )
  expect_error(
    object = scoreKozak(sequence = NA, position = 7),
    regexp = "must not be NA"
  )
  expect_error(
    object = scoreKozak(sequence = c("GCCACCATGGCC", "GCCACCATGGCC"),
                        position = 7),
    regexp = "must be a single sequence"
  )
  expect_error(
    object = scoreKozak(sequence = "GCCACCATGGNN", position = 7),
    regexp = "unambiguous nucleotides"
  )
  expect_error(
    object = scoreKozak(sequence = "", position = 7),
    regexp = "must not be an empty string"
  )
})
