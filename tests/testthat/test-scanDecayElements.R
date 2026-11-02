expected_decay_columns <- c(
  "element_type", "element_name", "start", "end", "width", "sequence",
  "distance_to_end", "is_canonical", "near_three_prime_end",
  "in_cluster", "distance_to_nearest_are"
)

# Overlapping motif count, computed independently of the package.
countOverlapping <- function(sequence, motif) {
  n_starts <- nchar(x = sequence) - nchar(x = motif) + 1L
  if (n_starts < 1L) {
    return(0L)
  }
  starts <- seq_len(length.out = n_starts)
  sum(
    substring(
      text = sequence, first = starts, last = starts + nchar(x = motif) - 1L
    ) == motif
  )
}

test_that("the result has the documented columns and types", {
  result <- scanDecayElements(
    utr3 = paste0(
      "CCCCCCCCCCATTTA", paste(rep(x = "C", times = 80), collapse = ""),
      "AATAAACCCCCCCCCCCCCCCCCCCC"
    )
  )

  expect_s3_class(object = result, class = "data.frame")
  expect_identical(
    object = names(x = result), expected = expected_decay_columns
  )
  expect_type(object = result$element_type, type = "character")
  expect_type(object = result$element_name, type = "character")
  expect_type(object = result$start, type = "integer")
  expect_type(object = result$end, type = "integer")
  expect_type(object = result$width, type = "integer")
  expect_type(object = result$sequence, type = "character")
  expect_type(object = result$distance_to_end, type = "integer")
  expect_type(object = result$is_canonical, type = "logical")
  expect_type(object = result$near_three_prime_end, type = "logical")
  expect_type(object = result$in_cluster, type = "logical")
  expect_type(object = result$distance_to_nearest_are, type = "integer")
})

test_that("a 3' UTR with no elements gives a zero-row data frame", {
  result <- scanDecayElements(
    utr3 = paste(rep(x = "CCGG", times = 20), collapse = ""),
    seeds = character(0)
  )

  expect_false(object = is.null(x = result))
  expect_s3_class(object = result, class = "data.frame")
  expect_identical(object = nrow(x = result), expected = 0L)
  expect_identical(
    object = names(x = result), expected = expected_decay_columns
  )

  populated <- scanDecayElements(
    utr3 = paste0(
      "CCCCCCCCCCATTTA", paste(rep(x = "C", times = 80), collapse = ""),
      "AATAAACCCCCCCCCCCCCCCCCCCC"
    )
  )
  expect_identical(
    object = vapply(
      X = result, FUN = typeof, FUN.VALUE = character(length = 1L)
    ),
    expected = vapply(
      X = populated, FUN = typeof, FUN.VALUE = character(length = 1L)
    )
  )
})

test_that("AU-rich pentamer counts match the planted truth", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- scanDecayElements(utr3 = row_i$utr3)
    pentamers <- result[result$element_type == "ARE_pentamer", ]

    expect_identical(
      object = nrow(x = pentamers),
      expected = row_i$n_are_pentamer,
      info = row_i$construct_id
    )
    # Overlapping pentamers are counted separately, as the dataset's
    # own truth column does.
    expect_identical(
      object = nrow(x = pentamers),
      expected = countOverlapping(sequence = row_i$utr3, motif = "ATTTA"),
      info = row_i$construct_id
    )
  }
})

test_that("canonical polyA signals match the planted truth", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- scanDecayElements(utr3 = row_i$utr3)

    has_canonical <- any(
      result$is_canonical & result$near_three_prime_end,
      na.rm = TRUE
    )
    expect_identical(
      object = has_canonical,
      expected = row_i$has_canonical_polya,
      info = row_i$construct_id
    )
  }
})

test_that("every element is reported where it actually occurs", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- scanDecayElements(utr3 = row_i$utr3)
    if (nrow(x = result) == 0L) {
      next
    }

    expect_identical(
      object = substring(
        text = row_i$utr3, first = result$start, last = result$end
      ),
      expected = result$sequence,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$width,
      expected = result$end - result$start + 1L,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$distance_to_end,
      expected = nchar(x = row_i$utr3) - result$end,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$start,
      expected = sort(x = result$start),
      info = row_i$construct_id
    )
  }
})

test_that("is_canonical applies to polyA signals only", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    result <- scanDecayElements(utr3 = simulatedConstructs$utr3[[i]])
    if (nrow(x = result) == 0L) {
      next
    }
    is_polya <- result$element_type == "polyA_signal"

    expect_false(object = anyNA(x = result$is_canonical[is_polya]))
    expect_true(object = all(is.na(x = result$is_canonical[!is_polya])))
    # The canonical flag follows the hexamer that matched.
    expect_identical(
      object = result$is_canonical[is_polya],
      expected = result$element_name[is_polya] == "AATAAA"
    )
  }
})

test_that("the ATTAAA variant is found and marked non-canonical", {
  utr3 <- paste0(
    paste(rep(x = "C", times = 60), collapse = ""),
    "ATTAAA",
    paste(rep(x = "C", times = 20), collapse = "")
  )
  result <- scanDecayElements(utr3 = utr3, seeds = character(0))
  signals <- result[result$element_type == "polyA_signal", ]

  expect_identical(object = nrow(x = signals), expected = 1L)
  expect_identical(object = signals$element_name, expected = "ATTAAA")
  expect_false(object = signals$is_canonical)
  expect_true(object = signals$near_three_prime_end)
})

test_that("clustering and the distance to the nearest pentamer agree", {
  # A single pentamer: nothing to cluster with.
  lone <- scanDecayElements(
    utr3 = paste0(
      "CCCCCCCCCCATTTA", paste(rep(x = "C", times = 80), collapse = ""),
      "AATAAACCCCCCCCCCCCCCCCCCCC"
    ),
    seeds = character(0)
  )
  lone_are <- lone[lone$element_type == "ARE_pentamer", ]
  expect_identical(object = nrow(x = lone_are), expected = 1L)
  expect_identical(object = lone_are$start, expected = 11L)
  expect_false(object = lone_are$in_cluster)
  expect_identical(
    object = lone_are$distance_to_nearest_are, expected = NA_integer_
  )

  # One run of overlapping pentamers, four nucleotides apart.
  run <- scanDecayElements(
    utr3 = paste0(
      paste(rep(x = "C", times = 30), collapse = ""),
      "TATTTATTTATTTATT",
      paste(rep(x = "C", times = 40), collapse = ""),
      "AATAAACCCCCCCCCCCCCCCCCCCC"
    ),
    seeds = character(0)
  )
  run_are <- run[run$element_type == "ARE_pentamer", ]
  expect_identical(object = nrow(x = run_are), expected = 3L)
  expect_identical(object = run_are$start, expected = c(32L, 36L, 40L))
  expect_true(object = all(run_are$in_cluster))
  expect_identical(
    object = run_are$distance_to_nearest_are, expected = c(4L, 4L, 4L)
  )

  # in_cluster is defined for pentamers only.
  expect_true(
    object = all(is.na(
      x = run$in_cluster[run$element_type != "ARE_pentamer"]
    ))
  )
  expect_true(
    object = all(is.na(
      x = run$distance_to_nearest_are[
        run$element_type != "ARE_pentamer"
      ]
    ))
  )
})

test_that("in_cluster follows are_cluster_fraction", {
  data(simulatedConstructs, package = "UTRunway")
  utr3 <- simulatedConstructs$utr3[[12]]

  # A fraction of zero makes the window zero nucleotides wide, so no
  # pentamer can cluster with any other.
  none <- scanDecayElements(utr3 = utr3, are_cluster_fraction = 0)
  expect_false(
    object = any(none$in_cluster, na.rm = TRUE)
  )

  # A fraction of one spans the region, so every pentamer clusters.
  all_clustered <- scanDecayElements(
    utr3 = utr3, are_cluster_fraction = 1
  )
  pentamers <- all_clustered[
    all_clustered$element_type == "ARE_pentamer",
  ]
  expect_true(object = all(pentamers$in_cluster))

  # The measured distances do not depend on the threshold.
  expect_identical(
    object = none$distance_to_nearest_are,
    expected = all_clustered$distance_to_nearest_are
  )
})

test_that("near_three_prime_end follows polya_end_window", {
  data(simulatedConstructs, package = "UTRunway")
  # The bundled constructs carry their polyA signal 20 nucleotides
  # from the 3' end, which is why the default window is absolute.
  utr3 <- simulatedConstructs$utr3[[1]]

  wide <- scanDecayElements(utr3 = utr3, polya_end_window = 30)
  narrow <- scanDecayElements(utr3 = utr3, polya_end_window = 10)
  signals_wide <- wide[wide$element_type == "polyA_signal", ]
  signals_narrow <- narrow[narrow$element_type == "polyA_signal", ]

  expect_identical(object = signals_wide$distance_to_end, expected = 20L)
  expect_true(object = signals_wide$near_three_prime_end)
  expect_false(object = signals_narrow$near_three_prime_end)

  # The flag is exactly the distance test, for every element.
  expect_identical(
    object = wide$near_three_prime_end,
    expected = wide$distance_to_end <= 30L
  )
})

test_that("the default miRNA seeds are the documented set", {
  seeds <- defaultMirnaSeeds()

  expect_length(object = seeds, n = 5L)
  expect_identical(
    object = seeds,
    expected = c(
      "miR-1/miR-206" = "ACATTCC",
      "let-7" = "CTACCTC",
      "miR-21" = "ATAAGCT",
      "miR-122" = "ACACTCC",
      "miR-155" = "AGCATTA"
    )
  )
  # Each is a 7mer-m8 site, so seven nucleotides of unambiguous DNA.
  expect_true(object = all(nchar(x = seeds) == 7L))
  expect_false(object = any(grepl(pattern = "[^ACGT]", x = seeds)))
})

test_that("a planted seed site is found and named", {
  utr3 <- paste0(
    paste(rep(x = "C", times = 40), collapse = ""),
    "AGCATTA",
    paste(rep(x = "C", times = 40), collapse = "")
  )
  result <- scanDecayElements(utr3 = utr3)
  seed_hits <- result[result$element_type == "miRNA_seed", ]

  expect_identical(object = nrow(x = seed_hits), expected = 1L)
  expect_identical(object = seed_hits$element_name, expected = "miR-155")
  expect_identical(object = seed_hits$start, expected = 41L)
  expect_identical(object = seed_hits$sequence, expected = "AGCATTA")
})

test_that("supplied seeds replace the defaults", {
  utr3 <- paste0(
    paste(rep(x = "C", times = 40), collapse = ""),
    "AGCATTA",
    paste(rep(x = "C", times = 40), collapse = "")
  )

  # A seed set that cannot match leaves no miRNA rows, even though
  # the default set would have matched here.
  other <- scanDecayElements(utr3 = utr3, seeds = c(mine = "GGGGGGG"))
  expect_identical(
    object = sum(other$element_type == "miRNA_seed"), expected = 0L
  )

  # IUPAC ambiguity codes are honoured: Y is C or T.
  ambiguous <- scanDecayElements(
    utr3 = "GGACTAACGGACTAATGG", seeds = c(amb = "ACTAAY")
  )
  seed_hits <- ambiguous[ambiguous$element_type == "miRNA_seed", ]
  expect_identical(object = nrow(x = seed_hits), expected = 2L)
  expect_identical(
    object = seed_hits$sequence, expected = c("ACTAAC", "ACTAAT")
  )
  expect_true(object = all(seed_hits$element_name == "amb"))
})

test_that("character(0) skips the miRNA scan only", {
  utr3 <- paste0(
    paste(rep(x = "C", times = 40), collapse = ""),
    "AGCATTAATTTA",
    paste(rep(x = "C", times = 30), collapse = ""),
    "AATAAACCCCCCCCCCCCCCCCCCCC"
  )

  with_seeds <- scanDecayElements(utr3 = utr3)
  without_seeds <- scanDecayElements(utr3 = utr3, seeds = character(0))

  expect_gt(
    object = sum(with_seeds$element_type == "miRNA_seed"), expected = 0
  )
  expect_identical(
    object = sum(without_seeds$element_type == "miRNA_seed"),
    expected = 0L
  )
  # The other element types are unaffected.
  expect_identical(
    object = sum(without_seeds$element_type == "polyA_signal"),
    expected = sum(with_seeds$element_type == "polyA_signal")
  )
  expect_identical(
    object = sum(without_seeds$element_type == "ARE_pentamer"),
    expected = sum(with_seeds$element_type == "ARE_pentamer")
  )
})

test_that("lower case and RNA input give the same result", {
  utr3 <- paste0(
    "CCCCCCCCCCATTTA", paste(rep(x = "C", times = 80), collapse = ""),
    "AATAAACCCCCCCCCCCCCCCCCCCC"
  )
  expect_identical(
    object = scanDecayElements(
      utr3 = tolower(x = gsub(pattern = "T", replacement = "U", x = utr3))
    ),
    expected = scanDecayElements(utr3 = utr3)
  )
})

test_that("invalid input is rejected with an informative message", {
  utr3 <- paste(rep(x = "ACGT", times = 30), collapse = "")

  expect_error(
    object = scanDecayElements(utr3 = NULL), regexp = "not NULL"
  )
  expect_error(
    object = scanDecayElements(utr3 = NA), regexp = "must not be NA"
  )
  expect_error(
    object = scanDecayElements(utr3 = "ACGTN"),
    regexp = "unambiguous nucleotides"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, polya_end_window = -1),
    regexp = "`polya_end_window` must be at least 0"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, polya_end_window = 2.5),
    regexp = "`polya_end_window` must be a whole number"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, are_cluster_fraction = 1.5),
    regexp = "must be a proportion between 0 and 1"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, are_cluster_fraction = "5%"),
    regexp = "must be a single number between 0 and 1"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, seeds = c("ACATTCC")),
    regexp = "`seeds` must be named"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, seeds = c(a = "AC", a = "GG")),
    regexp = "`seeds` must have unique names"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, seeds = c(bad = "AC[GT]")),
    regexp = "IUPAC nucleotide codes"
  )
  expect_error(
    object = scanDecayElements(utr3 = utr3, seeds = c(x = 42)),
    regexp = "`seeds` must be a named character vector"
  )
})
