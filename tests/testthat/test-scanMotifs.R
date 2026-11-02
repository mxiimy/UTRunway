expected_motif_columns <- c(
  "motif", "region", "start", "end", "width", "match", "pattern"
)

test_that("the result has the documented columns and types", {
  result <- scanMotifs(
    sequence = "CTTTTCTTTCTAAAAAAAAAA", region = "utr5"
  )

  expect_s3_class(object = result, class = "data.frame")
  expect_identical(
    object = names(x = result), expected = expected_motif_columns
  )
  expect_type(object = result$motif, type = "character")
  expect_type(object = result$region, type = "character")
  expect_type(object = result$start, type = "integer")
  expect_type(object = result$end, type = "integer")
  expect_type(object = result$width, type = "integer")
  expect_type(object = result$match, type = "character")
  expect_type(object = result$pattern, type = "character")
})

test_that("no match gives a zero-row data frame", {
  result <- scanMotifs(
    sequence = paste(rep(x = "A", times = 40), collapse = ""),
    region = "utr5"
  )

  expect_false(object = is.null(x = result))
  expect_s3_class(object = result, class = "data.frame")
  expect_identical(object = nrow(x = result), expected = 0L)
  expect_identical(
    object = names(x = result), expected = expected_motif_columns
  )

  populated <- scanMotifs(
    sequence = "CTTTTCTTTCTAAAAAAAAAA", region = "utr5"
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

test_that("built-in motifs match the planted truth", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- scanMotifs(sequence = row_i$utr5, region = "utr5")

    expect_identical(
      object = any(result$motif == "TOP"),
      expected = row_i$has_top_motif,
      info = row_i$construct_id
    )
    expect_identical(
      object = any(result$motif == "G_quadruplex"),
      expected = row_i$has_g_quadruplex,
      info = row_i$construct_id
    )
  }
})

test_that("reported matches really occur where they are reported", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    for (region in c("utr5", "cds", "utr3")) {
      sequence <- simulatedConstructs[[region]][[i]]
      result <- scanMotifs(sequence = sequence, region = region)
      if (nrow(x = result) == 0L) {
        next
      }

      expect_identical(
        object = substring(
          text = sequence, first = result$start, last = result$end
        ),
        expected = result$match,
        info = paste(simulatedConstructs$construct_id[[i]], region)
      )
      expect_identical(
        object = result$width,
        expected = result$end - result$start + 1L
      )
      expect_true(object = all(result$region == region))
      expect_identical(
        object = result$start, expected = sort(x = result$start)
      )
    }
  }
})

test_that("the TOP motif is anchored and only scanned in a 5' UTR", {
  sequence <- "CTTTTCTTTCTAAAAAAAAAA"

  as_utr5 <- scanMotifs(sequence = sequence, region = "utr5")
  expect_identical(object = as_utr5$motif, expected = "TOP")
  expect_identical(object = as_utr5$start, expected = 1L)
  expect_identical(object = as_utr5$match, expected = "CTTTTCTTTCT")

  # The same sequence is not scanned for TOP as another region.
  expect_identical(
    object = nrow(x = scanMotifs(sequence = sequence, region = "utr3")),
    expected = 0L
  )
  expect_identical(
    object = nrow(x = scanMotifs(sequence = sequence, region = "cds")),
    expected = 0L
  )

  # A pyrimidine stretch that does not start at position 1 is not a
  # TOP motif, because the motif is defined by its position.
  expect_identical(
    object = nrow(x = scanMotifs(
      sequence = "AACTTTTCTTTCTAAAAAAAA", region = "utr5"
    )),
    expected = 0L
  )

  # A C followed by too few pyrimidines does not qualify either.
  expect_identical(
    object = nrow(x = scanMotifs(
      sequence = "CTTAGGGAAAAAAAAAAAAA", region = "utr5"
    )),
    expected = 0L
  )
})

test_that("the G-quadruplex pattern is scanned in every region", {
  sequence <- "AAGGGAGGGAGGGAGGGAAA"

  for (region in c("utr5", "cds", "utr3")) {
    result <- scanMotifs(sequence = sequence, region = region)
    quad <- result[result$motif == "G_quadruplex", ]

    expect_identical(object = nrow(x = quad), expected = 1L)
    expect_identical(object = quad$start, expected = 3L)
    expect_identical(object = quad$end, expected = 17L)
    expect_identical(
      object = quad$match, expected = "GGGAGGGAGGGAGGG"
    )
  }

  # Three G-triplets are not enough.
  expect_identical(
    object = nrow(x = scanMotifs(
      sequence = "AAGGGAGGGAGGGAAACCCC", region = "utr5"
    )),
    expected = 0L
  )

  # A loop longer than seven nucleotides breaks the arrangement.
  expect_identical(
    object = nrow(x = scanMotifs(
      sequence = "GGGAAAAAAAAGGGAGGGAGGGAA", region = "utr5"
    )),
    expected = 0L
  )
})

test_that("supplied motifs replace the built-in set", {
  sequence <- "CTTTTCTTTCTATTTAAAAAAAAA"

  built_in <- scanMotifs(sequence = sequence, region = "utr5")
  expect_true(object = any(built_in$motif == "TOP"))

  supplied <- scanMotifs(
    sequence = sequence, region = "utr5", motifs = c(are = "ATTTA")
  )
  expect_true(object = all(supplied$motif == "are"))
  expect_false(object = any(supplied$motif == "TOP"))
})

test_that("supplied motifs honour IUPAC ambiguity codes", {
  # Y is C or T, so both forms of the site match.
  result <- scanMotifs(
    sequence = "GGACTAACGGACTAATGG",
    region = "utr3",
    motifs = c(y = "ACTAAY")
  )
  expect_identical(object = nrow(x = result), expected = 2L)
  expect_identical(
    object = result$match, expected = c("ACTAAC", "ACTAAT")
  )
  expect_true(object = all(result$pattern == "ACTAAY"))

  # N matches any base, so it finds the same two sites.
  any_base <- scanMotifs(
    sequence = "GGACTAACGGACTAATGG",
    region = "utr3",
    motifs = c(n = "ACTAAN")
  )
  expect_identical(object = any_base$start, expected = c(3L, 11L))

  # An unambiguous motif matches only its exact form.
  exact <- scanMotifs(
    sequence = "GGACTAACGGACTAATGG",
    region = "utr3",
    motifs = c(exact = "ACTAAC")
  )
  expect_identical(object = nrow(x = exact), expected = 1L)
  expect_identical(object = exact$start, expected = 3L)
})

test_that("overlapping matches are each reported", {
  result <- scanMotifs(
    sequence = "ATTTATTTATTTA",
    region = "utr3",
    motifs = c(are = "ATTTA")
  )

  expect_identical(object = nrow(x = result), expected = 3L)
  expect_identical(object = result$start, expected = c(1L, 5L, 9L))
})

test_that("several motifs are returned together in positional order", {
  result <- scanMotifs(
    sequence = "GGACTAACGGACTAATGG",
    region = "utr3",
    motifs = c(y = "ACTAAY", n = "ACTAAN")
  )

  expect_identical(object = nrow(x = result), expected = 4L)
  expect_identical(object = result$start, expected = c(3L, 3L, 11L, 11L))
  # Ties on position are broken by motif name, so the order is stable.
  expect_identical(
    object = result$motif, expected = c("n", "y", "n", "y")
  )
})

test_that("region is matched exactly, with no partial matching", {
  sequence <- "CTTTTCTTTCTAAAAAAAAAA"

  expect_error(
    object = scanMotifs(sequence = sequence, region = "utr"),
    regexp = "must be one of"
  )
  expect_error(
    object = scanMotifs(sequence = sequence, region = "UTR5"),
    regexp = "must be one of"
  )
  expect_error(
    object = scanMotifs(sequence = sequence, region = "five_prime"),
    regexp = "must be one of"
  )
  expect_error(
    object = scanMotifs(sequence = sequence, region = NULL),
    regexp = "not NULL"
  )
  expect_error(
    object = scanMotifs(sequence = sequence, region = c("utr5", "cds")),
    regexp = "must be a single character string"
  )
  expect_error(
    object = scanMotifs(sequence = sequence, region = 1),
    regexp = "must be a single character string"
  )
})

test_that("lower case and RNA input give the same result", {
  expect_identical(
    object = scanMotifs(sequence = "cuuuucuuucuaaaaaaaaaa", region = "utr5"),
    expected = scanMotifs(
      sequence = "CTTTTCTTTCTAAAAAAAAAA", region = "utr5"
    )
  )
  # A motif written as RNA finds the same sites as its DNA form.
  expect_identical(
    object = scanMotifs(
      sequence = "ATTTATTTATTTA", region = "utr3",
      motifs = c(are = "auuua")
    )$start,
    expected = c(1L, 5L, 9L)
  )
})

test_that("invalid motifs are rejected with an informative message", {
  sequence <- "CTTTTCTTTCTAAAAAAAAAA"

  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c("ATG", "GGG")
    ),
    regexp = "must be named"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5",
      motifs = c(a = "ATG", a = "GGG")
    ),
    regexp = "must have unique names"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(a = "ATG", "GGG")
    ),
    regexp = "must be named"
  )
  # A regular expression is not an IUPAC motif, and saying so is the
  # point: it would otherwise silently fail to match.
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(bad = "AT[GC]")
    ),
    regexp = "not regular expressions"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(bad = "G{3,}")
    ),
    regexp = "IUPAC nucleotide codes"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = character(0)
    ),
    regexp = "at least one motif"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(a = NA_character_)
    ),
    regexp = "must not contain NA"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(a = "")
    ),
    regexp = "must not contain empty motifs"
  )
  expect_error(
    object = scanMotifs(
      sequence = sequence, region = "utr5", motifs = c(a = 42)
    ),
    regexp = "must be a named character vector"
  )
})

test_that("a motif longer than the sequence matches nothing", {
  expect_identical(
    object = nrow(x = scanMotifs(
      sequence = "ACGT", region = "utr3", motifs = c(long = "ACGTACGTACGT")
    )),
    expected = 0L
  )
})
