# Independent GC count, to check the function's own figure.
countGcFraction <- function(sequence) {
  letters_in <- strsplit(x = sequence, split = "")[[1]]
  mean(letters_in %in% c("G", "C"))
}

test_that("the result has the documented columns and types", {
  result <- estimateStructure(
    sequence = paste(rep(x = "ACGT", times = 20), collapse = "")
  )

  expect_s3_class(object = result, class = "data.frame")
  expect_identical(
    object = names(x = result),
    expected = c(
      "window_start", "window_end", "gc_content", "max_stem_nt",
      "hairpin_flag"
    )
  )
  expect_type(object = result$window_start, type = "integer")
  expect_type(object = result$window_end, type = "integer")
  expect_type(object = result$gc_content, type = "double")
  expect_type(object = result$max_stem_nt, type = "integer")
  expect_type(object = result$hairpin_flag, type = "logical")
})

test_that("hairpin_flag is max_stem_nt measured against min_stem", {
  data(simulatedConstructs, package = "UTRunway")
  utr5 <- simulatedConstructs$utr5[[5]]

  for (stem in c(4L, 6L, 8L)) {
    result <- estimateStructure(sequence = utr5, min_stem = stem)
    expect_identical(
      object = result$hairpin_flag,
      expected = result$max_stem_nt >= stem,
      info = paste("min_stem", stem)
    )
  }

  # max_stem_nt is a measurement, so it does not move when the
  # threshold does.
  expect_identical(
    object = estimateStructure(sequence = utr5, min_stem = 4)$max_stem_nt,
    expected = estimateStructure(sequence = utr5, min_stem = 9)$max_stem_nt
  )
  expect_true(object = all(result$max_stem_nt >= 0L))
})

test_that("max_stem_nt finds the longest complementary run", {
  # Two eight nucleotide complementary runs around a three
  # nucleotide loop, and nothing longer.
  single <- estimateStructure(
    sequence = "GGGGCCCCAAAGGGGCCCC", window = 19
  )
  expect_identical(object = single$max_stem_nt, expected = 8L)
  expect_true(object = single$hairpin_flag)

  # No complementarity at all.
  flat <- estimateStructure(
    sequence = paste(rep(x = "A", times = 30), collapse = ""),
    window = 30
  )
  expect_identical(object = flat$max_stem_nt, expected = 0L)
  expect_false(object = flat$hairpin_flag)

  # A stem cannot exceed what two copies plus a loop allow.
  result <- estimateStructure(
    sequence = paste(rep(x = "ACGT", times = 10), collapse = ""),
    window = 40,
    min_loop = 3
  )
  expect_true(object = all(result$max_stem_nt <= (40L - 3L) %/% 2L))
})

test_that("windows follow the requested width and step", {
  sequence <- paste(rep(x = "ACGT", times = 25), collapse = "")

  result <- estimateStructure(sequence = sequence, window = 40, step = 10)
  expect_identical(
    object = result$window_start,
    expected = seq(from = 1L, to = nchar(x = sequence) - 40L + 1L, by = 10L)
  )
  expect_identical(
    object = result$window_end, expected = result$window_start + 39L
  )
  expect_true(object = all(result$window_end <= nchar(x = sequence)))

  # Non-overlapping windows when step equals window.
  tiled <- estimateStructure(sequence = sequence, window = 20, step = 20)
  expect_identical(object = tiled$window_start, expected = c(1L, 21L, 41L,
                                                             61L, 81L))
})

test_that("trailing bases that cannot fill a window are not covered", {
  # 25 nt with 10 nt windows stepping by 10 covers 1-10 and 11-20; the
  # last five bases cannot fill a window.
  result <- estimateStructure(
    sequence = paste(rep(x = "ACGTA", times = 5), collapse = ""),
    window = 10,
    step = 10,
    min_stem = 3,
    min_loop = 3
  )

  expect_identical(object = nrow(x = result), expected = 2L)
  expect_identical(object = max(result$window_end), expected = 20L)
})

test_that("GC content is counted correctly", {
  single <- estimateStructure(
    sequence = "GGCCAAAATT", window = 10, min_stem = 3, min_loop = 3
  )
  expect_equal(object = single$gc_content, expected = 0.4)

  expect_equal(
    object = estimateStructure(
      sequence = paste(rep(x = "GC", times = 10), collapse = ""),
      window = 20,
      min_stem = 3,
      min_loop = 3
    )$gc_content,
    expected = 1
  )
  expect_equal(
    object = estimateStructure(
      sequence = paste(rep(x = "AT", times = 10), collapse = ""),
      window = 20,
      min_stem = 3,
      min_loop = 3
    )$gc_content,
    expected = 0
  )

  # Every window of a real sequence, against an independent count.
  data(simulatedConstructs, package = "UTRunway")
  utr5 <- simulatedConstructs$utr5[[5]]
  result <- estimateStructure(sequence = utr5)
  expected_gc <- vapply(
    X = seq_len(length.out = nrow(x = result)),
    FUN = function(i) {
      countGcFraction(sequence = substr(
        x = utr5,
        start = result$window_start[[i]],
        stop = result$window_end[[i]]
      ))
    },
    FUN.VALUE = numeric(length = 1L)
  )
  expect_equal(object = result$gc_content, expected = expected_gc)
  expect_true(
    object = all(result$gc_content >= 0 & result$gc_content <= 1)
  )
})

test_that("a designed hairpin is flagged and a flat sequence is not", {
  # Two 5 nt complementary runs around a 3 nt loop.
  hairpin <- estimateStructure(
    sequence = "GGGGCCCCAAAGGGGCCCCTTTTTTTTTTTTTTTTTTTT",
    window = 39
  )
  expect_true(object = hairpin$hairpin_flag)

  # No complementary runs at all.
  flat <- estimateStructure(
    sequence = paste(rep(x = "A", times = 50), collapse = ""),
    window = 25,
    step = 25
  )
  expect_false(object = any(flat$hairpin_flag))
})

test_that("the flagged stem and partner really are complementary", {
  data(simulatedConstructs, package = "UTRunway")
  utr5 <- simulatedConstructs$utr5[[5]]
  # A five nucleotide stem is common enough that flagged windows
  # exist to check; the default of eight is deliberately stricter.
  result <- estimateStructure(sequence = utr5, min_stem = 5)
  flagged <- which(result$hairpin_flag)
  skip_if(
    condition = length(x = flagged) == 0L,
    message = "no window was flagged in this sequence"
  )

  # Re-find a pair independently of the function, and confirm the
  # partner is the exact reverse complement at a legal distance.
  for (i in flagged) {
    window_sequence <- substr(
      x = utr5,
      start = result$window_start[[i]],
      stop = result$window_end[[i]]
    )
    found_pair <- FALSE
    for (stem_start in seq_len(length.out = nchar(x = window_sequence) - 4L)) {
      stem <- substr(
        x = window_sequence, start = stem_start, stop = stem_start + 4L
      )
      partner <- as.character(x = Biostrings::reverseComplement(
        x = Biostrings::DNAString(x = stem)
      ))
      search_from <- stem_start + 5L + 3L
      if (search_from > nchar(x = window_sequence) - 4L) {
        next
      }
      downstream <- substr(
        x = window_sequence,
        start = search_from,
        stop = nchar(x = window_sequence)
      )
      if (grepl(pattern = partner, x = downstream, fixed = TRUE)) {
        found_pair <- TRUE
        break
      }
    }
    expect_true(object = found_pair, info = paste("window", i))
  }
})

test_that("min_stem and min_loop change what counts as a hairpin", {
  # A 5 nt stem pair separated by a 2 nt loop.
  sequence <- "GGGGGAACCCCCAAAAAA"

  expect_false(
    object = estimateStructure(
      sequence = sequence, window = 18, min_stem = 5, min_loop = 3
    )$hairpin_flag
  )
  expect_true(
    object = estimateStructure(
      sequence = sequence, window = 18, min_stem = 5, min_loop = 2
    )$hairpin_flag
  )

  # The same shift shows up in the measurement, not just the flag.
  # min_stem is set low only so that the window is wide enough to
  # hold two stems and a loop, which avoids an unrelated warning.
  expect_identical(
    object = estimateStructure(
      sequence = sequence, window = 18, min_stem = 4, min_loop = 3
    )$max_stem_nt,
    expected = 4L
  )
  expect_identical(
    object = estimateStructure(
      sequence = sequence, window = 18, min_stem = 4, min_loop = 2
    )$max_stem_nt,
    expected = 5L
  )

  # Demanding a longer stem can only ever remove flags, never add
  # them, since every long stem contains a shorter one.
  data(simulatedConstructs, package = "UTRunway")
  utr5 <- simulatedConstructs$utr5[[5]]
  lenient <- estimateStructure(sequence = utr5, min_stem = 4)
  strict <- estimateStructure(sequence = utr5, min_stem = 8)

  expect_identical(
    object = nrow(x = strict), expected = nrow(x = lenient)
  )
  expect_true(object = all(strict$hairpin_flag <= lenient$hairpin_flag))
  expect_true(object = any(lenient$hairpin_flag))
})

test_that("a sequence shorter than the window warns and gives one row", {
  # Two warnings: the sequence is shorter than the default window, and
  # the window it falls back to cannot hold a stem-loop.
  warnings_raised <- capture_warnings(
    code = estimateStructure(sequence = "ACGTACGTAC")
  )
  expect_length(object = warnings_raised, n = 2L)
  expect_match(
    object = warnings_raised,
    regexp = "only 10 nucleotides long",
    all = FALSE
  )
  expect_match(
    object = warnings_raised, regexp = "narrower than", all = FALSE
  )

  result <- suppressWarnings(
    estimateStructure(sequence = "ACGTACGTAC")
  )
  expect_identical(object = nrow(x = result), expected = 1L)
  expect_identical(object = result$window_start, expected = 1L)
  expect_identical(object = result$window_end, expected = 10L)
})

test_that("a window too narrow to hold a stem-loop warns", {
  expect_warning(
    object = estimateStructure(
      sequence = paste(rep(x = "ACGT", times = 20), collapse = ""),
      window = 8,
      step = 8
    ),
    regexp = "narrower than the 19 nucleotides"
  )
  result <- suppressWarnings(
    estimateStructure(
      sequence = paste(rep(x = "ACGT", times = 20), collapse = ""),
      window = 8,
      step = 8
    )
  )
  expect_false(object = any(result$hairpin_flag))
})

test_that("every construct in the dataset scans without complaint", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    for (region in c("utr5", "cds", "utr3")) {
      expect_no_warning(
        object = estimateStructure(
          sequence = simulatedConstructs[[region]][[i]]
        )
      )
      result <- estimateStructure(
        sequence = simulatedConstructs[[region]][[i]]
      )
      expect_gt(object = nrow(x = result), expected = 0)
      expect_false(object = anyNA(x = result))
    }
  }
})

test_that("lower case and RNA input give the same result", {
  expect_identical(
    object = estimateStructure(
      sequence = "ggggccccaaaggggccccuuuuuuuuuuuuuuuuuuuu", window = 39
    ),
    expected = estimateStructure(
      sequence = "GGGGCCCCAAAGGGGCCCCTTTTTTTTTTTTTTTTTTTT", window = 39
    )
  )
})

test_that("invalid input is rejected with an informative message", {
  sequence <- paste(rep(x = "ACGT", times = 20), collapse = "")

  expect_error(
    object = estimateStructure(sequence = NULL), regexp = "not NULL"
  )
  expect_error(
    object = estimateStructure(sequence = NA), regexp = "must not be NA"
  )
  expect_error(
    object = estimateStructure(sequence = "ACGTN"),
    regexp = "unambiguous nucleotides"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, window = 0),
    regexp = "`window` must be at least 1"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, step = -3),
    regexp = "`step` must be at least 1"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, window = "wide"),
    regexp = "`window` must be a single number"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, step = 2.5),
    regexp = "`step` must be a whole number"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, min_stem = 1),
    regexp = "`min_stem` must be at least 2"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, min_loop = 0),
    regexp = "`min_loop` must be at least 1"
  )
  expect_error(
    object = estimateStructure(sequence = sequence, window = c(10, 20)),
    regexp = "`window` must be a single number"
  )
})
