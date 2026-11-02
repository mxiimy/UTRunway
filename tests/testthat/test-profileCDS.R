expected_cds_fields <- c(
  "cds_length", "n_codons", "gc_content", "gc3_content",
  "starts_with_atg", "ends_with_stop", "is_multiple_of_three",
  "has_premature_stop", "premature_stop_positions",
  "premature_stop_codons", "n_out_of_frame_atg",
  "out_of_frame_atg_positions", "out_of_frame_atg_frames"
)

# GC fraction, counted independently of the package.
gcFraction <- function(sequence) {
  mean(strsplit(x = sequence, split = "")[[1]] %in% c("G", "C"))
}

# Every third base of each complete codon, as one string.
thirdPositions <- function(sequence) {
  positions <- seq(
    from = 3L, by = 3L, length.out = nchar(x = sequence) %/% 3L
  )
  paste(
    substring(text = sequence, first = positions, last = positions),
    collapse = ""
  )
}

test_that("the result has the documented fields and types", {
  result <- profileCDS(cds = "ATGAAATAACCCGGGTAA")

  expect_type(object = result, type = "list")
  expect_identical(
    object = names(x = result), expected = expected_cds_fields
  )
  expect_type(object = result$cds_length, type = "integer")
  expect_type(object = result$n_codons, type = "integer")
  expect_type(object = result$gc_content, type = "double")
  expect_type(object = result$gc3_content, type = "double")
  expect_type(object = result$starts_with_atg, type = "logical")
  expect_type(object = result$ends_with_stop, type = "logical")
  expect_type(object = result$is_multiple_of_three, type = "logical")
  expect_type(object = result$has_premature_stop, type = "logical")
  expect_type(object = result$premature_stop_positions, type = "integer")
  expect_type(object = result$premature_stop_codons, type = "character")
  expect_type(object = result$n_out_of_frame_atg, type = "integer")
  expect_type(
    object = result$out_of_frame_atg_positions, type = "integer"
  )
  expect_type(object = result$out_of_frame_atg_frames, type = "integer")
})

test_that("a worked example gives the values computed by hand", {
  # ATG AAA TAA CCC GGG TAA: a premature stop at codon three, and a
  # terminal stop.
  result <- profileCDS(cds = "ATGAAATAACCCGGGTAA")

  expect_identical(object = result$cds_length, expected = 18L)
  expect_identical(object = result$n_codons, expected = 6L)
  expect_equal(object = result$gc_content, expected = 7 / 18)
  expect_equal(object = result$gc3_content, expected = 0.5)
  expect_true(object = result$starts_with_atg)
  expect_true(object = result$ends_with_stop)
  expect_true(object = result$is_multiple_of_three)
  expect_true(object = result$has_premature_stop)
  expect_identical(
    object = result$premature_stop_positions, expected = 7L
  )
  expect_identical(
    object = result$premature_stop_codons, expected = "TAA"
  )
  expect_identical(object = result$n_out_of_frame_atg, expected = 0L)
  expect_length(object = result$out_of_frame_atg_positions, n = 0L)
  expect_length(object = result$out_of_frame_atg_frames, n = 0L)
})

test_that("GC content and GC3 match an independent count", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- profileCDS(cds = row_i$cds)

    expect_equal(
      object = result$gc_content,
      expected = gcFraction(sequence = row_i$cds),
      info = row_i$construct_id
    )
    expect_equal(
      object = result$gc3_content,
      expected = gcFraction(
        sequence = thirdPositions(sequence = row_i$cds)
      ),
      info = row_i$construct_id
    )
    expect_true(
      object = result$gc_content >= 0 && result$gc_content <= 1,
      info = row_i$construct_id
    )
    expect_true(
      object = result$gc3_content >= 0 && result$gc3_content <= 1,
      info = row_i$construct_id
    )
  }
})

test_that("GC content and GC3 can differ, and are read separately", {
  # ATG AAG AAG AAG AAG TAA. Every third base is a G except in the
  # terminal TAA, so GC3 is five sixths while overall GC is five
  # eighteenths: the two figures come apart.
  recoded <- profileCDS(cds = "ATGAAGAAGAAGAAGTAA")
  expect_equal(object = recoded$gc3_content, expected = 5 / 6)
  expect_equal(object = recoded$gc_content, expected = 5 / 18)
  expect_gt(
    object = recoded$gc3_content, expected = recoded$gc_content
  )

  # No G or C anywhere.
  flat <- profileCDS(cds = "ATAAAATTTAAATTATAA")
  expect_equal(object = flat$gc_content, expected = 0)
  expect_equal(object = flat$gc3_content, expected = 0)
})

test_that("structural facts match the dataset's coding sequences", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- profileCDS(cds = row_i$cds)

    expect_true(object = result$starts_with_atg, info = row_i$construct_id)
    expect_true(object = result$ends_with_stop, info = row_i$construct_id)
    expect_true(
      object = result$is_multiple_of_three, info = row_i$construct_id
    )
    expect_identical(
      object = result$cds_length,
      expected = nchar(x = row_i$cds),
      info = row_i$construct_id
    )
    expect_identical(
      object = result$n_codons,
      expected = nchar(x = row_i$cds) %/% 3L,
      info = row_i$construct_id
    )
  }
})

test_that("premature stop codons match the planted truth", {
  data(simulatedConstructs, package = "UTRunway")
  stop_codons <- c("TAA", "TAG", "TGA")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- profileCDS(cds = row_i$cds)

    expect_identical(
      object = result$has_premature_stop,
      expected = row_i$has_premature_stop,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$has_premature_stop,
      expected = length(x = result$premature_stop_positions) > 0L,
      info = row_i$construct_id
    )
    if (!result$has_premature_stop) {
      next
    }

    # Each reported position really holds a stop codon, in frame, and
    # before the final codon.
    expect_identical(
      object = substring(
        text = row_i$cds,
        first = result$premature_stop_positions,
        last = result$premature_stop_positions + 2L
      ),
      expected = result$premature_stop_codons,
      info = row_i$construct_id
    )
    expect_true(
      object = all(result$premature_stop_codons %in% stop_codons),
      info = row_i$construct_id
    )
    expect_true(
      object = all((result$premature_stop_positions - 1L) %% 3L == 0L),
      info = row_i$construct_id
    )
    expect_true(
      object = all(
        result$premature_stop_positions <= nchar(x = row_i$cds) - 5L
      ),
      info = row_i$construct_id
    )
  }
})

test_that("the terminal stop codon is not counted as premature", {
  # ATG AAA CCC TAA: the only stop is the terminator.
  result <- profileCDS(cds = "ATGAAACCCTAA")

  expect_true(object = result$ends_with_stop)
  expect_false(object = result$has_premature_stop)
  expect_length(object = result$premature_stop_positions, n = 0L)

  # A sequence that does not end in a stop still treats its last
  # codon as the terminator rather than a premature stop.
  no_terminator <- profileCDS(cds = "ATGAAACCCGGG")
  expect_false(object = no_terminator$ends_with_stop)
  expect_false(object = no_terminator$has_premature_stop)

  # Two premature stops are both reported, in order.
  two_stops <- profileCDS(cds = "ATGTAACCCTGAAAATAA")
  expect_identical(
    object = two_stops$premature_stop_positions, expected = c(4L, 10L)
  )
  expect_identical(
    object = two_stops$premature_stop_codons, expected = c("TAA", "TGA")
  )
})

test_that("out-of-frame ATGs are found and in-frame ones are not", {
  # ATG at position 1 is the start codon; the ATG at position 5 is
  # one base out of frame.
  result <- profileCDS(cds = "ATGCATGAAACCCGGTAA")

  expect_identical(object = result$n_out_of_frame_atg, expected = 1L)
  expect_identical(
    object = result$out_of_frame_atg_positions, expected = 5L
  )
  expect_identical(object = result$out_of_frame_atg_frames, expected = 1L)

  # An in-frame internal ATG is not out of frame, so it is not
  # reported here.
  in_frame <- profileCDS(cds = "ATGATGAAACCCGGGTAA")
  expect_identical(object = in_frame$n_out_of_frame_atg, expected = 0L)

  # Both offsets are distinguished.
  both <- profileCDS(cds = "ATGCATGCCATGAAACCCTAA")
  expect_true(
    object = all(both$out_of_frame_atg_frames %in% c(1L, 2L))
  )
})

test_that("out-of-frame ATG positions are consistent in the dataset", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- profileCDS(cds = row_i$cds)

    expect_identical(
      object = result$n_out_of_frame_atg,
      expected = length(x = result$out_of_frame_atg_positions),
      info = row_i$construct_id
    )
    if (result$n_out_of_frame_atg == 0L) {
      next
    }

    # Each position holds an ATG and is genuinely off the codon grid.
    expect_true(
      object = all(
        substring(
          text = row_i$cds,
          first = result$out_of_frame_atg_positions,
          last = result$out_of_frame_atg_positions + 2L
        ) == "ATG"
      ),
      info = row_i$construct_id
    )
    expect_identical(
      object = result$out_of_frame_atg_frames,
      expected = as.integer(
        x = (result$out_of_frame_atg_positions - 1L) %% 3L
      ),
      info = row_i$construct_id
    )
    expect_false(
      object = any(result$out_of_frame_atg_frames == 0L),
      info = row_i$construct_id
    )
  }
})

test_that("a length that is not a multiple of three warns", {
  expect_warning(
    object = profileCDS(cds = "ATGAAACCCTAAA"),
    regexp = "not a multiple of three"
  )
  result <- suppressWarnings(profileCDS(cds = "ATGAAACCCTAAA"))

  expect_false(object = result$is_multiple_of_three)
  expect_identical(object = result$cds_length, expected = 13L)
  # The trailing incomplete codon is ignored.
  expect_identical(object = result$n_codons, expected = 4L)
  expect_true(object = result$ends_with_stop)

  expect_no_warning(object = profileCDS(cds = "ATGAAACCCTAA"))
})

test_that("lower case and RNA input give the same result", {
  expect_identical(
    object = profileCDS(cds = "augaaauaacccgggtaa"),
    expected = profileCDS(cds = "ATGAAATAACCCGGGTAA")
  )
})

test_that("invalid input is rejected with an informative message", {
  expect_error(object = profileCDS(cds = NULL), regexp = "not NULL")
  expect_error(object = profileCDS(cds = NA), regexp = "must not be NA")
  expect_error(
    object = profileCDS(cds = ""), regexp = "must not be an empty string"
  )
  expect_error(
    object = profileCDS(cds = "ACGTN"), regexp = "unambiguous nucleotides"
  )
  expect_error(
    object = profileCDS(cds = 42), regexp = "must be a character string"
  )
  expect_error(
    object = profileCDS(cds = c("ATGTAA", "ATGTAA")),
    regexp = "must be a single sequence"
  )
  expect_error(
    object = profileCDS(cds = "AT"),
    regexp = "at least one complete codon"
  )
  # One complete codon is the smallest acceptable coding sequence.
  expect_silent(object = profileCDS(cds = "ATG"))
})
