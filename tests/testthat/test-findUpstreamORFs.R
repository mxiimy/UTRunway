expected_uorf_columns <- c(
  "aug_position", "distance_to_cds", "frame", "in_frame_with_cds",
  "opens_uorf", "stop_position", "uorf_length_nt", "extends_into_cds",
  "out_of_frame_overlap", "kozak_score", "kozak_strength"
)

test_that("the result has the documented columns and types", {
  result <- findUpstreamORFs(utr5 = "GCCACCATGAAATAACCGGTTACGTACGT")

  expect_s3_class(object = result, class = "data.frame")
  expect_identical(
    object = names(x = result), expected = expected_uorf_columns
  )
  expect_type(object = result$aug_position, type = "integer")
  expect_type(object = result$distance_to_cds, type = "integer")
  expect_type(object = result$frame, type = "integer")
  expect_type(object = result$in_frame_with_cds, type = "logical")
  expect_type(object = result$opens_uorf, type = "logical")
  expect_type(object = result$stop_position, type = "integer")
  expect_type(object = result$uorf_length_nt, type = "integer")
  expect_type(object = result$extends_into_cds, type = "logical")
  expect_type(object = result$out_of_frame_overlap, type = "logical")
  expect_type(object = result$kozak_score, type = "double")
  expect_type(object = result$kozak_strength, type = "character")
})

test_that("a 5' UTR with no AUG gives a zero-row data frame", {
  result <- findUpstreamORFs(utr5 = "CCCCACCCCACCCCACCCCA")

  expect_false(object = is.null(x = result))
  expect_s3_class(object = result, class = "data.frame")
  expect_identical(object = nrow(x = result), expected = 0L)
  expect_identical(
    object = names(x = result), expected = expected_uorf_columns
  )

  # The empty result must carry the same types as a populated one, so
  # that callers can rbind or subset either without special-casing.
  populated <- findUpstreamORFs(
    utr5 = "GCCACCATGAAATAACCGGTTACGTACGT"
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

test_that("upstream AUG and uORF counts match the planted truth", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- findUpstreamORFs(utr5 = row_i$utr5, cds = row_i$cds)

    expect_identical(
      object = nrow(x = result),
      expected = row_i$n_upstream_aug,
      info = row_i$construct_id
    )
    expect_identical(
      object = sum(result$opens_uorf),
      expected = row_i$n_upstream_orf,
      info = row_i$construct_id
    )
  }
})

test_that("positions, frames and distances are internally consistent", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- findUpstreamORFs(utr5 = row_i$utr5, cds = row_i$cds)
    if (nrow(x = result) == 0L) {
      next
    }

    # Each reported position really is an ATG.
    expect_true(
      object = all(
        substring(
          text = row_i$utr5,
          first = result$aug_position,
          last = result$aug_position + 2L
        ) == "ATG"
      ),
      info = row_i$construct_id
    )

    # Rows come back in positional order.
    expect_identical(
      object = result$aug_position,
      expected = sort(x = result$aug_position),
      info = row_i$construct_id
    )

    expect_identical(
      object = result$distance_to_cds,
      expected = nchar(x = row_i$utr5) + 1L - result$aug_position,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$frame,
      expected = as.integer(x = result$distance_to_cds %% 3L),
      info = row_i$construct_id
    )
    expect_identical(
      object = result$in_frame_with_cds,
      expected = result$frame == 0L,
      info = row_i$construct_id
    )

    # The three cases are mutually exclusive and exhaustive.
    expect_identical(
      object = result$extends_into_cds,
      expected = !result$opens_uorf,
      info = row_i$construct_id
    )
    expect_identical(
      object = result$out_of_frame_overlap,
      expected = !result$opens_uorf & result$frame != 0L,
      info = row_i$construct_id
    )
  }
})

test_that("a reported uORF really ends at an in-frame stop codon", {
  data(simulatedConstructs, package = "UTRunway")
  stop_codons <- c("TAA", "TAG", "TGA")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    result <- findUpstreamORFs(utr5 = row_i$utr5, cds = row_i$cds)
    uorfs <- result[result$opens_uorf, ]
    if (nrow(x = uorfs) == 0L) {
      next
    }

    expect_true(
      object = all(
        substring(
          text = row_i$utr5,
          first = uorfs$stop_position,
          last = uorfs$stop_position + 2L
        ) %in% stop_codons
      ),
      info = row_i$construct_id
    )
    # The stop codon is in frame with its own AUG, and lies wholly
    # inside the 5' UTR.
    expect_true(
      object = all((uorfs$stop_position - uorfs$aug_position) %% 3L == 0L),
      info = row_i$construct_id
    )
    expect_true(
      object = all(uorfs$stop_position + 2L <= nchar(x = row_i$utr5)),
      info = row_i$construct_id
    )
    expect_identical(
      object = uorfs$uorf_length_nt,
      expected = uorfs$stop_position + 3L - uorfs$aug_position,
      info = row_i$construct_id
    )
    expect_true(
      object = all(uorfs$uorf_length_nt %% 3L == 0L),
      info = row_i$construct_id
    )

    # An AUG that does not terminate inside the UTR has no stop.
    overlapping <- result[!result$opens_uorf, ]
    expect_true(
      object = all(is.na(x = overlapping$stop_position)),
      info = row_i$construct_id
    )
    expect_true(
      object = all(is.na(x = overlapping$uorf_length_nt)),
      info = row_i$construct_id
    )
  }
})

test_that("overlapping and nested starts each get their own row", {
  # Two AUGs, at 1 and 4, in the same frame and sharing the stop codon
  # at position 10. The second uORF is nested inside the first.
  result <- findUpstreamORFs(utr5 = "ATGATGAAATAACCGGTTACGTACGTACG")

  expect_identical(object = nrow(x = result), expected = 2L)
  expect_identical(object = result$aug_position, expected = c(1L, 4L))
  expect_true(object = all(result$opens_uorf))
  expect_identical(object = result$stop_position, expected = c(10L, 10L))
  expect_identical(object = result$uorf_length_nt, expected = c(12L, 9L))

  # Out-of-frame overlap: a second AUG one base along is in a
  # different frame and reaches a different stop codon.
  staggered <- findUpstreamORFs(utr5 = "ATGCATGAAACTAAGTAACCGGTTACGTA")
  expect_identical(object = nrow(x = staggered), expected = 2L)
  expect_identical(object = staggered$aug_position, expected = c(1L, 5L))
  expect_false(
    object = staggered$frame[[1]] == staggered$frame[[2]]
  )
})

test_that("an AUG with no in-frame stop extends into the CDS", {
  # ATG at position 7 of a 12 nt UTR: 6 nt from the main start codon,
  # so in frame with it, and no stop codon inside the UTR.
  in_frame <- findUpstreamORFs(utr5 = "CCCCCCATGCCC")
  expect_identical(object = in_frame$aug_position, expected = 7L)
  expect_identical(object = in_frame$distance_to_cds, expected = 6L)
  expect_identical(object = in_frame$frame, expected = 0L)
  expect_true(object = in_frame$in_frame_with_cds)
  expect_false(object = in_frame$opens_uorf)
  expect_true(object = in_frame$extends_into_cds)
  # In frame with the main start, so it is an N-terminal extension
  # rather than an out-of-frame overlap.
  expect_false(object = in_frame$out_of_frame_overlap)

  # One more base at the 3' end shifts it out of frame.
  out_of_frame <- findUpstreamORFs(utr5 = "CCCCCCATGCCCC")
  expect_identical(object = out_of_frame$frame, expected = 1L)
  expect_false(object = out_of_frame$in_frame_with_cds)
  expect_true(object = out_of_frame$extends_into_cds)
  expect_true(object = out_of_frame$out_of_frame_overlap)
})

test_that("a stop codon just past the UTR end does not count", {
  # TAA straddles the UTR boundary, so the uORF does not terminate
  # inside the 5' UTR and the AUG overlaps the main start instead.
  result <- findUpstreamORFs(utr5 = "CCCCCCATGAAAT", cds = "AAATGCCCTAA")

  expect_identical(object = nrow(x = result), expected = 1L)
  expect_false(object = result$opens_uorf)
  expect_identical(object = result$stop_position, expected = NA_integer_)
  expect_true(object = result$extends_into_cds)
})

test_that("the CDS only completes the Kozak context window", {
  data(simulatedConstructs, package = "UTRunway")
  row_i <- simulatedConstructs[
    simulatedConstructs$construct_id == "CONSTRUCT_03",
  ]

  with_cds <- findUpstreamORFs(utr5 = row_i$utr5, cds = row_i$cds)
  without_cds <- findUpstreamORFs(utr5 = row_i$utr5)

  # Which AUGs and uORFs are found must not depend on the CDS.
  feature_columns <- setdiff(
    x = expected_uorf_columns,
    y = c("kozak_score", "kozak_strength")
  )
  expect_identical(
    object = with_cds[, feature_columns],
    expected = without_cds[, feature_columns]
  )

  # An AUG flush against the end of the 5' UTR has no +4 base within
  # the UTR, so its context window can only be completed from the CDS.
  flush_utr5 <- "GCCACCGCCACCATG"
  scored_with_cds <- findUpstreamORFs(
    utr5 = flush_utr5, cds = "GCCCTGTAA"
  )
  scored_alone <- findUpstreamORFs(utr5 = flush_utr5)

  # The AUG occupies the last three bases of the UTR, so its A sits
  # three nucleotides before the A of the main start codon.
  expect_identical(object = scored_with_cds$distance_to_cds, expected = 3L)
  expect_identical(object = scored_with_cds$frame, expected = 0L)
  expect_false(object = is.na(x = scored_with_cds$kozak_strength))
  expect_identical(
    object = scored_with_cds$kozak_strength, expected = "strong"
  )
  expect_true(object = is.na(x = scored_alone$kozak_strength))
  expect_true(object = is.na(x = scored_alone$kozak_score))

  # The AUG itself is found either way.
  expect_identical(object = scored_alone$aug_position, expected = 13L)
})

test_that("Kozak columns are NA near the 5' end, without warning", {
  # The AUG is at position 1, so it has no context window at all.
  expect_no_warning(
    object = findUpstreamORFs(utr5 = "ATGAAATAACCGGTTACGTACGTACG")
  )
  result <- findUpstreamORFs(utr5 = "ATGAAATAACCGGTTACGTACGTACG")
  expect_identical(object = result$kozak_score, expected = NA_real_)
  expect_identical(
    object = result$kozak_strength, expected = NA_character_
  )

  # With a full window the columns are populated. Here +4 is a G and
  # -3 a purine, so this upstream start is a strong one.
  scored <- findUpstreamORFs(utr5 = "GCCACCATGGAATAACCGGTTACGTACGT")
  expect_false(object = is.na(x = scored$kozak_score))
  expect_identical(object = scored$kozak_strength, expected = "strong")
  expect_true(object = scored$opens_uorf)
})

test_that("lower case and RNA input give the same result", {
  expect_identical(
    object = findUpstreamORFs(utr5 = "augaaauaaccgguuacguacguacg"),
    expected = findUpstreamORFs(utr5 = "ATGAAATAACCGGTTACGTACGTACG")
  )
})

test_that("invalid input is rejected with an informative message", {
  expect_error(
    object = findUpstreamORFs(utr5 = NULL), regexp = "not NULL"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = ""),
    regexp = "must not be an empty string"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = NA), regexp = "must not be NA"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = 42),
    regexp = "must be a character string"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = c("ATGC", "ATGC")),
    regexp = "must be a single sequence"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = "ATGXYZ"),
    regexp = "unambiguous nucleotides"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = "ATGCCC", cds = 7),
    regexp = "`cds` must be a character string"
  )
  expect_error(
    object = findUpstreamORFs(utr5 = "ATGCCC", cds = "ATGXXX"),
    regexp = "`cds` must contain only unambiguous nucleotides"
  )
})
