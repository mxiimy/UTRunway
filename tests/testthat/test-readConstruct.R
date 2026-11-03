# Temporary directory for FASTA fixtures. Everything this file
# writes stays inside tempdir().
makeTestDir <- function() {
  path <- tempfile(pattern = "utrunway-test-")
  dir.create(path = path)
  path
}

goodConstructFrame <- function() {
  data.frame(
    construct_id = c("alpha", "beta"),
    utr5 = c("GCCACCGCCACC", "TTTTTTTTTTTT"),
    cds = c("ATGGCCTAA", "ATGGCCGGGTAG"),
    utr3 = c("AATAAACCCCC", "ATTAAACCCCC"),
    stringsAsFactors = FALSE
  )
}

test_that("a validated object has the documented class and columns", {
  result <- readConstruct(x = goodConstructFrame())

  expect_s3_class(object = result, class = "utrunway_constructs")
  expect_s3_class(object = result, class = "data.frame")
  expect_identical(
    object = class(x = result),
    expected = c("utrunway_constructs", "data.frame")
  )
  expect_identical(
    object = names(x = result),
    expected = c("construct_id", "utr5", "cds", "utr3", "has_ambiguous")
  )
  expect_identical(object = nrow(x = result), expected = 2L)
  expect_type(object = result$construct_id, type = "character")
  expect_type(object = result$utr5, type = "character")
  expect_type(object = result$has_ambiguous, type = "logical")
  expect_false(object = any(result$has_ambiguous))
})

test_that("the shipped FASTA round-trips to the dataset sequences", {
  data(simulatedConstructs, package = "UTRunway")
  fasta_path <- system.file(
    "extdata", "exampleConstructs.fasta",
    package = "UTRunway"
  )
  skip_if(
    condition = !nzchar(x = fasta_path),
    message = "example FASTA is not installed"
  )

  result <- readConstruct(x = fasta_path)

  expect_s3_class(object = result, class = "utrunway_constructs")
  expect_identical(object = nrow(x = result), expected = 3L)
  expect_identical(
    object = result$construct_id,
    expected = c("CONSTRUCT_01", "CONSTRUCT_06", "CONSTRUCT_12")
  )

  # Each region must match the dataset row it was written from.
  expected_rows <- simulatedConstructs[
    match(
      x = result$construct_id, table = simulatedConstructs$construct_id
    ),
  ]
  for (region in c("utr5", "cds", "utr3")) {
    expect_identical(
      object = result[[region]],
      expected = expected_rows[[region]],
      info = region
    )
  }
})

test_that("a data frame carries its extra columns through", {
  data(simulatedConstructs, package = "UTRunway")
  result <- readConstruct(x = simulatedConstructs)

  expect_s3_class(object = result, class = "utrunway_constructs")
  expect_identical(object = nrow(x = result), expected = 12L)
  # The required columns come first, then everything else.
  expect_identical(
    object = names(x = result)[1:5],
    expected = c("construct_id", "utr5", "cds", "utr3", "has_ambiguous")
  )
  expect_true(
    object = all(
      c("n_upstream_aug", "has_premature_stop", "n_are_pentamer") %in%
        names(x = result)
    )
  )
  expect_identical(
    object = result$n_upstream_aug,
    expected = simulatedConstructs$n_upstream_aug
  )
})

test_that("an already validated object is returned unchanged", {
  once <- readConstruct(x = goodConstructFrame())
  twice <- readConstruct(x = once)

  expect_identical(object = twice, expected = once)
})

test_that("lower case and RNA input are normalised to upper-case DNA", {
  result <- readConstruct(
    x = data.frame(
      construct_id = "rna_input",
      utr5 = "gccaccgccacc",
      cds = "auggccuaa",
      utr3 = "aauaaaccc",
      stringsAsFactors = FALSE
    )
  )

  expect_identical(object = result$utr5, expected = "GCCACCGCCACC")
  expect_identical(object = result$cds, expected = "ATGGCCTAA")
  expect_identical(object = result$utr3, expected = "AATAAACCC")
})

test_that("every problem is reported in one error", {
  broken <- data.frame(
    construct_id = c("bad_frame", "no_start", "no_stop"),
    utr5 = c("GCCACC", "GCCACC", "GCCACC"),
    cds = c("ATGGCCTA", "GGGGCCTAA", "ATGGCCGGG"),
    utr3 = c("AATAAA", "AATAAA", "AATAAA"),
    stringsAsFactors = FALSE
  )

  # One error, naming all three constructs, rather than stopping at
  # the first.
  expect_error(object = readConstruct(x = broken), regexp = "bad_frame")
  expect_error(object = readConstruct(x = broken), regexp = "no_start")
  expect_error(object = readConstruct(x = broken), regexp = "no_stop")
  expect_error(
    object = readConstruct(x = broken), regexp = "3 problems found"
  )

  # A single problem is reported in the singular.
  one_problem <- goodConstructFrame()
  one_problem$cds[[1]] <- "ATGGCCGGG"
  expect_error(
    object = readConstruct(x = one_problem), regexp = "1 problem found"
  )
})

test_that("each coding sequence rule is checked and named", {
  frame_error <- goodConstructFrame()
  frame_error$cds[[2]] <- "ATGGCCTA"
  expect_error(
    object = readConstruct(x = frame_error),
    regexp = "not a multiple of three"
  )

  start_error <- goodConstructFrame()
  start_error$cds[[1]] <- "GGGGCCTAA"
  expect_error(
    object = readConstruct(x = start_error),
    regexp = "rather than the start codon ATG"
  )

  stop_error <- goodConstructFrame()
  stop_error$cds[[1]] <- "ATGGCCGGG"
  expect_error(
    object = readConstruct(x = stop_error),
    regexp = "rather than a stop codon"
  )

  # All three stop codons are accepted.
  for (stop_codon in c("TAA", "TAG", "TGA")) {
    accepted <- goodConstructFrame()
    accepted$cds[[1]] <- paste0("ATGGCC", stop_codon)
    expect_s3_class(
      object = readConstruct(x = accepted),
      class = "utrunway_constructs"
    )
  }
})

test_that("a frame error does not also report a stop codon error", {
  # The last codon is undefined when the length does not divide by
  # three, so only the frame problem should be reported.
  frame_only <- goodConstructFrame()
  frame_only$cds[[1]] <- "ATGGCCTA"

  expect_error(
    object = readConstruct(x = frame_only), regexp = "1 problem found"
  )
  message <- tryCatch(
    expr = readConstruct(x = frame_only),
    error = function(e) conditionMessage(e)
  )
  expect_match(object = message, regexp = "not a multiple of three")
  expect_no_match(object = message, regexp = "stop codon")
})

test_that("bad alphabets and empty sequences are rejected by name", {
  bad_letter <- goodConstructFrame()
  bad_letter$utr5[[2]] <- "GCCA!CGCCACC"
  expect_error(
    object = readConstruct(x = bad_letter),
    regexp = "construct \"beta\", utr5 contains the letter !"
  )

  ambiguous_not_allowed <- goodConstructFrame()
  ambiguous_not_allowed$utr3[[1]] <- "AATAAACCCRR"
  expect_error(
    object = readConstruct(x = ambiguous_not_allowed),
    regexp = "only A, C, G, T, U and N are allowed"
  )

  empty_region <- goodConstructFrame()
  empty_region$utr3[[2]] <- ""
  expect_error(
    object = readConstruct(x = empty_region),
    regexp = "construct \"beta\", utr3 is empty"
  )
})

test_that("identifiers must be present and unique", {
  duplicated_ids <- goodConstructFrame()
  duplicated_ids$construct_id <- c("same", "same")
  expect_error(
    object = readConstruct(x = duplicated_ids),
    regexp = "construct_id must be unique"
  )

  empty_id <- goodConstructFrame()
  empty_id$construct_id[[2]] <- ""
  expect_error(
    object = readConstruct(x = empty_id),
    regexp = "missing or empty construct_id"
  )
})

test_that("an ambiguous base warns and is flagged, not rejected", {
  with_n <- goodConstructFrame()
  with_n$utr5[[1]] <- "GCCANCGCCACC"

  expect_warning(
    object = readConstruct(x = with_n), regexp = "Ambiguous base N"
  )
  expect_warning(
    object = readConstruct(x = with_n), regexp = "\"alpha\""
  )

  result <- suppressWarnings(readConstruct(x = with_n))
  expect_s3_class(object = result, class = "utrunway_constructs")
  expect_identical(
    object = result$has_ambiguous, expected = c(TRUE, FALSE)
  )

  # An N in any region flags the construct.
  in_cds <- goodConstructFrame()
  in_cds$cds[[2]] <- "ATGGNCGGGTAG"
  flagged <- suppressWarnings(readConstruct(x = in_cds))
  expect_identical(
    object = flagged$has_ambiguous, expected = c(FALSE, TRUE)
  )

  expect_no_warning(object = readConstruct(x = goodConstructFrame()))
})

test_that("missing data frame columns are named", {
  expect_error(
    object = readConstruct(
      x = data.frame(construct_id = "a", utr5 = "GCC", cds = "ATGTAA")
    ),
    regexp = "`utr3`"
  )
  expect_error(
    object = readConstruct(x = data.frame(utr5 = "GCC")),
    regexp = "these are missing"
  )
  expect_error(
    object = readConstruct(
      x = data.frame(
        construct_id = character(0), utr5 = character(0),
        cds = character(0), utr3 = character(0)
      )
    ),
    regexp = "at least one construct"
  )
})

test_that("malformed FASTA input is rejected informatively", {
  temporary <- makeTestDir()
  on.exit(expr = unlink(x = temporary, recursive = TRUE), add = TRUE)

  # A file with no header line at all has no records.
  no_records <- file.path(temporary, "no-records.fasta")
  writeLines(text = c("ACGT", "ACGT"), con = no_records)
  expect_error(
    object = readConstruct(x = no_records),
    regexp = "No FASTA records found"
  )

  empty_file <- file.path(temporary, "empty.fasta")
  writeLines(text = character(0), con = empty_file)
  expect_error(
    object = readConstruct(x = empty_file),
    regexp = "No FASTA records found"
  )

  # Sequence before the first header is a different fault: there are
  # records, but the file does not begin with one.
  leading_sequence <- file.path(temporary, "leading-sequence.fasta")
  writeLines(
    text = c(
      "ACGTACGT",
      ">CONSTRUCT_01|region=utr5", "GCCACC",
      ">CONSTRUCT_01|region=cds", "ATGGCCTAA",
      ">CONSTRUCT_01|region=utr3", "AATAAA"
    ),
    con = leading_sequence
  )
  expect_error(
    object = readConstruct(x = leading_sequence),
    regexp = "not a FASTA header"
  )

  no_region <- file.path(temporary, "no-region.fasta")
  writeLines(text = c(">CONSTRUCT_01", "ACGT"), con = no_region)
  expect_error(
    object = readConstruct(x = no_region),
    regexp = "does not name a region"
  )

  bad_region <- file.path(temporary, "bad-region.fasta")
  writeLines(
    text = c(">CONSTRUCT_01|region=promoter", "ACGT"), con = bad_region
  )
  expect_error(
    object = readConstruct(x = bad_region),
    regexp = "which is not one of"
  )

  missing_region <- file.path(temporary, "missing-region.fasta")
  writeLines(
    text = c(
      ">CONSTRUCT_01|region=utr5", "GCCACC",
      ">CONSTRUCT_01|region=cds", "ATGGCCTAA"
    ),
    con = missing_region
  )
  expect_error(
    object = readConstruct(x = missing_region),
    regexp = "has no utr3 record"
  )

  repeated_region <- file.path(temporary, "repeated-region.fasta")
  writeLines(
    text = c(
      ">CONSTRUCT_01|region=utr5", "GCCACC",
      ">CONSTRUCT_01|region=utr5", "GCCACC",
      ">CONSTRUCT_01|region=cds", "ATGGCCTAA",
      ">CONSTRUCT_01|region=utr3", "AATAAA"
    ),
    con = repeated_region
  )
  expect_error(
    object = readConstruct(x = repeated_region),
    regexp = "more than one utr5 record"
  )
})

test_that("both header forms are accepted and sequences are joined", {
  temporary <- makeTestDir()
  on.exit(expr = unlink(x = temporary, recursive = TRUE), add = TRUE)
  short_form <- file.path(temporary, "short-form.fasta")
  writeLines(
    text = c(
      ">CONSTRUCT_01|utr5", "GCCACC", "GCCACC",
      ">CONSTRUCT_01|region=cds", "ATGGCC", "TAA",
      ">CONSTRUCT_01|UTR3", "AATAAA"
    ),
    con = short_form
  )

  result <- readConstruct(x = short_form)
  expect_identical(object = nrow(x = result), expected = 1L)
  # Sequence lines are concatenated, and the region name is matched
  # without regard to case.
  expect_identical(object = result$utr5, expected = "GCCACCGCCACC")
  expect_identical(object = result$cds, expected = "ATGGCCTAA")
  expect_identical(object = result$utr3, expected = "AATAAA")
})

test_that("invalid x is rejected with an informative message", {
  expect_error(object = readConstruct(x = NULL), regexp = "not NULL")
  expect_error(
    object = readConstruct(x = 42),
    regexp = "must be a FASTA file path or a data frame"
  )
  expect_error(
    object = readConstruct(x = c("one.fasta", "two.fasta")),
    regexp = "must be a single FASTA file path"
  )
  expect_error(
    object = readConstruct(x = ""), regexp = "not an empty string"
  )
  expect_error(
    object = readConstruct(x = "no-such-file.fasta"),
    regexp = "not an existing file"
  )
})
