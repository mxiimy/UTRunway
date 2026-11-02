# The planted features in `simulatedConstructs` are re-derived here
# straight from the stored sequences and compared with the recorded
# truth columns. These checks are deliberately independent of
# data-raw/makeSimulatedConstructs.R: they are written against the
# definition of each feature rather than by reusing the generator's
# helpers, so that the committed dataset and the truth it claims
# cannot drift apart unnoticed.

stop_codons <- c("TAA", "TAG", "TGA")

# Count occurrences of `motif`, counting overlapping hits separately.
countMotifHits <- function(sequence, motif) {
  n_starts <- nchar(x = sequence) - nchar(x = motif) + 1L
  if (n_starts < 1L) {
    return(0L)
  }
  starts <- seq_len(length.out = n_starts)
  sum(
    substring(
      text = sequence,
      first = starts,
      last = starts + nchar(x = motif) - 1L
    ) == motif
  )
}

# Split a sequence into codons from `from`, dropping any trailing
# partial codon.
splitCodons <- function(sequence, from = 1L) {
  last_start <- nchar(x = sequence) - 2L
  if (from > last_start) {
    return(character(length = 0L))
  }
  starts <- seq(from = from, to = last_start, by = 3L)
  substring(text = sequence, first = starts, last = starts + 2L)
}

# Every 3-mer of `sequence`, indexed by its start position.
splitCodonsAtEveryOffset <- function(sequence) {
  n_starts <- nchar(x = sequence) - 2L
  if (n_starts < 1L) {
    return(character(length = 0L))
  }
  starts <- seq_len(length.out = n_starts)
  substring(text = sequence, first = starts, last = starts + 2L)
}

# Upstream AUGs whose in-frame stop codon also lies inside the 5' UTR.
countCompleteUorfs <- function(utr5) {
  aug_starts <- which(
    splitCodonsAtEveryOffset(sequence = utr5) == "ATG"
  )
  n_orf <- 0L
  for (start_i in aug_starts) {
    codons <- splitCodons(sequence = utr5, from = start_i)
    if (any(codons %in% stop_codons)) {
      n_orf <- n_orf + 1L
    }
  }
  n_orf
}

test_that("the dataset has the expected shape and column types", {
  data(simulatedConstructs, package = "UTRunway")

  expect_s3_class(object = simulatedConstructs, class = "data.frame")
  expect_identical(object = nrow(x = simulatedConstructs), expected = 12L)
  expect_identical(
    object = names(x = simulatedConstructs),
    expected = c(
      "construct_id", "utr5", "cds", "utr3",
      "n_upstream_aug", "n_upstream_orf",
      "has_strong_kozak", "has_weak_kozak",
      "has_g_quadruplex", "has_top_motif",
      "n_are_pentamer", "has_canonical_polya",
      "has_premature_stop"
    )
  )

  for (column in c("construct_id", "utr5", "cds", "utr3")) {
    expect_type(object = simulatedConstructs[[column]], type = "character")
  }
  for (column in c(
    "n_upstream_aug", "n_upstream_orf", "n_are_pentamer"
  )) {
    expect_type(object = simulatedConstructs[[column]], type = "integer")
  }
  for (column in c(
    "has_strong_kozak", "has_weak_kozak", "has_g_quadruplex",
    "has_top_motif", "has_canonical_polya", "has_premature_stop"
  )) {
    expect_type(object = simulatedConstructs[[column]], type = "logical")
  }

  expect_false(object = anyNA(x = simulatedConstructs))
  expect_identical(
    object = anyDuplicated(x = simulatedConstructs$construct_id),
    expected = 0L
  )
})

test_that("sequences are unambiguous DNA of the documented lengths", {
  data(simulatedConstructs, package = "UTRunway")

  all_sequences <- c(
    simulatedConstructs$utr5,
    simulatedConstructs$cds,
    simulatedConstructs$utr3
  )
  expect_false(
    object = any(grepl(pattern = "[^ACGT]", x = all_sequences))
  )

  utr5_nt <- nchar(x = simulatedConstructs$utr5)
  cds_nt <- nchar(x = simulatedConstructs$cds)
  utr3_nt <- nchar(x = simulatedConstructs$utr3)

  expect_true(object = all(utr5_nt >= 80L & utr5_nt <= 200L))
  expect_true(object = all(cds_nt >= 300L & cds_nt <= 600L))
  expect_true(object = all(utr3_nt >= 200L & utr3_nt <= 500L))
})

test_that("every CDS is a valid reading frame", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    cds_i <- simulatedConstructs$cds[[i]]
    codons <- splitCodons(sequence = cds_i)

    expect_identical(
      object = nchar(x = cds_i) %% 3L,
      expected = 0L,
      info = simulatedConstructs$construct_id[[i]]
    )
    expect_identical(
      object = codons[[1]],
      expected = "ATG",
      info = simulatedConstructs$construct_id[[i]]
    )
    expect_true(
      object = codons[[length(x = codons)]] %in% stop_codons,
      info = simulatedConstructs$construct_id[[i]]
    )
  }
})

test_that("planted premature in-frame stops match the truth column", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    codons <- splitCodons(sequence = simulatedConstructs$cds[[i]])
    internal_codons <- codons[-length(x = codons)]

    expect_identical(
      object = any(internal_codons %in% stop_codons),
      expected = simulatedConstructs$has_premature_stop[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
  }
})

test_that("upstream AUG and uORF counts match the truth columns", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    utr5_i <- simulatedConstructs$utr5[[i]]

    expect_identical(
      object = countMotifHits(sequence = utr5_i, motif = "ATG"),
      expected = simulatedConstructs$n_upstream_aug[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
    expect_identical(
      object = countCompleteUorfs(utr5 = utr5_i),
      expected = simulatedConstructs$n_upstream_orf[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
  }

  # A uORF needs an upstream AUG, so uORFs can never outnumber them.
  expect_true(
    object = all(
      simulatedConstructs$n_upstream_orf <=
        simulatedConstructs$n_upstream_aug
    )
  )
})

test_that("Kozak context columns match the sequence around the start", {
  data(simulatedConstructs, package = "UTRunway")

  # A construct cannot be both strong and weak at once.
  expect_false(
    object = any(
      simulatedConstructs$has_strong_kozak &
        simulatedConstructs$has_weak_kozak
    )
  )

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    utr5_i <- simulatedConstructs$utr5[[i]]
    cds_i <- simulatedConstructs$cds[[i]]

    # Position -3 relative to the A of the AUG, and position +4,
    # which is the first base of the second codon.
    base_minus_3 <- substr(
      x = utr5_i,
      start = nchar(x = utr5_i) - 2L,
      stop = nchar(x = utr5_i) - 2L
    )
    base_plus_4 <- substr(x = cds_i, start = 4L, stop = 4L)

    if (simulatedConstructs$has_strong_kozak[[i]]) {
      expect_true(
        object = base_minus_3 %in% c("A", "G"),
        info = simulatedConstructs$construct_id[[i]]
      )
      expect_identical(
        object = base_plus_4,
        expected = "G",
        info = simulatedConstructs$construct_id[[i]]
      )
      expect_true(
        object = grepl(
          pattern = "GCCACCATGG",
          x = paste0(utr5_i, cds_i),
          fixed = TRUE
        ),
        info = simulatedConstructs$construct_id[[i]]
      )
    }

    if (simulatedConstructs$has_weak_kozak[[i]]) {
      expect_true(
        object = base_minus_3 %in% c("C", "T"),
        info = simulatedConstructs$construct_id[[i]]
      )
      expect_false(
        object = base_plus_4 == "G",
        info = simulatedConstructs$construct_id[[i]]
      )
    }
  }
})

test_that("structure and TOP motif columns match the 5' UTR", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    utr5_i <- simulatedConstructs$utr5[[i]]

    # Four G-triplets, the arrangement scored as quadruplex-prone.
    has_g_run <- grepl(
      pattern = "G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",
      x = utr5_i
    )
    expect_identical(
      object = has_g_run,
      expected = simulatedConstructs$has_g_quadruplex[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )

    # A TOP motif is a C followed by a pyrimidine run at the 5' end.
    has_top <- grepl(pattern = "^C[CT]{4,}", x = utr5_i)
    expect_identical(
      object = has_top,
      expected = simulatedConstructs$has_top_motif[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
  }
})

test_that("3' UTR element columns match the sequence", {
  data(simulatedConstructs, package = "UTRunway")

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    utr3_i <- simulatedConstructs$utr3[[i]]

    # Overlapping AUUUA pentamers are counted separately.
    expect_identical(
      object = countMotifHits(sequence = utr3_i, motif = "ATTTA"),
      expected = simulatedConstructs$n_are_pentamer[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
    expect_identical(
      object = countMotifHits(sequence = utr3_i, motif = "AATAAA") > 0L,
      expected = simulatedConstructs$has_canonical_polya[[i]],
      info = simulatedConstructs$construct_id[[i]]
    )
  }
})

test_that("planted features are spread as the dataset design states", {
  data(simulatedConstructs, package = "UTRunway")

  expect_identical(
    object = sum(simulatedConstructs$n_upstream_aug > 0L),
    expected = 3L
  )
  expect_identical(
    object = sum(simulatedConstructs$has_strong_kozak),
    expected = 3L
  )
  expect_identical(
    object = sum(simulatedConstructs$has_weak_kozak),
    expected = 3L
  )
  expect_identical(
    object = sum(simulatedConstructs$has_g_quadruplex),
    expected = 2L
  )
  expect_identical(
    object = sum(simulatedConstructs$has_top_motif),
    expected = 2L
  )
  expect_identical(
    object = sum(simulatedConstructs$n_are_pentamer > 0L),
    expected = 3L
  )
  expect_identical(
    object = sum(!simulatedConstructs$has_canonical_polya),
    expected = 2L
  )
  expect_identical(
    object = sum(simulatedConstructs$has_premature_stop),
    expected = 2L
  )
})

test_that("the shipped FASTA agrees with the dataset", {
  data(simulatedConstructs, package = "UTRunway")

  fasta_path <- system.file(
    "extdata", "exampleConstructs.fasta",
    package = "UTRunway"
  )
  skip_if(
    condition = !nzchar(x = fasta_path),
    message = "example FASTA is not installed"
  )

  fasta_lines <- readLines(con = fasta_path)
  header_at <- grep(pattern = "^>", x = fasta_lines)
  expect_identical(object = length(x = header_at), expected = 9L)

  block_end <- c(header_at[-1] - 1L, length(x = fasta_lines))

  for (i in seq_along(along.with = header_at)) {
    header <- sub(
      pattern = "^>", replacement = "", x = fasta_lines[[header_at[[i]]]]
    )
    parts <- strsplit(x = header, split = "|", fixed = TRUE)[[1]]
    construct_id <- parts[[1]]
    region <- sub(pattern = "^region=", replacement = "", x = parts[[2]])

    expect_true(object = region %in% c("utr5", "cds", "utr3"))
    expect_true(
      object = construct_id %in% simulatedConstructs$construct_id
    )

    sequence <- paste(
      fasta_lines[(header_at[[i]] + 1L):block_end[[i]]],
      collapse = ""
    )
    expect_identical(
      object = sequence,
      expected = simulatedConstructs[[region]][
        simulatedConstructs$construct_id == construct_id
      ],
      info = header
    )
  }
})
