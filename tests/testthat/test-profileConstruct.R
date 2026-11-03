# Profiling all twelve constructs takes a few seconds, so it is done
# once for the whole file.
data(simulatedConstructs, package = "UTRunway")
all_features <- profileConstruct(constructs = simulatedConstructs)

# The region summary features each region must always report. These
# are a contract: every plotting function and rankCandidates() reads
# the table by these names, so a rename has to fail here rather than
# quietly produce an empty plot downstream.
expected_scalar_features <- list(
  utr5 = c(
    "length", "gc_content", "n_upstream_aug", "n_upstream_orf",
    "n_out_of_frame_overlap", "kozak_score", "kozak_score_scaled",
    "kozak_strength_ordinal", "max_stem_nt", "max_window_gc",
    "n_motif_TOP", "n_motif_G_quadruplex"
  ),
  cds = c(
    "length", "gc_content", "gc3_content", "n_premature_stop",
    "n_out_of_frame_atg", "max_stem_nt", "max_window_gc",
    "n_motif_G_quadruplex"
  ),
  utr3 = c(
    "length", "gc_content", "n_are_pentamer", "n_are_clustered",
    "n_polya_signal", "n_polya_canonical", "n_mirna_seed",
    "max_stem_nt", "max_window_gc", "n_motif_G_quadruplex"
  )
)

# Positional feature names each region may report. A name outside
# this set means a new feature arrived without being documented.
expected_positional_features <- list(
  utr5 = c(
    "gc_window", "stem_window", "upstream_aug", "upstream_orf",
    "motif_TOP", "motif_G_quadruplex"
  ),
  cds = c(
    "gc_window", "stem_window", "premature_stop", "out_of_frame_atg",
    "motif_G_quadruplex"
  ),
  utr3 = c(
    "gc_window", "stem_window", "are_pentamer",
    "polya_signal_canonical", "polya_signal_variant",
    "motif_G_quadruplex"
  )
)

featureValue <- function(features, construct_id, feature) {
  features$value[
    features$construct_id == construct_id & features$feature == feature
  ]
}

test_that("the table has the documented columns and types", {
  expect_s3_class(object = all_features, class = "data.frame")
  expect_identical(
    object = names(x = all_features),
    expected = c(
      "construct_id", "region", "feature", "value", "position"
    )
  )
  expect_type(object = all_features$construct_id, type = "character")
  expect_type(object = all_features$region, type = "character")
  expect_type(object = all_features$feature, type = "character")
  expect_type(object = all_features$value, type = "double")
  expect_type(object = all_features$position, type = "integer")
  expect_gt(object = nrow(x = all_features), expected = 0)
})

test_that("every construct and region is represented", {
  expect_identical(
    object = sort(x = unique(x = all_features$construct_id)),
    expected = sort(x = simulatedConstructs$construct_id)
  )
  expect_identical(
    object = sort(x = unique(x = all_features$region)),
    expected = c("cds", "utr3", "utr5")
  )

  # Rows are grouped by construct in input order.
  first_appearance <- vapply(
    X = simulatedConstructs$construct_id,
    FUN = function(id) which(all_features$construct_id == id)[[1]],
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  )
  expect_identical(
    object = first_appearance, expected = sort(x = first_appearance)
  )
})

test_that("each region reports exactly its documented summaries", {
  for (region in names(x = expected_scalar_features)) {
    for (construct_id in simulatedConstructs$construct_id) {
      rows <- all_features[
        all_features$construct_id == construct_id &
          all_features$region == region &
          is.na(all_features$position),
      ]

      expect_identical(
        object = sort(x = rows$feature),
        expected = sort(x = expected_scalar_features[[region]]),
        info = paste(construct_id, region)
      )
      # One row per summary feature, never duplicated.
      expect_identical(
        object = anyDuplicated(x = rows$feature),
        expected = 0L,
        info = paste(construct_id, region)
      )
    }
  }
})

test_that("positional features stay within the documented set", {
  for (region in names(x = expected_positional_features)) {
    found <- unique(
      x = all_features$feature[
        all_features$region == region & !is.na(all_features$position)
      ]
    )
    unexpected <- setdiff(
      x = found,
      y = expected_positional_features[[region]]
    )

    # miRNA seed sites carry their family in the name, so they are
    # matched by prefix rather than listed.
    unexpected <- unexpected[
      !grepl(pattern = "^mirna_seed_", x = unexpected)
    ]
    expect_length(object = unexpected, n = 0L)
  }
})

test_that("summary features carry no position and positional ones do", {
  is_summary <- all_features$feature %in%
    unlist(x = expected_scalar_features, use.names = FALSE)

  expect_true(object = all(is.na(all_features$position[is_summary])))
  expect_false(object = anyNA(x = all_features$position[!is_summary]))
})

test_that("the coding sequence is separable from the scoreable rows", {
  scoreable <- all_features[all_features$region != "cds", ]
  cds_rows <- all_features[all_features$region == "cds", ]

  expect_gt(object = nrow(x = scoreable), expected = 0)
  expect_gt(object = nrow(x = cds_rows), expected = 0)
  expect_identical(
    object = nrow(x = scoreable) + nrow(x = cds_rows),
    expected = nrow(x = all_features)
  )
  expect_identical(
    object = sort(x = unique(x = scoreable$region)),
    expected = c("utr3", "utr5")
  )

  # The CDS-only features are absent once the CDS rows are dropped.
  expect_false(object = any(scoreable$feature == "gc3_content"))
  expect_false(object = any(scoreable$feature == "n_premature_stop"))
})

test_that("feature values match the planted truth", {
  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    id <- row_i$construct_id
    value <- function(feature) {
      featureValue(
        features = all_features, construct_id = id, feature = feature
      )
    }

    expect_identical(
      object = value(feature = "n_upstream_aug"),
      expected = as.numeric(x = row_i$n_upstream_aug),
      info = id
    )
    expect_identical(
      object = value(feature = "n_upstream_orf"),
      expected = as.numeric(x = row_i$n_upstream_orf),
      info = id
    )
    expect_identical(
      object = value(feature = "n_are_pentamer"),
      expected = as.numeric(x = row_i$n_are_pentamer),
      info = id
    )
    expect_identical(
      object = value(feature = "n_polya_canonical") > 0,
      expected = row_i$has_canonical_polya,
      info = id
    )
    expect_identical(
      object = value(feature = "n_motif_TOP") > 0,
      expected = row_i$has_top_motif,
      info = id
    )
    expect_identical(
      object = all_features$value[
        all_features$construct_id == id &
          all_features$region == "utr5" &
          all_features$feature == "n_motif_G_quadruplex"
      ] > 0,
      expected = row_i$has_g_quadruplex,
      info = id
    )
    expect_identical(
      object = value(feature = "n_premature_stop") > 0,
      expected = row_i$has_premature_stop,
      info = id
    )

    if (row_i$has_strong_kozak) {
      expect_identical(
        object = value(feature = "kozak_strength_ordinal"),
        expected = 2,
        info = id
      )
    }
    if (row_i$has_weak_kozak) {
      expect_identical(
        object = value(feature = "kozak_strength_ordinal"),
        expected = 0,
        info = id
      )
    }
  }
})

test_that("kozak_strength_ordinal matches scoreKozak's own verdict", {
  encoding <- c(weak = 0, adequate = 1, strong = 2)

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    direct <- suppressWarnings(scoreKozak(
      sequence = paste0(row_i$utr5, row_i$cds),
      position = nchar(x = row_i$utr5) + 1L
    ))

    expect_equal(
      object = featureValue(
        features = all_features,
        construct_id = row_i$construct_id,
        feature = "kozak_score"
      ),
      expected = direct$score,
      info = row_i$construct_id
    )
    expect_identical(
      object = featureValue(
        features = all_features,
        construct_id = row_i$construct_id,
        feature = "kozak_strength_ordinal"
      ),
      expected = unname(obj = encoding[[direct$strength]]),
      info = row_i$construct_id
    )
  }
})

test_that("length and GC content match the sequences profiled", {
  gcFraction <- function(sequence) {
    mean(strsplit(x = sequence, split = "")[[1]] %in% c("G", "C"))
  }

  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    for (region in c("utr5", "cds", "utr3")) {
      rows <- all_features[
        all_features$construct_id == row_i$construct_id &
          all_features$region == region,
      ]

      expect_identical(
        object = rows$value[rows$feature == "length"],
        expected = as.numeric(x = nchar(x = row_i[[region]])),
        info = paste(row_i$construct_id, region)
      )
      expect_equal(
        object = rows$value[rows$feature == "gc_content"],
        expected = gcFraction(sequence = row_i[[region]]),
        info = paste(row_i$construct_id, region)
      )
    }
  }
})

test_that("positions fall inside the region they belong to", {
  for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
    row_i <- simulatedConstructs[i, ]
    for (region in c("utr5", "cds", "utr3")) {
      positions <- all_features$position[
        all_features$construct_id == row_i$construct_id &
          all_features$region == region &
          !is.na(all_features$position)
      ]

      expect_true(
        object = all(positions >= 1L),
        info = paste(row_i$construct_id, region)
      )
      expect_true(
        object = all(positions <= nchar(x = row_i[[region]])),
        info = paste(row_i$construct_id, region)
      )
    }
  }
})

test_that("window rows follow the sliding window settings", {
  row_one <- simulatedConstructs[1, ]
  features <- profileConstruct(constructs = row_one)

  for (region in c("utr5", "cds", "utr3")) {
    starts <- features$position[
      features$region == region & features$feature == "gc_window"
    ]
    expected_starts <- seq(
      from = 1L,
      to = nchar(x = row_one[[region]]) - 40L + 1L,
      by = 10L
    )

    expect_identical(
      object = starts,
      expected = as.integer(x = expected_starts),
      info = region
    )
    # The two window profiles cover the same windows.
    expect_identical(
      object = features$position[
        features$region == region & features$feature == "stem_window"
      ],
      expected = starts,
      info = region
    )
  }
})

test_that("counts are reported as zero rather than left out", {
  # CONSTRUCT_01 has no TOP motif, no G-quadruplex and no premature
  # stop, so each count must be present and zero.
  for (feature in c("n_motif_TOP", "n_motif_G_quadruplex")) {
    value <- all_features$value[
      all_features$construct_id == "CONSTRUCT_01" &
        all_features$region == "utr5" &
        all_features$feature == feature
    ]
    expect_identical(object = value, expected = 0, info = feature)
  }
  expect_identical(
    object = featureValue(
      features = all_features,
      construct_id = "CONSTRUCT_01",
      feature = "n_premature_stop"
    ),
    expected = 0
  )
  # And no positional rows exist for a feature counted as zero.
  expect_identical(
    object = sum(all_features$feature == "motif_TOP" &
                   all_features$construct_id == "CONSTRUCT_01"),
    expected = 0L
  )
})

test_that("the three input forms give the same table", {
  from_frame <- profileConstruct(
    constructs = simulatedConstructs[1, ]
  )
  from_object <- profileConstruct(
    constructs = readConstruct(x = simulatedConstructs[1, ])
  )
  expect_identical(object = from_object, expected = from_frame)

  fasta_path <- system.file(
    "extdata", "exampleConstructs.fasta",
    package = "UTRunway"
  )
  skip_if(
    condition = !nzchar(x = fasta_path),
    message = "example FASTA is not installed"
  )
  from_fasta <- profileConstruct(constructs = fasta_path)
  expect_identical(
    object = sort(x = unique(x = from_fasta$construct_id)),
    expected = c("CONSTRUCT_01", "CONSTRUCT_06", "CONSTRUCT_12")
  )
  # The shared construct profiles identically from either source.
  expect_identical(
    object = from_fasta[from_fasta$construct_id == "CONSTRUCT_01", ],
    expected = from_frame
  )
})

test_that("arguments in ... reach the feature functions", {
  row_one <- simulatedConstructs[1, ]

  # A narrower window gives more windows over the same sequence.
  default_windows <- sum(
    all_features$feature == "gc_window" &
      all_features$construct_id == "CONSTRUCT_01"
  )
  narrow <- profileConstruct(
    constructs = row_one, window = 20, step = 20, min_stem = 4
  )
  expect_false(
    object = sum(narrow$feature == "gc_window") == default_windows
  )

  # Skipping the miRNA scan zeroes that count.
  no_seeds <- profileConstruct(
    constructs = simulatedConstructs[12, ], seeds = character(0)
  )
  expect_identical(
    object = no_seeds$value[no_seeds$feature == "n_mirna_seed"],
    expected = 0
  )
  expect_identical(
    object = sum(grepl(pattern = "^mirna_seed_", x = no_seeds$feature)),
    expected = 0L
  )

  # Supplied motifs replace the built-ins, including in the counts.
  supplied <- profileConstruct(
    constructs = row_one, motifs = c(mine = "ATTTA")
  )
  expect_true(object = any(supplied$feature == "n_motif_mine"))
  expect_false(object = any(supplied$feature == "n_motif_TOP"))
  expect_false(
    object = any(supplied$feature == "n_motif_G_quadruplex")
  )
})

test_that("a mistyped argument in ... is an error, not ignored", {
  row_one <- simulatedConstructs[1, ]

  expect_error(
    object = profileConstruct(constructs = row_one, windwo = 20),
    regexp = "Unknown argument `windwo`"
  )
  expect_error(
    object = profileConstruct(constructs = row_one, min_step = 3),
    regexp = "Unknown argument"
  )
  expect_error(
    object = profileConstruct(constructs = row_one, windwo = 20),
    regexp = "The names accepted are"
  )
  expect_error(
    object = profileConstruct(constructs = row_one, 40),
    regexp = "must be named"
  )
})

test_that("ambiguous constructs are refused by name", {
  # Regions are kept wider than the default 40 nucleotide window, so
  # that profiling the clean construct raises no window warnings.
  filler <- paste(rep(x = "C", times = 53), collapse = "")
  clean_utr5 <- paste0(filler, "CGCCACC")
  with_n <- suppressWarnings(readConstruct(
    x = data.frame(
      construct_id = c("clean", "ambiguous"),
      utr5 = c(clean_utr5, paste0(filler, "NGCCACC")),
      cds = rep(
        x = paste0(
          "ATG", paste(rep(x = "GCC", times = 19), collapse = ""), "TAA"
        ),
        times = 2
      ),
      utr3 = rep(
        x = paste0(
          paste(rep(x = "C", times = 50), collapse = ""), "AATAAACCCC"
        ),
        times = 2
      ),
      stringsAsFactors = FALSE
    )
  ))

  expect_error(
    object = profileConstruct(constructs = with_n),
    regexp = "\"ambiguous\""
  )
  expect_error(
    object = profileConstruct(constructs = with_n),
    regexp = "require unambiguous sequence"
  )

  # Dropping the flagged construct makes it profileable again.
  cleaned <- with_n[!with_n$has_ambiguous, ]
  expect_s3_class(
    object = profileConstruct(constructs = cleaned),
    class = "data.frame"
  )
})

test_that("invalid constructs input is rejected", {
  expect_error(
    object = profileConstruct(constructs = NULL), regexp = "not NULL"
  )
  expect_error(
    object = profileConstruct(constructs = 42),
    regexp = "must be a FASTA file path or a data frame"
  )
  expect_error(
    object = profileConstruct(
      constructs = data.frame(construct_id = "a", utr5 = "GCC")
    ),
    regexp = "missing"
  )
})
