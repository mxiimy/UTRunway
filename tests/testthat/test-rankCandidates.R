# Profiling is the slow step, so it is done once for the whole file.
data(simulatedConstructs, package = "UTRunway")
all_features <- profileConstruct(constructs = simulatedConstructs)
six_features <- all_features[
  all_features$construct_id %in% simulatedConstructs$construct_id[1:6],
]

# The scored feature set and its directions are a contract: the
# ranking's meaning depends on them, so a change has to be deliberate
# and show up here rather than silently reordering candidates.
expected_directions <- data.frame(
  region = c(
    rep(x = "utr5", times = 7L), rep(x = "utr3", times = 5L)
  ),
  feature = c(
    "n_upstream_aug", "n_upstream_orf", "n_out_of_frame_overlap",
    "kozak_score_scaled", "gc_content", "max_stem_nt",
    "n_motif_G_quadruplex",
    "gc_content", "n_are_pentamer", "n_are_clustered",
    "n_polya_canonical", "n_mirna_seed"
  ),
  direction = c(-1L, -1L, -1L, 1L, -1L, -1L, -1L, 1L, -1L, -1L, 1L, -1L),
  stringsAsFactors = FALSE
)

test_that("the direction table is the documented set", {
  directions <- featureDirections()

  expect_s3_class(object = directions, class = "data.frame")
  expect_identical(
    object = directions[, c("region", "feature", "direction")],
    expected = expected_directions
  )
  # Every direction carries a stated reason, so no judgement is
  # made silently.
  expect_true(object = all(nzchar(x = directions$rationale)))
  expect_true(object = all(directions$direction %in% c(-1L, 1L)))

  # GC content is scored in both UTRs, in opposite directions.
  gc_rows <- directions[directions$feature == "gc_content", ]
  expect_identical(object = nrow(x = gc_rows), expected = 2L)
  expect_identical(
    object = gc_rows$direction[gc_rows$region == "utr5"], expected = -1L
  )
  expect_identical(
    object = gc_rows$direction[gc_rows$region == "utr3"], expected = 1L
  )
})

test_that("features with no asserted direction are not scored", {
  directions <- featureDirections()

  # Duplicate representations of a quantity already scored, which
  # would otherwise count it two or three times over.
  expect_false(object = "kozak_score" %in% directions$feature)
  expect_false(
    object = "kozak_strength_ordinal" %in% directions$feature
  )
  expect_false(object = "max_window_gc" %in% directions$feature)

  # Features whose direction is a design choice, not a quality.
  expect_false(object = "length" %in% directions$feature)
  expect_false(object = "n_motif_TOP" %in% directions$feature)
  expect_false(object = "n_polya_signal" %in% directions$feature)

  # 3' UTR structure, where neither direction is safe to assert.
  utr3_features <- directions$feature[directions$region == "utr3"]
  expect_false(object = "max_stem_nt" %in% utr3_features)
  expect_false(object = "n_motif_G_quadruplex" %in% utr3_features)
})

test_that("the returned object has the documented shape", {
  ranking <- rankCandidates(feature_table = six_features)

  expect_s3_class(object = ranking, class = "utrunway_ranking")
  expect_identical(
    object = names(x = ranking),
    expected = c(
      "ranking", "contributions", "weights", "directions",
      "n_candidates"
    )
  )
  expect_identical(
    object = names(x = ranking$ranking),
    expected = c("construct_id", "score", "rank")
  )
  expect_identical(
    object = names(x = ranking$contributions),
    expected = c(
      "construct_id", "region", "feature", "value", "normalised",
      "direction", "weight", "contribution"
    )
  )
  expect_identical(object = ranking$n_candidates, expected = 6L)
  expect_identical(object = nrow(x = ranking$ranking), expected = 6L)
  expect_identical(
    object = nrow(x = ranking$contributions), expected = 6L * 12L
  )
})

test_that("no coding sequence feature can reach a score", {
  ranking <- rankCandidates(feature_table = all_features)

  expect_identical(
    object = sort(x = unique(x = ranking$contributions$region)),
    expected = c("utr3", "utr5")
  )
  expect_false(object = any(ranking$contributions$region == "cds"))

  # CDS-only features are absent by name as well as by region.
  for (feature in c("gc3_content", "n_premature_stop",
                    "n_out_of_frame_atg")) {
    expect_false(
      object = feature %in% ranking$contributions$feature,
      info = feature
    )
  }

  # Removing the CDS rows entirely changes nothing, which is the
  # proof that they never contributed.
  without_cds <- all_features[all_features$region != "cds", ]
  expect_identical(
    object = rankCandidates(feature_table = without_cds)$ranking,
    expected = ranking$ranking
  )
})

test_that("positional rows do not reach a score either", {
  # Ranking uses the per-region summaries; the counts among them
  # already summarise the positional rows.
  summaries_only <- all_features[is.na(all_features$position), ]

  expect_identical(
    object = rankCandidates(feature_table = summaries_only)$ranking,
    expected = rankCandidates(feature_table = all_features)$ranking
  )
})

test_that("scores are a weighted mean of the contributions", {
  ranking <- rankCandidates(feature_table = all_features)

  # Each candidate's contributions sum exactly to its score.
  for (construct_id in ranking$ranking$construct_id) {
    rows <- ranking$contributions[
      ranking$contributions$construct_id == construct_id,
    ]
    expect_equal(
      object = sum(rows$contribution),
      expected = ranking$ranking$score[
        ranking$ranking$construct_id == construct_id
      ],
      info = construct_id
    )
    expect_equal(
      object = rows$contribution,
      expected = rows$normalised * rows$weight,
      info = construct_id
    )
  }

  expect_equal(object = sum(ranking$weights), expected = 1)
  expect_true(
    object = all(
      ranking$ranking$score >= 0 & ranking$ranking$score <= 1
    )
  )
  expect_true(object = all(ranking$contributions$normalised >= 0))
  expect_true(object = all(ranking$contributions$normalised <= 1))
})

test_that("normalisation orients every feature so higher is better", {
  ranking <- rankCandidates(feature_table = all_features)
  contributions <- ranking$contributions

  for (i in seq_len(length.out = nrow(x = expected_directions))) {
    feature <- expected_directions$feature[[i]]
    region <- expected_directions$region[[i]]
    rows <- contributions[
      contributions$feature == feature & contributions$region == region,
    ]
    if (length(x = unique(x = rows$value)) < 2L) {
      next
    }

    best_value <- rows$value[which.max(x = rows$normalised)]
    worst_value <- rows$value[which.min(x = rows$normalised)]

    if (expected_directions$direction[[i]] == 1L) {
      expect_gt(
        object = best_value, expected = worst_value,
        label = paste(region, feature)
      )
    } else {
      expect_lt(
        object = best_value, expected = worst_value,
        label = paste(region, feature)
      )
    }
  }
})

test_that("ranks follow the scores, best first", {
  ranking <- rankCandidates(feature_table = all_features)$ranking

  expect_identical(
    object = ranking$rank, expected = sort(x = ranking$rank)
  )
  expect_identical(
    object = ranking$score, expected = sort(x = ranking$score,
                                            decreasing = TRUE)
  )
  expect_identical(object = ranking$rank[[1]], expected = 1L)
  expect_identical(
    object = ranking$rank,
    expected = rank(x = -ranking$score, ties.method = "min")
  )
})

test_that("the ranking separates the planted good from the bad", {
  ranking <- rankCandidates(feature_table = all_features)$ranking
  rank_of <- function(construct_id) {
    ranking$rank[ranking$construct_id == construct_id]
  }

  # CONSTRUCT_03 carries three upstream AUGs and two uORFs;
  # CONSTRUCT_06 carries none and a strong Kozak context.
  expect_lt(
    object = rank_of("CONSTRUCT_06"),
    expected = rank_of("CONSTRUCT_03")
  )
  # CONSTRUCT_12 carries nine AU-rich pentamers.
  expect_lt(
    object = rank_of("CONSTRUCT_04"),
    expected = rank_of("CONSTRUCT_12")
  )

  # Ranking on upstream ORFs alone must put the constructs carrying
  # them last.
  uorf_only <- rankCandidates(
    feature_table = all_features,
    weights = c(n_upstream_orf = 1)
  )$ranking
  with_uorfs <- simulatedConstructs$construct_id[
    simulatedConstructs$n_upstream_orf > 0
  ]
  expect_true(
    object = all(
      uorf_only$rank[uorf_only$construct_id %in% with_uorfs] >
        max(uorf_only$rank[!uorf_only$construct_id %in% with_uorfs])
    )
  )
})

test_that("equal weighting is the default and is recorded", {
  ranking <- rankCandidates(feature_table = all_features)

  expected_keys <- paste(
    expected_directions$region, expected_directions$feature, sep = "."
  )

  expect_length(object = ranking$weights, n = 12L)
  expect_true(
    object = all(
      abs(ranking$weights - 1 / 12) < sqrt(x = .Machine$double.eps)
    )
  )
  # Weight names are region-qualified, because gc_content is scored
  # in both UTRs and a bare name would not identify which.
  expect_identical(
    object = sort(x = names(x = ranking$weights)),
    expected = sort(x = expected_keys)
  )
  expect_identical(
    object = anyDuplicated(x = names(x = ranking$weights)),
    expected = 0L
  )
  expect_identical(
    object = ranking$directions[expected_keys],
    expected = stats::setNames(
      object = expected_directions$direction, nm = expected_keys
    )
  )
})

test_that("the weights recorded are the weights applied", {
  # These must agree: a ranking is only traceable if the weights it
  # reports are the ones that produced it.
  for (weights in list(
    NULL,
    c(gc_content = 1),
    c(utr3.gc_content = 1),
    c(kozak_score_scaled = 2, n_upstream_orf = 1)
  )) {
    ranking <- rankCandidates(
      feature_table = six_features, weights = weights
    )
    applied <- unique(
      ranking$contributions[, c("region", "feature", "weight")]
    )
    applied$key <- paste(
      applied$region, applied$feature, sep = "."
    )

    expect_equal(
      object = applied$weight[order(applied$key)],
      expected = as.numeric(
        x = ranking$weights[sort(x = names(x = ranking$weights))]
      )
    )
  }
})

test_that("a bare name covers every region a feature is scored in", {
  both <- rankCandidates(
    feature_table = six_features, weights = c(gc_content = 1)
  )
  expect_identical(
    object = names(x = both$weights[both$weights > 0]),
    expected = c("utr5.gc_content", "utr3.gc_content")
  )
  expect_equal(
    object = unname(obj = both$weights[both$weights > 0]),
    expected = c(0.5, 0.5)
  )

  # Naming the region picks out one of them.
  one <- rankCandidates(
    feature_table = six_features, weights = c(utr3.gc_content = 1)
  )
  expect_identical(
    object = names(x = one$weights[one$weights > 0]),
    expected = "utr3.gc_content"
  )

  # The two give different rankings, so the distinction is real.
  expect_false(
    object = identical(
      x = both$ranking$score, y = one$ranking$score
    )
  )

  # Setting the same feature twice is refused rather than silently
  # letting one win.
  expect_error(
    object = rankCandidates(
      feature_table = six_features,
      weights = c(gc_content = 1, utr3.gc_content = 2)
    ),
    regexp = "more than once"
  )
})

test_that("a supplied weight vector gives unnamed features none", {
  ranking <- rankCandidates(
    feature_table = all_features,
    weights = c(kozak_score_scaled = 2, n_upstream_orf = 1)
  )

  expect_identical(object = sum(ranking$weights > 0), expected = 2L)
  expect_equal(
    object = unname(obj = ranking$weights[["utr5.kozak_score_scaled"]]),
    expected = 2 / 3
  )
  expect_equal(
    object = unname(obj = ranking$weights[["utr5.n_upstream_orf"]]),
    expected = 1 / 3
  )
  expect_equal(object = sum(ranking$weights), expected = 1)

  # A feature with no weight contributes nothing to any score.
  unweighted <- ranking$contributions[
    ranking$contributions$weight == 0,
  ]
  expect_true(object = all(unweighted$contribution == 0))
})

test_that("weights change the ordering and are recorded with it", {
  by_kozak <- rankCandidates(
    feature_table = all_features, weights = c(kozak_score_scaled = 1)
  )
  by_are <- rankCandidates(
    feature_table = all_features, weights = c(n_are_pentamer = 1)
  )

  expect_false(
    object = identical(
      x = by_kozak$ranking$construct_id,
      y = by_are$ranking$construct_id
    )
  )
  # The weights returned are the ones that produced the ordering,
  # which is what makes a ranking reproducible.
  expect_identical(
    object = names(x = by_kozak$weights[by_kozak$weights > 0]),
    expected = "utr5.kozak_score_scaled"
  )
  expect_identical(
    object = rankCandidates(
      feature_table = all_features, weights = by_kozak$weights
    )$ranking,
    expected = by_kozak$ranking
  )
})

test_that("a feature with no spread is neutral for everyone", {
  # Two copies of one construct differ in nothing, so no feature can
  # tell them apart.
  twins <- simulatedConstructs[c(4, 4), ]
  twins$construct_id <- c("twin_a", "twin_b")
  ranking <- rankCandidates(
    feature_table = profileConstruct(constructs = twins)
  )

  expect_true(object = all(ranking$contributions$normalised == 0.5))
  expect_equal(object = ranking$ranking$score, expected = c(0.5, 0.5))
  expect_identical(object = ranking$ranking$rank, expected = c(1L, 1L))
})

test_that("at least two candidates are needed to compare", {
  one_candidate <- all_features[
    all_features$construct_id == "CONSTRUCT_01",
  ]
  expect_error(
    object = rankCandidates(feature_table = one_candidate),
    regexp = "at least two candidates"
  )
})

test_that("missing default features warn rather than pass unnoticed", {
  # Replacing the built-in motifs removes n_motif_G_quadruplex, which
  # the default weighting expects.
  replaced <- profileConstruct(
    constructs = simulatedConstructs[1:3, ], motifs = c(mine = "ATTTA")
  )

  expect_warning(
    object = rankCandidates(feature_table = replaced),
    regexp = "absent from the feature table"
  )
  ranking <- suppressWarnings(
    rankCandidates(feature_table = replaced)
  )
  expect_false(
    object = "n_motif_G_quadruplex" %in% names(x = ranking$weights)
  )
  expect_equal(object = sum(ranking$weights), expected = 1)
})

test_that("invalid weights are rejected with an informative message", {
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = c(nonsense = 1)
    ),
    regexp = "cannot be weighted"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = c(gc3_content = 1)
    ),
    regexp = "not scored or not present"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = c(2, 1)
    ),
    regexp = "must be named"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features,
      weights = c(gc_content = 1, gc_content = 2)
    ),
    regexp = "unique names"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = c(gc_content = NA_real_)
    ),
    regexp = "finite numbers"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = c(gc_content = 0)
    ),
    regexp = "at least one positive weight"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features, weights = "equal"
    ),
    regexp = "must be a named numeric vector"
  )
})

test_that("negative weights are refused, since the score is a mean", {
  # A negative weight can make the weights sum to zero, which would
  # give an infinite score, and otherwise pushes scores outside the
  # documented 0 to 1 range and can invert the ordering.
  expect_error(
    object = rankCandidates(
      feature_table = all_features,
      weights = c(gc_content = 1, n_are_pentamer = -1)
    ),
    regexp = "must not be negative"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features,
      weights = c(gc_content = 1, n_are_pentamer = -1)
    ),
    regexp = "weighted mean"
  )
})

test_that("an unusable feature table is rejected", {
  expect_error(
    object = rankCandidates(feature_table = NULL), regexp = "not NULL"
  )
  expect_error(
    object = rankCandidates(feature_table = data.frame(a = 1)),
    regexp = "missing"
  )
  expect_error(
    object = rankCandidates(feature_table = "features"),
    regexp = "must be a data frame"
  )
  expect_error(
    object = rankCandidates(
      feature_table = all_features[0, ]
    ),
    regexp = "no rows"
  )
})

test_that("printing a ranking states what it is and is not", {
  ranking <- rankCandidates(feature_table = six_features)
  printed <- paste(
    capture.output(expr = print(x = ranking)), collapse = " "
  )

  expect_match(object = printed, regexp = "6 candidate constructs")
  expect_match(object = printed, regexp = "Not a prediction")
  expect_match(object = printed, regexp = "Weights used")
  expect_invisible(call = print(x = ranking))
})
