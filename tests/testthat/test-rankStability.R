data(simulatedConstructs, package = "UTRunway")
stability_features <- profileConstruct(
  constructs = simulatedConstructs[1:6, ]
)

test_that("the returned object has the documented shape", {
  set.seed(seed = 1)
  stability <- rankStability(
    feature_table = stability_features, n = 40
  )

  expect_s3_class(object = stability, class = "utrunway_rank_stability")
  expect_identical(
    object = names(x = stability),
    expected = c(
      "distribution", "summary", "weights", "n_replicates",
      "concentration"
    )
  )
  expect_identical(
    object = names(x = stability$distribution),
    expected = c("construct_id", "rank", "n_times", "proportion")
  )
  expect_identical(
    object = names(x = stability$summary),
    expected = c(
      "construct_id", "base_rank", "median_rank", "min_rank",
      "max_rank", "proportion_best"
    )
  )
  expect_identical(object = stability$n_replicates, expected = 40L)
  expect_identical(object = stability$concentration, expected = 50)
  expect_identical(object = nrow(x = stability$summary), expected = 6L)
  # Every candidate against every rank, including combinations never
  # seen, so the result can be drawn as a complete grid.
  expect_identical(
    object = nrow(x = stability$distribution), expected = 6L * 6L
  )
})

test_that("each replicate produces a valid ranking", {
  set.seed(seed = 2)
  n_replicates <- 50L
  stability <- rankStability(
    feature_table = stability_features, n = n_replicates
  )
  n_candidates <- nrow(x = stability$summary)

  # Each candidate holds exactly one rank per replicate.
  per_candidate <- tapply(
    X = stability$distribution$n_times,
    INDEX = stability$distribution$construct_id,
    FUN = sum
  )
  expect_true(object = all(as.numeric(x = per_candidate) == n_replicates))

  # And each rank is taken by exactly one candidate per replicate, so
  # every replicate is a permutation rather than an arbitrary set of
  # numbers.
  per_rank <- tapply(
    X = stability$distribution$n_times,
    INDEX = stability$distribution$rank,
    FUN = sum
  )
  expect_true(object = all(as.numeric(x = per_rank) == n_replicates))

  expect_identical(
    object = sort(x = unique(x = stability$distribution$rank)),
    expected = seq_len(length.out = n_candidates)
  )
  proportions <- tapply(
    X = stability$distribution$proportion,
    INDEX = stability$distribution$construct_id,
    FUN = sum
  )
  expect_equal(
    object = as.numeric(x = proportions),
    expected = rep(x = 1, times = n_candidates)
  )
})

test_that("the summary agrees with the distribution it came from", {
  set.seed(seed = 3)
  stability <- rankStability(
    feature_table = stability_features, n = 60
  )

  for (construct_id in stability$summary$construct_id) {
    summary_row <- stability$summary[
      stability$summary$construct_id == construct_id,
    ]
    seen <- stability$distribution[
      stability$distribution$construct_id == construct_id &
        stability$distribution$n_times > 0,
    ]

    expect_identical(
      object = summary_row$min_rank,
      expected = min(seen$rank),
      info = construct_id
    )
    expect_identical(
      object = summary_row$max_rank,
      expected = max(seen$rank),
      info = construct_id
    )
    expect_equal(
      object = summary_row$proportion_best,
      expected = stability$distribution$proportion[
        stability$distribution$construct_id == construct_id &
          stability$distribution$rank == 1L
      ],
      info = construct_id
    )
    expect_true(
      object = summary_row$median_rank >= summary_row$min_rank &&
        summary_row$median_rank <= summary_row$max_rank,
      info = construct_id
    )
  }

  # The summary is ordered by median rank, best first.
  expect_identical(
    object = stability$summary$median_rank,
    expected = sort(x = stability$summary$median_rank)
  )
})

test_that("base_rank is the ranking under the unperturbed weights", {
  set.seed(seed = 4)
  stability <- rankStability(
    feature_table = stability_features, n = 30
  )
  base <- rankCandidates(feature_table = stability_features)$ranking

  expect_identical(
    object = stability$summary$base_rank[
      order(stability$summary$construct_id)
    ],
    expected = base$rank[order(base$construct_id)]
  )
  # With no weights supplied, the draws are centred on the default
  # equal weighting, and that is what gets recorded.
  expect_identical(
    object = stability$weights,
    expected = rankCandidates(
      feature_table = stability_features
    )$weights
  )
  expect_equal(object = sum(stability$weights), expected = 1)
})

test_that("the weights centred on are the ones recorded", {
  set.seed(seed = 5)
  supplied <- c(kozak_score_scaled = 3, n_upstream_orf = 1)
  stability <- rankStability(
    feature_table = stability_features, n = 30, weights = supplied
  )

  expect_equal(object = sum(stability$weights), expected = 1)
  expect_equal(
    object = unname(
      obj = stability$weights[["utr5.kozak_score_scaled"]]
    ),
    expected = 0.75
  )
  expect_identical(
    object = sum(stability$weights > 0), expected = 2L
  )
  # A feature excluded by the weights stays excluded in every
  # replicate, since dropping it is a decision, not a quantity to
  # perturb.
  expect_identical(
    object = stability$summary$base_rank,
    expected = rankCandidates(
      feature_table = stability_features, weights = supplied
    )$ranking$rank[
      match(
        x = stability$summary$construct_id,
        table = rankCandidates(
          feature_table = stability_features, weights = supplied
        )$ranking$construct_id
      )
    ]
  )
})

test_that("results are reproducible under a set seed", {
  set.seed(seed = 11)
  first <- rankStability(feature_table = stability_features, n = 25)
  set.seed(seed = 11)
  second <- rankStability(feature_table = stability_features, n = 25)

  expect_identical(object = first, expected = second)

  # A different seed generally gives a different sample, so the
  # result really does depend on the draw.
  set.seed(seed = 12)
  third <- rankStability(feature_table = stability_features, n = 25)
  expect_false(
    object = identical(x = first$distribution, y = third$distribution)
  )
})

test_that("concentration controls how far the ranks move", {
  set.seed(seed = 21)
  tight <- rankStability(
    feature_table = stability_features, n = 100, concentration = 500
  )
  set.seed(seed = 21)
  loose <- rankStability(
    feature_table = stability_features, n = 100, concentration = 1
  )

  tight_spread <- mean(tight$summary$max_rank - tight$summary$min_rank)
  loose_spread <- mean(loose$summary$max_rank - loose$summary$min_rank)
  expect_lt(object = tight_spread, expected = loose_spread)

  # Weights barely perturbed reproduce the unperturbed ranking.
  expect_identical(
    object = tight$summary$median_rank,
    expected = as.numeric(x = tight$summary$base_rank)
  )
})

test_that("a clearly better candidate holds its rank", {
  # Over the whole set, the construct with no upstream AUGs, no
  # AU-rich elements, a strong Kozak context and a canonical polyA
  # signal should lead under almost any weighting.
  all_features <- profileConstruct(constructs = simulatedConstructs)
  set.seed(seed = 31)
  stability <- rankStability(feature_table = all_features, n = 100)

  leader <- rankCandidates(
    feature_table = all_features
  )$ranking$construct_id[[1]]
  leader_row <- stability$summary[
    stability$summary$construct_id == leader,
  ]

  expect_identical(object = leader_row$base_rank, expected = 1L)
  expect_gt(object = leader_row$proportion_best, expected = 0.5)
})

test_that("invalid arguments are rejected with a clear message", {
  expect_error(
    object = rankStability(feature_table = stability_features, n = 1),
    regexp = "`n` must be at least 2"
  )
  expect_error(
    object = rankStability(
      feature_table = stability_features, n = 2.5
    ),
    regexp = "`n` must be a whole number"
  )
  expect_error(
    object = rankStability(
      feature_table = stability_features, concentration = 0
    ),
    regexp = "`concentration` must be a single positive number"
  )
  expect_error(
    object = rankStability(
      feature_table = stability_features, concentration = -5
    ),
    regexp = "single positive number"
  )
  expect_error(
    object = rankStability(
      feature_table = stability_features, concentration = c(1, 2)
    ),
    regexp = "single positive number"
  )
  expect_error(
    object = rankStability(feature_table = NULL), regexp = "not NULL"
  )
  expect_error(
    object = rankStability(
      feature_table = stability_features, weights = c(nonsense = 1)
    ),
    regexp = "cannot be weighted"
  )
})

test_that("printing a stability summary states what it is and is not", {
  set.seed(seed = 41)
  stability <- rankStability(
    feature_table = stability_features, n = 20
  )
  printed <- paste(
    capture.output(expr = print(x = stability)), collapse = " "
  )

  expect_match(object = printed, regexp = "20 perturbed weightings")
  expect_match(object = printed, regexp = "concentration")
  expect_match(object = printed, regexp = "Not a prediction")
  expect_invisible(call = print(x = stability))
})
