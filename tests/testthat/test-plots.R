# Plot tests that can fail.
#
# Checking that a plotting function returns a ggplot object is close
# to true by construction: a function that builds the wrong figure
# entirely still returns one. So these tests build each plot and read
# the values that were actually drawn, out of the built layers. The
# score breakdown's reference rule once sat at 0.083 instead of 0.583
# while every other test passed; the assertions below are the kind
# that would have caught it.

data(simulatedConstructs, package = "UTRunway")
plot_features <- profileConstruct(
  constructs = simulatedConstructs[c(1, 3, 5, 6, 10, 12), ]
)
plot_ranking <- rankCandidates(feature_table = plot_features)

# The data actually drawn by a built plot, as one data frame per
# layer.
builtLayers <- function(plot) {
  ggplot2::ggplot_build(plot)$data
}

# The first built layer carrying a given column.
layerWith <- function(plot, column) {
  layers <- builtLayers(plot = plot)
  found <- Filter(
    f = function(layer) column %in% names(x = layer), x = layers
  )
  if (length(x = found) == 0L) {
    return(NULL)
  }
  found[[1]]
}

# Axis labels a built plot will actually print.
axisLabels <- function(plot, axis) {
  built <- ggplot2::ggplot_build(plot)
  built$layout$panel_params[[1]][[axis]]$get_labels()
}

test_that("every plot builds without error on the bundled data", {
  set.seed(seed = 1)
  stability <- rankStability(feature_table = plot_features, n = 20)

  plots <- list(
    map = plotConstructMap(
      feature_table = plot_features, construct_id = "CONSTRUCT_03"
    ),
    profile = plotRegionProfile(
      feature_table = plot_features, construct_id = "CONSTRUCT_05"
    ),
    comparison = plotCandidateComparison(
      feature_table = plot_features
    ),
    breakdown = plotScoreBreakdown(ranking = plot_ranking),
    ranking = plotCandidateRanking(stability = stability)
  )

  for (name in names(x = plots)) {
    # Building is what catches an error that construction does not:
    # a scale with no matching data, an aesthetic of the wrong
    # length, a factor level that is never used.
    expect_no_error(
      object = ggplot2::ggplot_build(plots[[name]]),
      message = name
    )
    expect_s3_class(object = plots[[name]], class = "ggplot")
  }
  expect_s3_class(object = plots$map, class = "patchwork")
})

test_that("map positions are region offset plus region position", {
  # The map lays the regions end to end, so a 3' UTR feature at
  # position 40 must be drawn at 40 plus the lengths before it. An
  # off-by-one or a missing offset would put features in the wrong
  # region entirely.
  rows <- plot_features[
    plot_features$construct_id == "CONSTRUCT_12", ,
    drop = FALSE
  ]
  lengths <- vapply(
    X = c("utr5", "cds", "utr3"),
    FUN = function(region) {
      rows$value[rows$region == region & rows$feature == "length"]
    },
    FUN.VALUE = numeric(length = 1L)
  )
  utr3_start <- lengths[["utr5"]] + lengths[["cds"]]

  map <- plotConstructMap(
    feature_table = plot_features, construct_id = "CONSTRUCT_12"
  )
  drawn <- layerWith(plot = map[[2]], column = "shape")

  pentamers <- rows[rows$feature == "are_pentamer", , drop = FALSE]
  expect_gt(object = nrow(x = pentamers), expected = 0)
  expected_x <- utr3_start + pentamers$position + (5 - 1) / 2

  expect_true(
    object = all(round(x = sort(x = expected_x), digits = 6) %in%
                   round(x = sort(x = drawn$x), digits = 6))
  )
  # Every mark lies inside the construct.
  expect_true(object = all(drawn$x >= 0))
  expect_true(object = all(drawn$x <= sum(lengths)))
})

test_that("the map omits CDS rows unless they are asked for", {
  without <- plotConstructMap(
    feature_table = plot_features, construct_id = "CONSTRUCT_03"
  )
  with_context <- plotConstructMap(
    feature_table = plot_features, construct_id = "CONSTRUCT_03",
    show_cds_context = TRUE
  )

  rows_without <- axisLabels(plot = without[[2]], axis = "y")
  rows_with <- axisLabels(plot = with_context[[2]], axis = "y")

  expect_false(
    object = any(grepl(pattern = "CDS context", x = rows_without))
  )
  expect_true(
    object = any(grepl(pattern = "CDS context", x = rows_with))
  )
  # Asking for context adds rows and removes none.
  expect_true(object = all(rows_without %in% rows_with))
  expect_gt(
    object = length(x = rows_with), expected = length(x = rows_without)
  )
})

test_that("coding sequence marks are drawn as context, not as ink", {
  # Two of the context markers use plotting symbols that take no
  # fill, so their colour has to come from the outline. Setting only
  # the fill would leave them drawn in ink, as prominent as a scored
  # feature.
  palette <- utrunwayPalette()
  map <- plotConstructMap(
    feature_table = plot_features, construct_id = "CONSTRUCT_03",
    show_cds_context = TRUE
  )
  drawn <- layerWith(plot = map[[2]], column = "shape")
  unfillable <- drawn[drawn$shape %in% c(3L, 4L), , drop = FALSE]

  expect_gt(object = nrow(x = unfillable), expected = 0)
  expect_true(
    object = all(unfillable$colour == palette[["context"]])
  )
  expect_false(object = any(unfillable$colour == palette[["ink"]]))
})

test_that("every map row is named, so no mark is colour-only", {
  map <- plotConstructMap(
    feature_table = plot_features, construct_id = "CONSTRUCT_03",
    show_cds_context = TRUE
  )
  labels <- axisLabels(plot = map[[2]], axis = "y")

  expect_gt(object = length(x = labels), expected = 0)
  expect_false(object = any(is.na(x = labels)))
  expect_true(object = all(nzchar(x = labels)))
  # An upstream AUG's initiation context belongs in the row name
  # rather than in a shade of one colour.
  expect_true(
    object = any(grepl(pattern = "^uAUG, ", x = labels))
  )
})

test_that("the region profile draws the GC values it was given", {
  profile <- plotRegionProfile(
    feature_table = plot_features, construct_id = "CONSTRUCT_05"
  )
  rows <- plot_features[
    plot_features$construct_id == "CONSTRUCT_05" &
      plot_features$feature == "gc_window",
    ,
    drop = FALSE
  ]
  # The line layer carries both panels' series; flipped_aes marks it
  # out from the shading, the boundary rules and the GC3 annotation,
  # all of which also have a group column.
  lines <- layerWith(plot = profile, column = "flipped_aes")
  gc_panel <- lines[lines$PANEL == levels(x = lines$PANEL)[[1]], ,
                    drop = FALSE]

  expect_identical(
    object = nrow(x = gc_panel), expected = nrow(x = rows)
  )
  expect_equal(
    object = sort(x = round(x = gc_panel$y, digits = 8)),
    expected = sort(x = round(x = rows$value, digits = 8))
  )
  # Both panels keep their values inside the ranges they can take.
  expect_true(object = all(gc_panel$y >= 0 & gc_panel$y <= 1))

  stem_panel <- lines[lines$PANEL == levels(x = lines$PANEL)[[2]], ,
                      drop = FALSE]
  expect_true(object = all(stem_panel$y >= 0))
})

test_that("comparison scales features to span zero to one", {
  comparison <- plotCandidateComparison(
    feature_table = plot_features
  )
  drawn <- layerWith(plot = comparison, column = "y")
  drawn <- drawn[is.finite(x = drawn$y), , drop = FALSE]

  expect_gt(object = nrow(x = drawn), expected = 0)
  expect_true(object = all(drawn$y >= 0))
  expect_true(object = all(drawn$y <= 1))
  # At least one feature must actually reach both ends, or the
  # scaling is not doing anything.
  expect_true(object = any(drawn$y == 0))
  expect_true(object = any(drawn$y == 1))
})

test_that("highlighted candidates are drawn and others are not", {
  # Colours have to be gathered from every layer: the first one
  # carrying a colour column is the background, not the highlight.
  allColours <- function(plot) {
    unique(x = unlist(x = lapply(
      X = builtLayers(plot = plot),
      FUN = function(layer) {
        if ("colour" %in% names(x = layer)) layer$colour else NULL
      }
    )))
  }

  highlighted <- plotCandidateComparison(
    feature_table = plot_features,
    highlight = c("CONSTRUCT_03", "CONSTRUCT_06")
  )
  used <- intersect(
    x = allColours(plot = highlighted),
    y = candidateHighlightColours()
  )
  expect_length(object = used, n = 2L)

  plain <- plotCandidateComparison(feature_table = plot_features)
  expect_length(
    object = intersect(
      x = allColours(plot = plain), y = candidateHighlightColours()
    ),
    n = 0L
  )
})

test_that("the breakdown's ceiling is the region's total weight", {
  # The regression this file exists for. The rule was once computed
  # as sum(unique(weights)), which with equal weighting collapsed
  # twelve identical values into one and drew the rule at 0.083.
  breakdown <- plotScoreBreakdown(ranking = plot_ranking)
  weights <- plot_ranking$weights

  expected <- c(
    utr5 = sum(weights[grepl(pattern = "^utr5\\.", x = names(weights))]),
    utr3 = sum(weights[grepl(pattern = "^utr3\\.", x = names(weights))])
  )
  rule <- layerWith(plot = breakdown, column = "xintercept")

  expect_false(object = is.null(x = rule))
  expect_equal(
    object = sort(x = round(x = rule$xintercept, digits = 6)),
    expected = sort(x = round(x = unname(obj = expected), digits = 6))
  )
  # A single weight is not the ceiling, which is what the bug drew.
  expect_false(
    object = any(abs(rule$xintercept - min(weights)) < 1e-9)
  )
})

test_that("breakdown bar lengths add up to the composite scores", {
  breakdown <- plotScoreBreakdown(ranking = plot_ranking)
  bars <- layerWith(plot = breakdown, column = "fill")
  bars <- bars[!is.na(bars$xmin) & !is.na(bars$xmax), , drop = FALSE]

  drawn_total <- sum(bars$xmax - bars$xmin)
  expect_equal(
    object = drawn_total,
    expected = sum(plot_ranking$ranking$score),
    tolerance = 1e-6
  )
})

test_that("the ranking plot draws the stability figures given", {
  set.seed(seed = 2)
  stability <- rankStability(feature_table = plot_features, n = 40)
  plot <- plotCandidateRanking(stability = stability)

  segments <- layerWith(plot = plot, column = "xend")
  expect_equal(
    object = sort(x = segments$x),
    expected = sort(x = as.numeric(x = stability$summary$min_rank))
  )
  expect_equal(
    object = sort(x = segments$xend),
    expected = sort(x = as.numeric(x = stability$summary$max_rank))
  )

  # A candidate that never moved must have no span to draw.
  settled <- stability$summary[
    stability$summary$min_rank == stability$summary$max_rank, ,
    drop = FALSE
  ]
  if (nrow(x = settled) > 0L) {
    zero_length <- segments[segments$x == segments$xend, , drop = FALSE]
    expect_identical(
      object = nrow(x = zero_length), expected = nrow(x = settled)
    )
  }
})

test_that("only candidates that led are labelled as having led", {
  set.seed(seed = 3)
  stability <- rankStability(feature_table = plot_features, n = 40)
  plot <- plotCandidateRanking(stability = stability)
  labels <- layerWith(plot = plot, column = "label")

  n_led <- sum(stability$summary$proportion_best > 0)
  expect_identical(object = nrow(x = labels), expected = n_led)
  expect_false(
    object = any(grepl(pattern = "^0% first$", x = labels$label))
  )
})

test_that("plots reject input they cannot draw", {
  expect_error(
    object = plotConstructMap(
      feature_table = plot_features, construct_id = "NO_SUCH"
    ),
    regexp = "is not in `feature_table`"
  )
  expect_error(
    object = plotConstructMap(
      feature_table = plot_features,
      construct_id = c("CONSTRUCT_01", "CONSTRUCT_03")
    ),
    regexp = "single character string"
  )
  expect_error(
    object = plotConstructMap(
      feature_table = plot_features, construct_id = "CONSTRUCT_01",
      show_cds_context = "yes"
    ),
    regexp = "must be TRUE or FALSE"
  )
  expect_error(
    object = plotRegionProfile(
      feature_table = plot_features, construct_id = NULL
    ),
    regexp = "must name one construct"
  )
  expect_error(
    object = plotCandidateComparison(
      feature_table = plot_features,
      highlight = c("CONSTRUCT_01", "CONSTRUCT_03", "CONSTRUCT_05",
                    "CONSTRUCT_06")
    ),
    regexp = "at most three candidates"
  )
  expect_error(
    object = plotCandidateComparison(
      feature_table = plot_features, highlight = "NO_SUCH"
    ),
    regexp = "not in `feature_table`"
  )
  expect_error(
    object = plotScoreBreakdown(ranking = plot_features),
    regexp = "must be an object from rankCandidates"
  )
  expect_error(
    object = plotCandidateRanking(stability = plot_ranking),
    regexp = "must be an object from rankStability"
  )
  expect_error(
    object = plotConstructMap(
      feature_table = data.frame(a = 1), construct_id = "x"
    ),
    regexp = "missing"
  )
})

test_that("a construct with nothing to annotate still draws", {
  # A design with no detected feature must give a figure that says
  # so, rather than failing or drawing an empty panel with no
  # explanation.
  # An AGC repeat carries nothing the package detects: it has no
  # ATG, no four G-triplets, no AUUUA, no AAUAAA and no seed site,
  # and it does not open with a pyrimidine run, which would be a TOP
  # motif.
  plain <- data.frame(
    construct_id = "plain",
    utr5 = paste(rep(x = "AGC", times = 30L), collapse = ""),
    cds = paste0(
      "ATG", paste(rep(x = "GCA", times = 40L), collapse = ""), "TAA"
    ),
    utr3 = paste(rep(x = "AGC", times = 70L), collapse = ""),
    stringsAsFactors = FALSE
  )
  features <- profileConstruct(constructs = plain)

  map <- plotConstructMap(
    feature_table = features, construct_id = "plain"
  )
  expect_no_error(object = ggplot2::ggplot_build(map[[2]]))
  expect_match(
    object = map[[2]]$labels$title, regexp = "none found"
  )
})

test_that("the palette gives every class a shape and a label", {
  palette <- utrunwayPalette()
  shapes <- utrunwayShapes()
  labels <- utrunwayLabels()

  feature_keys <- c(
    "upstream_aug_weak", "upstream_aug_adequate",
    "upstream_aug_strong", "top_motif", "g_quadruplex",
    "polya_signal", "are_pentamer", "mirna_seed"
  )
  expect_true(object = all(feature_keys %in% names(x = palette)))
  expect_false(object = anyNA(x = palette[feature_keys]))

  shape_keys <- c(
    "upstream_aug", "top_motif", "g_quadruplex",
    "polya_signal_canonical", "polya_signal_variant",
    "are_pentamer", "mirna_seed"
  )
  expect_true(object = all(shape_keys %in% names(x = shapes)))
  # Only symbols that accept a fill, or a mark loses its colour.
  expect_true(object = all(shapes[shape_keys] >= 21L))
  expect_true(object = all(shapes[shape_keys] <= 25L))

  label_keys <- c(
    "upstream_aug_weak", "upstream_aug_adequate",
    "upstream_aug_strong", "upstream_orf", "top_motif",
    "g_quadruplex", "polya_signal_canonical",
    "polya_signal_variant", "are_pentamer", "mirna_seed"
  )
  expect_true(object = all(label_keys %in% names(x = labels)))
  expect_true(object = all(nzchar(x = labels[label_keys])))
})

test_that("no two rows look alike in both symbol and colour", {
  # Every row of the annotation track is drawn together, so there
  # are more classes than fillable symbols and a symbol has to be
  # reused. What must hold is that no two classes share a symbol and
  # a colour: a square is either the aqua TOP motif or the green
  # polyA variant, never ambiguous.
  spec <- annotationRowSpec()
  appearance <- paste(
    annotationShapes()[spec$shape_key],
    utrunwayPalette()[spec$colour_key]
  )

  expect_identical(
    object = anyDuplicated(x = appearance), expected = 0L
  )
  # A symbol genuinely is reused, so the check above is doing work
  # rather than passing because everything happens to be distinct.
  expect_gt(
    object = length(x = annotationShapes()[spec$shape_key]),
    expected = length(
      x = unique(x = annotationShapes()[spec$shape_key])
    )
  )
})
