#' Scoring direction of every ranked feature
#'
#' @description
#' The table that decides which features are scored and which way
#' round. A reader can see every judgement the ranking rests on here,
#' in one place, and change it by editing this function.
#'
#' \code{direction} is \code{1} when a larger value is better for the
#' design's prospects and \code{-1} when a larger value is worse.
#' Features are keyed by region as well as name, so the same
#' measurement can be scored in one region and left alone in another.
#'
#' A feature absent from this table is never scored. That includes
#' every coding sequence feature, by design, and the features listed
#' in the \code{unscored} section of \code{\link{rankCandidates}}'s
#' documentation, each of which lacks a direction this package is
#' willing to assert.
#'
#' @return A data frame with columns \code{region}, \code{feature},
#'   \code{direction} and \code{rationale}.
#' @noRd
featureDirections <- function() {
  directions <- rbind(
    data.frame(
      region = "utr5", feature = "n_upstream_aug", direction = -1L,
      rationale = paste(
        "An upstream AUG captures scanning ribosomes and reduces",
        "initiation at the main start codon."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "n_upstream_orf", direction = -1L,
      rationale = paste(
        "A complete upstream ORF is the more repressive form of an",
        "upstream AUG."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "n_out_of_frame_overlap",
      direction = -1L,
      rationale = paste(
        "An upstream AUG reading through the main start codon out of",
        "frame produces no functional protein from that event."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "kozak_score_scaled", direction = 1L,
      rationale = paste(
        "A stronger initiation context raises the proportion of",
        "scanning ribosomes that initiate at the main start codon."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "gc_content", direction = -1L,
      rationale = paste(
        "A GC-rich 5' UTR tends to be more structured, which impedes",
        "scanning."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "max_stem_nt", direction = -1L,
      rationale = paste(
        "A longer complementary run is a longer possible stem, and",
        "stable structure near the cap impedes initiation."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr5", feature = "n_motif_G_quadruplex",
      direction = -1L,
      rationale = paste(
        "A quadruplex-prone run in the 5' UTR is generally",
        "inhibitory to scanning."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr3", feature = "gc_content", direction = 1L,
      rationale = paste(
        "An AU-rich 3' UTR is the destabilising case, so higher GC",
        "means less of the composition that recruits decay factors.",
        "Correlated with the AU-rich element counts below, which",
        "measure the same tendency through a specific motif."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr3", feature = "n_are_pentamer", direction = -1L,
      rationale = paste(
        "AUUUA pentamers recruit destabilising factors and shorten",
        "transcript half-life."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr3", feature = "n_are_clustered", direction = -1L,
      rationale = paste(
        "Clustered and overlapping pentamers drive decay far more",
        "strongly than isolated ones."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr3", feature = "n_polya_canonical", direction = 1L,
      rationale = paste(
        "A canonical AAUAAA hexamer near the 3' end supports",
        "efficient cleavage and polyadenylation."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      region = "utr3", feature = "n_mirna_seed", direction = -1L,
      rationale = paste(
        "Each seed match site is a place a miRNA could repress the",
        "transcript."
      ),
      stringsAsFactors = FALSE
    )
  )

  rownames(x = directions) <- NULL
  directions
}

#' Build the oriented, normalised score matrix
#'
#' Normalisation and orientation do not depend on the weights, so this
#' is computed once and reused for every weighting, which is what
#' makes \code{\link{rankStability}} affordable.
#'
#' @param feature_table A tidy feature table.
#' @return A list with the raw and normalised matrices, the features
#'   used and their directions.
#' @noRd
buildScoreMatrix <- function(feature_table) {
  directions <- featureDirections()

  # Scoring uses the untranslated regions only. This is the filter
  # that enforces it.
  scoreable <- feature_table[feature_table$region != "cds", , drop = FALSE]
  scoreable <- scoreable[is.na(scoreable$position), , drop = FALSE]

  wanted <- paste(directions$region, directions$feature)
  present <- paste(scoreable$region, scoreable$feature)
  scoreable <- scoreable[present %in% wanted, , drop = FALSE]

  construct_ids <- unique(x = feature_table$construct_id)
  if (length(x = construct_ids) < 2L) {
    stop(
      "Ranking is a comparison within the submitted set, so it needs ",
      "at least two candidates, but ", length(x = construct_ids),
      " was supplied.",
      call. = FALSE
    )
  }

  available <- directions[
    paste(directions$region, directions$feature) %in%
      unique(x = paste(scoreable$region, scoreable$feature)),
    ,
    drop = FALSE
  ]
  missing_features <- setdiff(
    x = directions$feature, y = available$feature
  )
  if (length(x = missing_features) > 0L) {
    warning(
      "These features are scored by default but are absent from the ",
      "feature table, so they take no part in the ranking: ",
      paste(missing_features, collapse = ", "),
      ". This happens when profileConstruct() was given a `motifs` ",
      "or `seeds` argument that replaced the defaults.",
      call. = FALSE
    )
  }
  if (nrow(x = available) == 0L) {
    stop(
      "None of the features that can be scored are present in ",
      "`feature_table`. Was it produced by profileConstruct()?",
      call. = FALSE
    )
  }

  # Features are keyed by region as well as name, because a name
  # alone is not unique: gc_content is scored in both untranslated
  # regions, in opposite directions. A duplicated name would make
  # weight lookup silently resolve to whichever came first.
  available$key <- paste(available$region, available$feature, sep = ".")

  raw <- matrix(
    data = NA_real_,
    nrow = length(x = construct_ids),
    ncol = nrow(x = available),
    dimnames = list(construct_ids, available$key)
  )

  for (j in seq_len(length.out = nrow(x = available))) {
    rows <- scoreable[
      scoreable$region == available$region[[j]] &
        scoreable$feature == available$feature[[j]],
      ,
      drop = FALSE
    ]
    if (anyDuplicated(x = rows$construct_id) != 0L) {
      stop(
        "`feature_table` has more than one ",
        available$feature[[j]], " value for a single construct.",
        call. = FALSE
      )
    }
    raw[, j] <- rows$value[match(x = construct_ids, table = rows$construct_id)]
  }

  normalised <- raw
  for (j in seq_len(length.out = ncol(x = raw))) {
    column <- raw[, j]
    finite <- column[is.finite(x = column)]
    span <- if (length(x = finite) == 0L) {
      c(NA_real_, NA_real_)
    } else {
      range(finite)
    }

    scaled <- if (length(x = finite) == 0L || span[[1]] == span[[2]]) {
      # Nothing to tell the candidates apart on this feature, so it
      # is neutral for all of them rather than arbitrarily decisive.
      rep(x = 0.5, times = length(x = column))
    } else {
      (column - span[[1]]) / (span[[2]] - span[[1]])
    }

    if (available$direction[[j]] == -1L) {
      scaled <- 1 - scaled
    }
    # An unmeasurable feature is neutral, not a penalty.
    scaled[!is.finite(x = scaled)] <- 0.5
    normalised[, j] <- scaled
  }

  list(
    raw = raw,
    normalised = normalised,
    features = available,
    construct_ids = construct_ids
  )
}

#' Resolve the weights to use
#'
#' @param weights User-supplied weights, or NULL.
#' @param features The features available to score.
#' @return A named numeric vector over `features$feature`.
#' @noRd
resolveWeights <- function(weights, features) {
  if (is.null(x = weights)) {
    return(stats::setNames(
      object = rep(x = 1, times = nrow(x = features)),
      nm = features$key
    ))
  }

  if (!is.numeric(x = weights)) {
    stop(
      "`weights` must be a named numeric vector, not an object of ",
      "class \"", class(x = weights)[[1]], "\".",
      call. = FALSE
    )
  }
  if (length(x = weights) == 0L) {
    stop(
      "`weights` must name at least one feature, or be NULL for ",
      "equal weighting.",
      call. = FALSE
    )
  }
  weight_names <- names(x = weights)
  if (is.null(x = weight_names) || any(!nzchar(x = weight_names))) {
    stop(
      "`weights` must be named, so that each weight is attached to a ",
      "feature. The features that can be weighted are: ",
      paste(features$feature, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(x = weight_names) != 0L) {
    stop(
      "`weights` must have unique names, but these are repeated: ",
      paste(
        unique(x = weight_names[duplicated(x = weight_names)]),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  if (anyNA(x = weights) || any(!is.finite(x = weights))) {
    stop("`weights` must all be finite numbers.", call. = FALSE)
  }
  if (any(weights < 0)) {
    stop(
      "`weights` must not be negative, because the score is a ",
      "weighted mean: a negative weight can make the weights sum to ",
      "zero, giving an infinite score, and otherwise pushes scores ",
      "outside 0 to 1 and can invert the whole ordering. Each ",
      "feature already has a direction, so weight 0 drops a feature ",
      "and a larger weight makes it count for more. To score a ",
      "feature the opposite way round, change its direction in ",
      "featureDirections() rather than negating its weight.",
      call. = FALSE
    )
  }
  if (all(weights == 0)) {
    stop("`weights` must include at least one positive weight.",
         call. = FALSE)
  }

  unknown <- weight_names[
    !weight_names %in% features$key & !weight_names %in% features$feature
  ]
  if (length(x = unknown) > 0L) {
    stop(
      if (length(x = unknown) == 1L) "This feature cannot be weighted"
      else "These features cannot be weighted",
      " because it is not scored or not present: ",
      paste(unknown, collapse = ", "),
      ". The features that can be weighted are: ",
      paste(features$key, collapse = ", "),
      ". A bare feature name may be used where it is unambiguous.",
      call. = FALSE
    )
  }

  # Decision: a scored feature left out of a supplied `weights`
  # vector gets a weight of zero, rather than keeping its default
  # weight. Naming a few features is therefore how a caller ranks on
  # just those, which is a common thing to want. The alternative,
  # topping up the unnamed features to equal weight, would mean a
  # caller asking for `c(kozak_score_scaled = 2)` silently got eleven
  # other features at weight 1 as well, and a ranking they did not
  # ask for. The returned `weights` record which features actually
  # carried weight, so either way the result stays traceable.
  resolved <- stats::setNames(
    object = rep(x = 0, times = nrow(x = features)),
    nm = features$key
  )
  already_set <- character(length = 0L)

  for (i in seq_along(along.with = weight_names)) {
    supplied_name <- weight_names[[i]]
    # A region-qualified key names one feature; a bare feature name
    # names every region that feature is scored in, which for
    # gc_content is both of them.
    keys <- if (supplied_name %in% features$key) {
      supplied_name
    } else {
      features$key[features$feature == supplied_name]
    }

    clash <- intersect(x = keys, y = already_set)
    if (length(x = clash) > 0L) {
      stop(
        "`weights` sets ", paste(clash, collapse = ", "),
        " more than once, because \"", supplied_name,
        "\" overlaps a weight already given. Name each feature once, ",
        "using the region-qualified form where a bare name covers ",
        "more than one region.",
        call. = FALSE
      )
    }

    resolved[keys] <- as.numeric(x = weights[[i]])
    already_set <- c(already_set, keys)
  }

  resolved
}

#' Composite scores from a normalised matrix and weights
#'
#' @param normalised The oriented, normalised score matrix.
#' @param weights A named numeric vector over its columns.
#' @return A numeric vector of scores between 0 and 1.
#' @noRd
scoreFromMatrix <- function(normalised, weights) {
  aligned <- weights[colnames(x = normalised)]
  as.numeric(x = normalised %*% aligned) / sum(aligned)
}

#' Rank Candidate Constructs Against One Another
#'
#' @description
#' Turns the tidy feature table into a ranked comparison of candidate
#' designs, with a composite score for each candidate, the
#' contribution every feature made to it, and the weights that
#' produced it.
#'
#' @section What this ranking is and is not:
#' \strong{This is a relative comparison within the submitted set,
#' under weights the caller chooses.} Each feature is normalised across
#' the candidates supplied, so a score says how a candidate compares
#' with the others in this call and nothing more. Adding or removing a
#' candidate changes everyone's scores, and changing the weights
#' changes the order.
#'
#' \strong{It does not predict that any candidate will perform well
#' experimentally.} No part of this package is fitted to measured
#' expression, translation or stability data; the directions it scores
#' on are drawn from published sequence determinants, and a weighted
#' sum of normalised sequence features is a design heuristic, not a
#' model of outcome. Use it to narrow a list of candidates to test,
#' and treat \code{\link{rankStability}} as the check on how much the
#' order depends on the weights you chose.
#'
#' @section Scoring excludes the coding sequence:
#' Rows where \code{region} is \code{"cds"} are filtered out before
#' anything is scored, so no coding sequence feature can reach a
#' score. This is deliberate and is the package's scope: UTRunway
#' compares untranslated region designs around a coding sequence that
#' has already been chosen. Coding sequence features are profiled and
#' plotted as context only. For codon usage and codon optimisation,
#' see the cubar package.
#'
#' Positional rows, such as individual motif occurrences and the
#' sliding window profiles, are also filtered out: ranking uses the
#' per-region summaries, and the counts among them already summarise
#' the positional rows.
#'
#' @section Features scored, and which way round:
#' \tabular{llr}{
#'   \strong{region} \tab \strong{feature} \tab \strong{better when} \cr
#'   utr5 \tab n_upstream_aug         \tab lower \cr
#'   utr5 \tab n_upstream_orf         \tab lower \cr
#'   utr5 \tab n_out_of_frame_overlap \tab lower \cr
#'   utr5 \tab kozak_score_scaled     \tab higher \cr
#'   utr5 \tab gc_content             \tab lower \cr
#'   utr5 \tab max_stem_nt            \tab lower \cr
#'   utr5 \tab n_motif_G_quadruplex   \tab lower \cr
#'   utr3 \tab gc_content             \tab higher \cr
#'   utr3 \tab n_are_pentamer         \tab lower \cr
#'   utr3 \tab n_are_clustered        \tab lower \cr
#'   utr3 \tab n_polya_canonical      \tab higher \cr
#'   utr3 \tab n_mirna_seed           \tab lower
#' }
#'
#' Note that \code{gc_content} is scored in opposite directions in the
#' two untranslated regions, because it stands for different things
#' there: structure that impedes scanning in the 5' UTR, and
#' AU-richness that recruits decay factors in the 3' UTR. Note too
#' that the 3' UTR GC content and the two AU-rich element counts
#' measure overlapping tendencies, the first through composition and
#' the others through a specific motif, so between them they carry
#' three of the twelve default weights.
#'
#' @section Features deliberately left unscored:
#' Some profiled features have no direction this package is willing to
#' assert, and scoring them would dress a guess as a measurement:
#'
#' \itemize{
#'   \item \code{length} in any region: neither longer nor shorter is
#'     better in general.
#'   \item \code{n_motif_TOP}: a TOP motif makes translation
#'     responsive to growth signalling. That is a design decision,
#'     not a quality.
#'   \item \code{n_polya_signal}: counting every hexamer including
#'     variants mixes a functional signal with alternative
#'     polyadenylation, which cuts both ways.
#'   \item \code{max_stem_nt} and \code{n_motif_G_quadruplex} in the
#'     3' UTR: structure there can protect a transcript or help
#'     destabilise it, so neither direction is safe to assert. The
#'     3' UTR's \code{gc_content} \emph{is} scored, as composition
#'     rather than as structure; see the table above.
#'   \item \code{kozak_score} and \code{kozak_strength_ordinal}: the
#'     same quantity as \code{kozak_score_scaled}. Scoring more than
#'     one of them would count initiation context two or three times
#'     over.
#'   \item \code{max_window_gc}: largely the same signal as
#'     \code{gc_content}, so scoring both would double count it.
#' }
#'
#' @section How a score is computed:
#' Each feature is scaled across the candidates to the range 0 to 1 by
#' its own minimum and maximum, then flipped where a lower value is
#' better, so that 1 is always the best of the set on that feature.
#' The score is the weighted mean of those values, so it also lies
#' between 0 and 1, and higher is better. Rank 1 is the highest score,
#' with ties sharing the better rank.
#'
#' Two consequences worth knowing. A feature on which every candidate
#' is identical has no spread to scale, so it is neutral at 0.5 for
#' all of them and does not discriminate. And because scaling uses the
#' minimum and maximum, one extreme candidate compresses the rest
#' towards each other on that feature.
#'
#' @param feature_table A tidy feature table from
#'   \code{\link{profileConstruct}}, holding at least two candidates.
#' @param weights Named numeric vector of weights, or \code{NULL} for
#'   equal weighting across every scored feature. Weights must be
#'   non-negative and are rescaled to sum to one. A scored feature
#'   left out of a supplied vector is given no weight, so naming a few
#'   features is how to rank on just those.
#'
#'   Names may be region-qualified, as \code{"utr5.gc_content"}, or a
#'   bare feature name. A bare name sets every region that feature is
#'   scored in, which matters only for \code{gc_content}, the one
#'   feature scored in both: \code{c(gc_content = 2)} weights it in
#'   both untranslated regions, while \code{c(utr3.gc_content = 2)}
#'   weights it in the 3' UTR alone. The returned weights are always
#'   region-qualified, so the result records exactly which features
#'   carried weight.
#'
#' @return An object of class \code{"utrunway_ranking"}: a list with
#'   \describe{
#'     \item{ranking}{Data frame of \code{construct_id}, \code{score}
#'       and \code{rank}, ordered best first.}
#'     \item{contributions}{Data frame with one row per candidate per
#'       scored feature: \code{construct_id}, \code{region},
#'       \code{feature}, \code{value} as profiled,
#'       \code{normalised} after scaling and orientation,
#'       \code{direction}, \code{weight} as used, and
#'       \code{contribution}, the part of the score this feature
#'       supplied.}
#'     \item{weights}{The weights used, rescaled to sum to one.}
#'     \item{directions}{The direction applied to each scored
#'       feature.}
#'     \item{n_candidates}{How many candidates were compared.}
#'   }
#'   The weights and directions are returned because a ranking is only
#'   reproducible if the choices that produced it travel with it.
#'
#' @references
#' Calvo, S. E., Pervouchine, D. D. and Mootha, V. K. (2009). Upstream
#' open reading frames cause widespread reduction of protein expression
#' and are polymorphic among humans. \emph{Proceedings of the National
#' Academy of Sciences} 106(18), 7507-7512.
#' \doi{10.1073/pnas.0810916106}
#'
#' Kozak, M. (1986). Point mutations define a sequence flanking the AUG
#' initiator codon that modulates translation by eukaryotic ribosomes.
#' \emph{Cell} 44(2), 283-292. \doi{10.1016/0092-8674(86)90762-2}
#'
#' Babendure, J. R., Babendure, J. L., Ding, J.-H. and Tsien, R. Y.
#' (2006). Control of mammalian translation by mRNA structure near
#' caps. \emph{RNA} 12(5), 851-861. \doi{10.1261/rna.2309906}
#'
#' Shaw, G. and Kamen, R. (1986). A conserved AU sequence from the 3'
#' untranslated region of GM-CSF mRNA mediates selective mRNA
#' degradation. \emph{Cell} 46(5), 659-667.
#' \doi{10.1016/0092-8674(86)90341-7}
#'
#' Bartel, D. P. (2009). MicroRNAs: target recognition and regulatory
#' functions. \emph{Cell} 136(2), 215-233.
#' \doi{10.1016/j.cell.2009.01.002}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[1:6, ]
#' )
#'
#' ranking <- rankCandidates(feature_table = features)
#' ranking
#'
#' # The ranked candidates, best first.
#' ranking$ranking
#'
#' # Why the leader leads: its largest feature contributions.
#' leader <- ranking$ranking$construct_id[[1]]
#' contributions <- ranking$contributions[
#'   ranking$contributions$construct_id == leader,
#' ]
#' contributions[
#'   order(contributions$contribution, decreasing = TRUE),
#'   c("feature", "value", "normalised", "contribution")
#' ]
#'
#' # The weights travel with the result, so a ranking can be traced
#' # back to the choices behind it.
#' ranking$weights
#'
#' # Rank on initiation context and upstream ORFs alone.
#' rankCandidates(
#'   feature_table = features,
#'   weights = c(kozak_score_scaled = 2, n_upstream_orf = 1)
#' )$ranking
#'
#' # GC content is scored in both UTRs, so it can be weighted in one
#' # of them by naming the region.
#' rankCandidates(
#'   feature_table = features,
#'   weights = c(utr3.gc_content = 1, n_are_pentamer = 1)
#' )$weights
#'
#' @export
rankCandidates <- function(feature_table, weights = NULL) {
  validateFeatureTable(feature_table = feature_table)

  prepared <- buildScoreMatrix(feature_table = feature_table)
  resolved <- resolveWeights(
    weights = weights, features = prepared$features
  )
  normalised_weights <- resolved / sum(resolved)

  score <- scoreFromMatrix(
    normalised = prepared$normalised, weights = normalised_weights
  )

  ranking <- data.frame(
    construct_id = prepared$construct_ids,
    score = score,
    rank = rank(x = -score, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(ranking$rank, ranking$construct_id), ]
  rownames(x = ranking) <- NULL

  n_features <- nrow(x = prepared$features)
  n_candidates <- length(x = prepared$construct_ids)
  contributions <- data.frame(
    construct_id = rep(x = prepared$construct_ids, times = n_features),
    region = rep(x = prepared$features$region, each = n_candidates),
    feature = rep(x = prepared$features$feature, each = n_candidates),
    value = as.numeric(x = prepared$raw),
    normalised = as.numeric(x = prepared$normalised),
    direction = rep(
      x = prepared$features$direction, each = n_candidates
    ),
    weight = rep(
      x = normalised_weights[prepared$features$key],
      each = n_candidates
    ),
    stringsAsFactors = FALSE
  )
  contributions$contribution <-
    contributions$normalised * contributions$weight
  rownames(x = contributions) <- NULL

  structure(
    .Data = list(
      ranking = ranking,
      contributions = contributions,
      weights = normalised_weights,
      directions = stats::setNames(
        object = prepared$features$direction,
        nm = prepared$features$key
      ),
      n_candidates = n_candidates
    ),
    class = "utrunway_ranking"
  )
}

#' Check that a feature table has the expected shape
#'
#' @param feature_table Value supplied by the caller.
#' @return Invisibly TRUE; called for the error.
#' @noRd
validateFeatureTable <- function(feature_table) {
  required <- c("construct_id", "region", "feature", "value", "position")

  if (is.null(x = feature_table)) {
    stop(
      "`feature_table` must be a feature table from ",
      "profileConstruct(), not NULL.",
      call. = FALSE
    )
  }
  if (!is.data.frame(x = feature_table)) {
    stop(
      "`feature_table` must be a data frame from profileConstruct(), ",
      "not an object of class \"", class(x = feature_table)[[1]], "\".",
      call. = FALSE
    )
  }
  missing_columns <- setdiff(x = required, y = names(x = feature_table))
  if (length(x = missing_columns) > 0L) {
    stop(
      "`feature_table` must have the columns ",
      paste0("`", required, "`", collapse = ", "), ", but ",
      if (length(x = missing_columns) == 1L) "this one is" else
        "these are",
      " missing: ",
      paste0("`", missing_columns, "`", collapse = ", "),
      ". Was it produced by profileConstruct()?",
      call. = FALSE
    )
  }
  if (nrow(x = feature_table) == 0L) {
    stop("`feature_table` has no rows.", call. = FALSE)
  }
  if (!is.numeric(x = feature_table$value)) {
    stop("`feature_table$value` must be numeric.", call. = FALSE)
  }

  invisible(x = TRUE)
}

#' Print a candidate ranking
#'
#' @param x An object of class \code{"utrunway_ranking"}.
#' @param ... Ignored, for consistency with the generic.
#' @return \code{x}, invisibly.
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[1:4, ]
#' )
#' print(x = rankCandidates(feature_table = features))
#'
#' @export
print.utrunway_ranking <- function(x, ...) {
  cat(
    "Ranking of", x$n_candidates, "candidate constructs on",
    length(x = x$weights), "UTR features\n"
  )
  cat("Relative within this set, under these weights. Not a ",
      "prediction of experimental performance.\n", sep = "")
  cat("\n")
  print(x$ranking)
  cat("\nWeights used:\n")
  print(round(x = x$weights, digits = 4))
  invisible(x = x)
}

#' How Stable a Ranking Is Under the Weights Chosen
#'
#' @description
#' Re-ranks the same candidates many times under randomly perturbed
#' weights and reports the spread of each candidate's rank. A
#' candidate that holds first place across most weightings is a
#' genuinely strong candidate within the set; one whose rank swings
#' widely was an artefact of the particular weights chosen.
#'
#' @details
#' Weights are drawn from a Dirichlet distribution centred on the
#' weights supplied, which keeps every draw non-negative and summing
#' to one, as a set of weights must. \code{concentration} sets how
#' tightly the draws cluster around the centre: a large value
#' perturbs the weights gently, a small value explores widely.
#' Features given a weight of zero stay at zero, since excluding a
#' feature is a decision rather than a quantity to perturb.
#'
#' Normalisation of the features themselves does not depend on the
#' weights, so it is done once and reused across every replicate.
#'
#' @section Reproducibility:
#' This function draws random numbers. Call \code{set.seed()} before
#' it if you need the same answer twice, and report the seed alongside
#' the result.
#'
#' @section What it does and does not establish:
#' This shows how sensitive the ordering is to the weights, which is
#' worth knowing before trusting any ranking built on chosen weights.
#' It says nothing about whether the features or their directions are
#' the right ones, and nothing about experimental performance. A
#' ranking can be perfectly stable and still rest on the wrong
#' assumptions.
#'
#' @param feature_table A tidy feature table from
#'   \code{\link{profileConstruct}}, holding at least two candidates.
#' @param n Single whole number, how many perturbed weightings to
#'   draw.
#' @param weights Named numeric vector of weights to centre the draws
#'   on, or \code{NULL} for equal weighting. Interpreted exactly as in
#'   \code{\link{rankCandidates}}.
#' @param concentration Single positive number setting how tightly the
#'   drawn weights cluster around \code{weights}. Larger is tighter.
#'
#' @return An object of class \code{"utrunway_rank_stability"}: a list
#'   with
#'   \describe{
#'     \item{distribution}{Data frame of \code{construct_id},
#'       \code{rank}, \code{n_times} and \code{proportion}, covering
#'       every candidate and rank combination, including those never
#'       seen, so it can be drawn as a complete grid.}
#'     \item{summary}{Data frame of \code{construct_id},
#'       \code{base_rank} under the unperturbed weights, and the
#'       \code{median_rank}, \code{min_rank}, \code{max_rank} and
#'       \code{proportion_best} across the replicates, ordered by
#'       median rank.}
#'     \item{weights}{The centre weights used, summing to one.}
#'     \item{n_replicates}{How many weightings were drawn.}
#'     \item{concentration}{The concentration used.}
#'   }
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[1:6, ]
#' )
#'
#' set.seed(seed = 1)
#' stability <- rankStability(feature_table = features, n = 50)
#' stability
#'
#' # How often each candidate took first place.
#' stability$summary[, c("construct_id", "base_rank", "proportion_best")]
#'
#' # Wider exploration of the weight space moves ranks around more.
#' set.seed(seed = 1)
#' rankStability(
#'   feature_table = features, n = 50, concentration = 2
#' )$summary
#'
#' @export
rankStability <- function(feature_table, n = 100, weights = NULL,
                          concentration = 50) {
  validateFeatureTable(feature_table = feature_table)
  n <- validateWholeNumber(value = n, arg_name = "n", min_value = 2L)
  if (!is.numeric(x = concentration) ||
        length(x = concentration) != 1L ||
        is.na(x = concentration) ||
        !is.finite(x = concentration) ||
        concentration <= 0) {
    stop(
      "`concentration` must be a single positive number, setting how ",
      "tightly the drawn weights cluster around `weights`.",
      call. = FALSE
    )
  }

  base_ranking <- rankCandidates(
    feature_table = feature_table, weights = weights
  )
  prepared <- buildScoreMatrix(feature_table = feature_table)
  centre <- base_ranking$weights
  construct_ids <- prepared$construct_ids
  n_candidates <- length(x = construct_ids)

  is_used <- centre > 0
  ranks <- matrix(
    data = NA_integer_, nrow = n, ncol = n_candidates,
    dimnames = list(NULL, construct_ids)
  )

  for (replicate in seq_len(length.out = n)) {
    drawn <- centre
    # A Dirichlet draw is independent gamma draws, normalised.
    gamma_draw <- stats::rgamma(
      n = sum(is_used), shape = centre[is_used] * concentration, rate = 1
    )
    if (sum(gamma_draw) == 0) {
      gamma_draw <- rep(x = 1, times = sum(is_used))
    }
    drawn[is_used] <- gamma_draw / sum(gamma_draw)

    score <- scoreFromMatrix(
      normalised = prepared$normalised, weights = drawn
    )
    ranks[replicate, ] <- rank(x = -score, ties.method = "min")
  }

  possible_ranks <- seq_len(length.out = n_candidates)
  distribution <- data.frame(
    construct_id = rep(x = construct_ids, each = n_candidates),
    rank = rep(x = possible_ranks, times = n_candidates),
    stringsAsFactors = FALSE
  )
  distribution$n_times <- vapply(
    X = seq_len(length.out = nrow(x = distribution)),
    FUN = function(i) {
      sum(
        ranks[, distribution$construct_id[[i]]] ==
          distribution$rank[[i]]
      )
    },
    FUN.VALUE = numeric(length = 1L)
  )
  distribution$proportion <- distribution$n_times / n

  base_rank <- base_ranking$ranking$rank[
    match(x = construct_ids, table = base_ranking$ranking$construct_id)
  ]
  summary_table <- data.frame(
    construct_id = construct_ids,
    base_rank = base_rank,
    median_rank = apply(X = ranks, MARGIN = 2L, FUN = stats::median),
    min_rank = apply(X = ranks, MARGIN = 2L, FUN = min),
    max_rank = apply(X = ranks, MARGIN = 2L, FUN = max),
    proportion_best = apply(
      X = ranks, MARGIN = 2L, FUN = function(column) mean(column == 1L)
    ),
    stringsAsFactors = FALSE
  )
  summary_table <- summary_table[
    order(summary_table$median_rank, summary_table$base_rank),
  ]
  rownames(x = summary_table) <- NULL
  rownames(x = distribution) <- NULL

  structure(
    .Data = list(
      distribution = distribution,
      summary = summary_table,
      weights = centre,
      n_replicates = n,
      concentration = concentration
    ),
    class = "utrunway_rank_stability"
  )
}

#' Print a rank stability summary
#'
#' @param x An object of class \code{"utrunway_rank_stability"}.
#' @param ... Ignored, for consistency with the generic.
#' @return \code{x}, invisibly.
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[1:4, ]
#' )
#' set.seed(seed = 1)
#' print(x = rankStability(feature_table = features, n = 20))
#'
#' @export
print.utrunway_rank_stability <- function(x, ...) {
  cat(
    "Rank stability over", x$n_replicates,
    "perturbed weightings (concentration", x$concentration, ")\n"
  )
  cat(
    "Sensitivity of the order to the weights only. Not a prediction ",
    "of experimental performance.\n",
    sep = ""
  )
  cat("\n")
  print(x$summary)
  invisible(x = x)
}
