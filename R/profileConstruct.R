#' Profile Every Construct Into One Tidy Feature Table
#'
#' @description
#' Runs every feature function in the package over every region of
#' every construct and returns the results as a single long table, one
#' row per feature value. This table is the one input that every
#' plotting function and \code{rankCandidates()} takes, so a construct
#' is profiled once and then drawn and scored from the same numbers.
#'
#' @details
#' Rows come grouped by construct in input order, then by region in
#' the order \code{utr5}, \code{cds}, \code{utr3}.
#'
#' \strong{Scoring excludes the coding sequence.} Every row carries a
#' \code{region}, so the rows that may contribute to a score are
#' exactly those where \code{region} is not \code{"cds"}:
#'
#' \preformatted{scoreable <- features[features$region != "cds", ]}
#'
#' Coding sequence rows are present so that a construct map can draw
#' the whole design and a comparison can show its context, and they
#' are deliberately left out of ranking, because the package's scope
#' is untranslated region design. See \code{\link{profileCDS}}.
#'
#' @section Feature naming:
#' \code{value} is always numeric, so that features can be compared
#' and scored arithmetically without parsing. What the number means
#' follows from the name:
#'
#' \itemize{
#'   \item A name beginning \code{n_} is a count, and is present with
#'     a value of zero when nothing was found, so that an absent
#'     feature is distinguishable from an unprofiled one.
#'   \item \code{length}, \code{gc_content}, \code{gc3_content},
#'     \code{kozak_score}, \code{kozak_score_scaled},
#'     \code{max_stem_nt} and \code{max_window_gc} are single summary
#'     values for the region.
#'   \item \code{kozak_strength_ordinal} encodes the categorical
#'     strength as \code{0} for weak, \code{1} for adequate and
#'     \code{2} for strong, and is \code{NA} when the context window
#'     does not fit. The ordering is what makes it scoreable.
#'   \item Every other name is a positional feature: one row per
#'     occurrence, with \code{position} giving its 1-based start
#'     within its own region and \code{value} its width in
#'     nucleotides. The exception is \code{upstream_orf}, whose value
#'     is the length of the open reading frame, and the two window
#'     profiles \code{gc_window} and \code{stem_window}, whose values
#'     are the window's GC content and longest complementary run.
#'   \item A motif is named \code{motif_} followed by the motif's
#'     name, and a miRNA seed site \code{mirna_seed_} followed by the
#'     family, so \code{grepl("^motif_", feature)} selects all motif
#'     occurrences.
#' }
#'
#' @section Positions are within a region, not the construct:
#' \code{position} is always relative to the start of its own region.
#' A feature at position 40 of the \code{cds} is 40 nucleotides into
#' the coding sequence, not 40 nucleotides into the construct. To lay
#' the regions end to end for a whole-construct map, add the lengths
#' of the preceding regions, which are themselves rows of this table.
#'
#' @param constructs Constructs to profile: an object from
#'   \code{\link{readConstruct}}, or anything it accepts, namely a
#'   FASTA file path or a data frame with the columns
#'   \code{construct_id}, \code{utr5}, \code{cds} and \code{utr3}.
#'   Input that has not been validated is passed through
#'   \code{\link{readConstruct}} first.
#' @param ... Optional named arguments passed on to the feature
#'   functions: \code{window}, \code{step}, \code{min_stem} and
#'   \code{min_loop} go to \code{\link{estimateStructure}};
#'   \code{seeds}, \code{polya_end_window} and
#'   \code{are_cluster_fraction} go to
#'   \code{\link{scanDecayElements}}; and \code{motifs} goes to
#'   \code{\link{scanMotifs}}, where it replaces the built-in motifs
#'   in every region. Any other name is an error rather than being
#'   ignored, so that a mistyped argument cannot pass unnoticed.
#'
#' @return A data frame with one row per feature value and these
#'   columns:
#'   \describe{
#'     \item{construct_id}{Character. Which construct the row
#'       describes.}
#'     \item{region}{Character. \code{"utr5"}, \code{"cds"} or
#'       \code{"utr3"}, the region the feature was computed from.}
#'     \item{feature}{Character. The feature's name, as described
#'       above.}
#'     \item{value}{Numeric. The feature's value.}
#'     \item{position}{Integer. 1-based start of a positional feature
#'       within its region, or \code{NA} for a region summary.}
#'   }
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # Profile two candidates.
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[c(1, 5), ]
#' )
#' head(x = features, n = 12)
#'
#' # The region summaries for one construct's 5' UTR.
#' summaries <- features[
#'   features$construct_id == "CONSTRUCT_01" &
#'     features$region == "utr5" &
#'     is.na(features$position),
#'   c("feature", "value")
#' ]
#' summaries
#'
#' # The rows that ranking may use: everything but the CDS.
#' table(features$region)
#' nrow(x = features[features$region != "cds", ])
#'
#' # One feature across candidates, which is what ranking compares.
#' features[
#'   features$feature == "n_upstream_orf",
#'   c("construct_id", "value")
#' ]
#'
#' # Positional features, ready to draw on a map.
#' head(
#'   x = features[
#'     !is.na(features$position) &
#'       features$feature != "gc_window" &
#'       features$feature != "stem_window",
#'   ],
#'   n = 8
#' )
#'
#' @importFrom Biostrings DNAString letterFrequency
#' @export
profileConstruct <- function(constructs, ...) {
  constructs <- readConstruct(x = constructs)

  extra_args <- list(...)
  checkProfileArgs(extra_args = extra_args)

  if (any(constructs$has_ambiguous)) {
    ambiguous_ids <- constructs$construct_id[constructs$has_ambiguous]
    stop(
      "The feature functions require unambiguous sequence, but ",
      if (length(x = ambiguous_ids) == 1L) "this construct contains"
      else "these constructs contain",
      " an N: ",
      paste0("\"", ambiguous_ids, "\"", collapse = ", "),
      ". Drop them, for example with ",
      "constructs[!constructs$has_ambiguous, ], or resolve the ",
      "ambiguous positions.",
      call. = FALSE
    )
  }

  rows <- list()
  for (i in seq_len(length.out = nrow(x = constructs))) {
    construct_id <- constructs$construct_id[[i]]

    rows[[length(x = rows) + 1L]] <- profileUtr5(
      construct_id = construct_id,
      utr5 = constructs$utr5[[i]],
      cds = constructs$cds[[i]],
      extra_args = extra_args
    )
    rows[[length(x = rows) + 1L]] <- profileCdsRegion(
      construct_id = construct_id,
      cds = constructs$cds[[i]],
      extra_args = extra_args
    )
    rows[[length(x = rows) + 1L]] <- profileUtr3(
      construct_id = construct_id,
      utr3 = constructs$utr3[[i]],
      extra_args = extra_args
    )
  }

  features <- do.call(what = rbind, args = rows)
  rownames(x = features) <- NULL
  features
}

#' Names accepted in profileConstruct()'s dots, by destination
#'
#' @return A named list of argument names for each feature function.
#' @noRd
profileArgNames <- function() {
  list(
    estimateStructure = c("window", "step", "min_stem", "min_loop"),
    scanDecayElements = c(
      "seeds", "polya_end_window", "are_cluster_fraction"
    ),
    scanMotifs = "motifs"
  )
}

#' Reject unnamed or unrecognised arguments passed through the dots
#'
#' @param extra_args The result of list(...).
#' @return Invisibly TRUE; called for the error.
#' @noRd
checkProfileArgs <- function(extra_args) {
  accepted <- unlist(x = profileArgNames(), use.names = FALSE)

  if (length(x = extra_args) == 0L) {
    return(invisible(x = TRUE))
  }

  arg_names <- names(x = extra_args)
  if (is.null(x = arg_names) || any(!nzchar(x = arg_names))) {
    stop(
      "Arguments passed through `...` must be named. The names ",
      "accepted are: ", paste(accepted, collapse = ", "), ".",
      call. = FALSE
    )
  }

  unknown <- setdiff(x = arg_names, y = accepted)
  if (length(x = unknown) > 0L) {
    stop(
      if (length(x = unknown) == 1L) "Unknown argument " else
        "Unknown arguments ",
      paste0("`", unknown, "`", collapse = ", "),
      " passed through `...`. The names accepted are: ",
      paste(accepted, collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(x = TRUE)
}

#' Call a feature function with only the dots it accepts
#'
#' @param fun The function to call.
#' @param fun_name Its name, to look up the arguments it accepts.
#' @param fixed_args Arguments always supplied, as a named list.
#' @param extra_args The result of list(...).
#' @return Whatever `fun` returns.
#' @noRd
callWithProfileArgs <- function(fun, fun_name, fixed_args, extra_args) {
  permitted <- profileArgNames()[[fun_name]]
  passed <- extra_args[intersect(x = names(x = extra_args), y = permitted)]
  do.call(what = fun, args = c(fixed_args, passed))
}

#' Build feature rows
#'
#' @param construct_id Construct the rows belong to.
#' @param region Region the rows were computed from.
#' @param feature Feature names.
#' @param value Feature values.
#' @param position Positions, or NA for region summaries.
#' @return A data frame of feature rows.
#' @noRd
featureRows <- function(construct_id, region, feature, value,
                        position = NA_integer_) {
  if (length(x = feature) == 0L || length(x = value) == 0L) {
    return(NULL)
  }
  if (length(x = feature) != length(x = value)) {
    stop(
      "Internal error building feature rows for \"", construct_id,
      "\": ", length(x = feature), " feature names but ",
      length(x = value), " values.",
      call. = FALSE
    )
  }

  data.frame(
    construct_id = rep(x = construct_id, times = length(x = feature)),
    region = rep(x = region, times = length(x = feature)),
    feature = as.character(x = feature),
    value = as.numeric(x = value),
    position = as.integer(x = position),
    stringsAsFactors = FALSE
  )
}

#' Prefix a possibly empty vector of names
#'
#' `paste0("motif_", character(0))` returns a single "motif_" rather
#' than an empty vector, which would otherwise produce a feature name
#' with no value to go with it. This keeps an empty input empty.
#'
#' @param prefix The prefix to add.
#' @param x A character vector, possibly empty.
#' @return A character vector the same length as `x`.
#' @noRd
prefixFeature <- function(prefix, x) {
  if (length(x = x) == 0L) {
    return(character(length = 0L))
  }
  paste0(prefix, x)
}

#' GC content of a whole sequence
#'
#' @param sequence An upper-case DNA character string of length one.
#' @return A single numeric proportion.
#' @noRd
sequenceGcContent <- function(sequence) {
  as.numeric(x = Biostrings::letterFrequency(
    x = Biostrings::DNAString(x = sequence),
    letters = "GC",
    as.prob = TRUE
  ))
}

#' Rows for the sliding window profile and its summaries
#'
#' @param construct_id Construct the rows belong to.
#' @param region Region being profiled.
#' @param sequence The region's sequence.
#' @param extra_args The result of list(...).
#' @return A list of feature row data frames.
#' @noRd
windowProfileRows <- function(construct_id, region, sequence,
                              extra_args) {
  structure_profile <- callWithProfileArgs(
    fun = estimateStructure,
    fun_name = "estimateStructure",
    fixed_args = list(sequence = sequence),
    extra_args = extra_args
  )

  list(
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = "max_stem_nt",
      value = max(structure_profile$max_stem_nt)
    ),
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = "max_window_gc",
      value = max(structure_profile$gc_content)
    ),
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = rep(
        x = "gc_window", times = nrow(x = structure_profile)
      ),
      value = structure_profile$gc_content,
      position = structure_profile$window_start
    ),
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = rep(
        x = "stem_window", times = nrow(x = structure_profile)
      ),
      value = structure_profile$max_stem_nt,
      position = structure_profile$window_start
    )
  )
}

#' Rows for the motifs found in a region
#'
#' Counts are emitted for every motif scanned, including those with no
#' occurrences, so that an absent motif reads as zero rather than as a
#' missing row.
#'
#' @param construct_id Construct the rows belong to.
#' @param region Region being scanned.
#' @param sequence The region's sequence.
#' @param extra_args The result of list(...).
#' @return A list of feature row data frames.
#' @noRd
motifRows <- function(construct_id, region, sequence, extra_args) {
  hits <- callWithProfileArgs(
    fun = scanMotifs,
    fun_name = "scanMotifs",
    fixed_args = list(sequence = sequence, region = region),
    extra_args = extra_args
  )

  scanned <- if (is.null(x = extra_args$motifs)) {
    if (region == "utr5") c("TOP", "G_quadruplex") else "G_quadruplex"
  } else {
    names(x = extra_args$motifs)
  }

  counts <- vapply(
    X = scanned,
    FUN = function(motif_name) sum(hits$motif == motif_name),
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  )

  list(
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = prefixFeature(prefix = "n_motif_", x = scanned),
      value = counts
    ),
    featureRows(
      construct_id = construct_id,
      region = region,
      feature = prefixFeature(prefix = "motif_", x = hits$motif),
      value = hits$width,
      position = hits$start
    )
  )
}

#' Feature rows for a 5' UTR
#'
#' @param construct_id Construct the rows belong to.
#' @param utr5 The 5' UTR sequence.
#' @param cds The coding sequence, used to complete the Kozak window.
#' @param extra_args The result of list(...).
#' @return A data frame of feature rows.
#' @noRd
profileUtr5 <- function(construct_id, utr5, cds, extra_args) {
  upstream <- findUpstreamORFs(utr5 = utr5, cds = cds)

  # The main start codon's context reaches into the coding sequence,
  # so it is scored on the two joined. It is reported against the
  # 5' UTR because initiation context is a 5' UTR design choice, and
  # so must stay available to scoring.
  kozak <- suppressWarnings(scoreKozak(
    sequence = paste0(utr5, cds), position = nchar(x = utr5) + 1L
  ))
  strength_ordinal <- c(weak = 0, adequate = 1, strong = 2)[
    kozak$strength
  ]

  uorfs <- upstream[upstream$opens_uorf, , drop = FALSE]

  rows <- c(
    list(
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "length", value = nchar(x = utr5)
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "gc_content",
        value = sequenceGcContent(sequence = utr5)
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "n_upstream_aug", value = nrow(x = upstream)
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "n_upstream_orf", value = sum(upstream$opens_uorf)
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "n_out_of_frame_overlap",
        value = sum(upstream$out_of_frame_overlap)
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "kozak_score", value = kozak$score
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "kozak_score_scaled", value = kozak$score_scaled
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = "kozak_strength_ordinal",
        value = if (length(x = strength_ordinal) == 0L) NA_real_ else
          strength_ordinal
      )
    ),
    windowProfileRows(
      construct_id = construct_id, region = "utr5",
      sequence = utr5, extra_args = extra_args
    ),
    motifRows(
      construct_id = construct_id, region = "utr5",
      sequence = utr5, extra_args = extra_args
    ),
    list(
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = rep(x = "upstream_aug", times = nrow(x = upstream)),
        value = rep(x = 3L, times = nrow(x = upstream)),
        position = upstream$aug_position
      ),
      featureRows(
        construct_id = construct_id, region = "utr5",
        feature = rep(x = "upstream_orf", times = nrow(x = uorfs)),
        value = uorfs$uorf_length_nt,
        position = uorfs$aug_position
      )
    )
  )

  do.call(what = rbind, args = rows)
}

#' Feature rows for a coding sequence, as context only
#'
#' @param construct_id Construct the rows belong to.
#' @param cds The coding sequence.
#' @param extra_args The result of list(...).
#' @return A data frame of feature rows.
#' @noRd
profileCdsRegion <- function(construct_id, cds, extra_args) {
  profile <- profileCDS(cds = cds)

  rows <- c(
    list(
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = "length", value = profile$cds_length
      ),
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = "gc_content", value = profile$gc_content
      ),
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = "gc3_content", value = profile$gc3_content
      ),
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = "n_premature_stop",
        value = length(x = profile$premature_stop_positions)
      ),
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = "n_out_of_frame_atg",
        value = profile$n_out_of_frame_atg
      )
    ),
    windowProfileRows(
      construct_id = construct_id, region = "cds",
      sequence = cds, extra_args = extra_args
    ),
    motifRows(
      construct_id = construct_id, region = "cds",
      sequence = cds, extra_args = extra_args
    ),
    list(
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = rep(
          x = "premature_stop",
          times = length(x = profile$premature_stop_positions)
        ),
        value = rep(
          x = 3L,
          times = length(x = profile$premature_stop_positions)
        ),
        position = profile$premature_stop_positions
      ),
      featureRows(
        construct_id = construct_id, region = "cds",
        feature = rep(
          x = "out_of_frame_atg",
          times = length(x = profile$out_of_frame_atg_positions)
        ),
        value = rep(
          x = 3L,
          times = length(x = profile$out_of_frame_atg_positions)
        ),
        position = profile$out_of_frame_atg_positions
      )
    )
  )

  do.call(what = rbind, args = rows)
}

#' Feature rows for a 3' UTR
#'
#' @param construct_id Construct the rows belong to.
#' @param utr3 The 3' UTR sequence.
#' @param extra_args The result of list(...).
#' @return A data frame of feature rows.
#' @noRd
profileUtr3 <- function(construct_id, utr3, extra_args) {
  elements <- callWithProfileArgs(
    fun = scanDecayElements,
    fun_name = "scanDecayElements",
    fixed_args = list(utr3 = utr3),
    extra_args = extra_args
  )

  is_pentamer <- elements$element_type == "ARE_pentamer"
  is_signal <- elements$element_type == "polyA_signal"
  is_seed <- elements$element_type == "miRNA_seed"

  # A positional name that keeps the element's identity: the two
  # polyA hexamers differ in whether they are canonical, and a seed
  # site is only interpretable with its miRNA family.
  element_feature <- character(length = nrow(x = elements))
  if (nrow(x = elements) > 0L) {
    element_feature[is_pentamer] <- "are_pentamer"
    element_feature[is_signal] <- ifelse(
      test = elements$is_canonical[is_signal],
      yes = "polya_signal_canonical",
      no = "polya_signal_variant"
    )
    element_feature[is_seed] <- prefixFeature(
      prefix = "mirna_seed_", x = elements$element_name[is_seed]
    )
  }

  rows <- c(
    list(
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "length", value = nchar(x = utr3)
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "gc_content",
        value = sequenceGcContent(sequence = utr3)
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "n_are_pentamer", value = sum(is_pentamer)
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "n_are_clustered",
        value = sum(elements$in_cluster, na.rm = TRUE)
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "n_polya_signal", value = sum(is_signal)
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "n_polya_canonical",
        value = sum(
          elements$is_canonical & elements$near_three_prime_end,
          na.rm = TRUE
        )
      ),
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = "n_mirna_seed", value = sum(is_seed)
      )
    ),
    windowProfileRows(
      construct_id = construct_id, region = "utr3",
      sequence = utr3, extra_args = extra_args
    ),
    motifRows(
      construct_id = construct_id, region = "utr3",
      sequence = utr3, extra_args = extra_args
    ),
    list(
      featureRows(
        construct_id = construct_id, region = "utr3",
        feature = element_feature,
        value = elements$width,
        position = elements$start
      )
    )
  )

  do.call(what = rbind, args = rows)
}
