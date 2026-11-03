# Theming and small pieces shared by every plot in the package, so
# that the figures look like one family and a reader learns the
# conventions once. Internal: none of this is exported.

#' A plain theme shared by every plot in the package
#'
#' Recessive grid and axes, text in ink rather than in any series
#' colour, and a caption style that can carry the warnings the
#' figures need to carry.
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @noRd
utrunwayTheme <- function(base_size = 11) {
  palette <- utrunwayPalette()

  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = palette[["surface"]], colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = palette[["surface"]], colour = NA
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = palette[["grid"]], linewidth = 0.3
      ),
      axis.text = ggplot2::element_text(colour = palette[["muted"]]),
      axis.title = ggplot2::element_text(
        colour = palette[["ink_secondary"]]
      ),
      strip.text = ggplot2::element_text(
        colour = palette[["ink"]], size = base_size - 1, hjust = 0
      ),
      plot.title = ggplot2::element_text(
        colour = palette[["ink"]], size = base_size, face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        colour = palette[["ink_secondary"]], size = base_size - 1
      ),
      plot.caption = ggplot2::element_text(
        colour = palette[["ink_secondary"]], size = base_size - 2,
        hjust = 0
      ),
      legend.text = ggplot2::element_text(
        colour = palette[["ink_secondary"]], size = base_size - 2
      ),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
}

#' Drop the x axis from a track that is not the bottom one
#'
#' @return A ggplot2 theme component.
#' @noRd
hideXAxis <- function() {
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
}

#' The three regions, in the order they occur in a transcript
#'
#' @return A character vector of length three.
#' @noRd
regionLevels <- function() {
  c("utr5", "cds", "utr3")
}

#' Facet labels naming each region and whether it is scored
#'
#' The coding sequence is labelled as context everywhere it appears,
#' so the exclusion is visible in the figure rather than only in the
#' documentation.
#'
#' @return A named character vector.
#' @noRd
regionFacetLabels <- function() {
  c(
    utr5 = "5' UTR (scored)",
    cds = "CDS (context, not scored)",
    utr3 = "3' UTR (scored)"
  )
}

#' The caption every figure carries about the coding sequence
#'
#' @return A character string.
#' @noRd
contextCaption <- function() {
  paste(
    "Coding sequence features are shown as context and take no part",
    "in any score or ranking:\nUTRunway compares untranslated region",
    "designs around a coding sequence already chosen."
  )
}

#' Colours for highlighting individual candidates
#'
#' Three, and no more, because only three of this palette's hues
#' separate from one another under simulated colour vision deficiency
#' when all pairs are in play, which is the case when lines overlap.
#' Highlighted candidates are also labelled directly, so the colour
#' is a help rather than the only clue.
#'
#' @return A character vector of three hex colours.
#' @noRd
candidateHighlightColours <- function() {
  palette <- utrunwayPalette()
  unname(obj = c(
    palette[["upstream_aug_adequate"]],
    palette[["top_motif"]],
    palette[["g_quadruplex"]]
  ))
}

#' A single-hue ramp for the segments of a stacked bar
#'
#' Stack segments are ordered parts of one whole, so they take one
#' hue light to dark rather than a set of categorical colours. The
#' endpoints are the documented sequential blue steps nearest the
#' permitted ends of the ramp; intermediate steps are interpolated,
#' which keeps adjacent steps at least 0.06 apart in lightness for up
#' to about eight segments.
#'
#' @param n How many segments are needed.
#' @return A character vector of n hex colours, light to dark.
#' @noRd
stackRamp <- function(n) {
  if (n <= 1L) {
    return("#4976AD")
  }
  grDevices::colorRampPalette(colors = c("#86b6ef", "#0d366b"))(n)
}

#' Scale a vector to 0 to 1 across the candidates supplied
#'
#' A feature on which every candidate is identical has no spread to
#' scale, so it sits at the midpoint for all of them rather than
#' being arbitrarily decisive. This matches how
#' \code{\link{rankCandidates}} normalises.
#'
#' @param values A numeric vector.
#' @return A numeric vector between 0 and 1.
#' @noRd
normaliseWithinSet <- function(values) {
  finite <- values[is.finite(x = values)]
  if (length(x = finite) == 0L || min(finite) == max(finite)) {
    return(rep(x = 0.5, times = length(x = values)))
  }
  scaled <- (values - min(finite)) / (max(finite) - min(finite))
  scaled[!is.finite(x = scaled)] <- 0.5
  scaled
}

#' Check that a ranking object came from rankCandidates()
#'
#' @param ranking Value supplied by the caller.
#' @param arg_name Name of the caller's argument.
#' @return Invisibly TRUE; called for the error.
#' @noRd
validateRankingObject <- function(ranking, arg_name = "ranking") {
  if (is.null(x = ranking)) {
    stop(
      "`", arg_name, "` must be an object from rankCandidates(), ",
      "not NULL.",
      call. = FALSE
    )
  }
  if (!inherits(x = ranking, what = "utrunway_ranking")) {
    stop(
      "`", arg_name, "` must be an object from rankCandidates(), not ",
      "an object of class \"", class(x = ranking)[[1]], "\".",
      call. = FALSE
    )
  }
  invisible(x = TRUE)
}

#' Check that a stability object came from rankStability()
#'
#' @param stability Value supplied by the caller.
#' @param arg_name Name of the caller's argument.
#' @return Invisibly TRUE; called for the error.
#' @noRd
validateStabilityObject <- function(stability,
                                    arg_name = "stability") {
  if (is.null(x = stability)) {
    stop(
      "`", arg_name, "` must be an object from rankStability(), ",
      "not NULL.",
      call. = FALSE
    )
  }
  if (!inherits(x = stability, what = "utrunway_rank_stability")) {
    stop(
      "`", arg_name, "` must be an object from rankStability(), not ",
      "an object of class \"", class(x = stability)[[1]], "\".",
      call. = FALSE
    )
  }
  invisible(x = TRUE)
}

#' The region summary rows of a feature table
#'
#' @param feature_table A validated feature table.
#' @return The rows with no position, as a data frame.
#' @noRd
summaryRows <- function(feature_table) {
  feature_table[is.na(feature_table$position), , drop = FALSE]
}

#' Describe a set of weights in one short line
#'
#' Listing twelve names and numbers runs off the side of a figure, so
#' equal weighting is stated as such, and an uneven set is listed only
#' as far as it fits.
#'
#' @param weights A named numeric vector, summing to one.
#' @return A character string of length one.
#' @noRd
describeWeights <- function(weights) {
  used <- weights[weights > 0]
  if (length(x = used) == 0L) {
    return("No feature carries any weight.")
  }

  spread <- max(used) - min(used)
  if (spread < sqrt(x = .Machine$double.eps)) {
    return(paste0(
      "Weights: equal across all ", length(x = used),
      " scored features."
    ))
  }

  ordered <- sort(x = used, decreasing = TRUE)
  shown <- utils::head(x = ordered, n = 4L)
  described <- paste(
    names(x = shown), round(x = shown, digits = 3),
    sep = " ", collapse = ", "
  )
  if (length(x = ordered) > length(x = shown)) {
    described <- paste0(
      described, ", and ", length(x = ordered) - length(x = shown),
      " more"
    )
  }
  paste0("Weights, largest first: ", described, ".")
}
