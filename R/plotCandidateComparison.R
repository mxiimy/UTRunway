#' Compare Every Candidate Across Every Feature At Once
#'
#' A function that draws each candidate as a line across all of its
#' region summary features, scaled so that features on different
#' units can share one axis. It answers which candidates are unusual,
#' and on what, before any weighting is chosen. Each feature is
#' scaled so 0 is the lowest value in the set and 1 the highest; the
#' scaling says nothing about which end is better, since direction
#' belongs to scoring. Candidates are drawn in a recessive grey
#' unless named, because one colour per candidate would be
#' unreadable and only three colours of this palette separate
#' reliably when all pairs are in play.
#'
#' @param feature_table A data frame produced by profileConstruct(),
#'    holding at least two candidates.
#'
#' @param highlight A character vector of at most three construct
#'    identifiers to bring forward in colour and label, or NULL to
#'    draw every candidate in grey. Default is NULL.
#'
#' @return Returns a ggplot2 object faceted by region, with the
#'    coding sequence panel shaded and labelled as context, which can
#'    be customised further.
#'
#' @examples
#' features <- UTRunway::profileConstruct(
#'    constructs = UTRunway::simulatedConstructs[1:6, ])
#'
#' # The whole field at once.
#' UTRunway::plotCandidateComparison(feature_table = features)
#'
#' # Two candidates brought forward against the rest.
#' UTRunway::plotCandidateComparison(
#'    feature_table = features,
#'    highlight = c("CONSTRUCT_03", "CONSTRUCT_06"))
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data
#' Analysis}. Springer-Verlag New York.
#' \href{https://ggplot2.tidyverse.org}{Link}
#'
#' @export
#' @importFrom ggplot2 .data
plotCandidateComparison <- function(feature_table, highlight = NULL) {
  validateFeatureTable(feature_table = feature_table)
  palette <- utrunwayPalette()

  summaries <- summaryRows(feature_table = feature_table)
  construct_ids <- unique(x = summaries$construct_id)
  if (length(x = construct_ids) < 2L) {
    stop(
      "Comparing candidates needs at least two of them, but ",
      length(x = construct_ids), " was supplied.",
      call. = FALSE
    )
  }
  highlight <- validateHighlight(
    highlight = highlight, construct_ids = construct_ids
  )

  # Scale each feature across the candidates, so that counts and
  # proportions can share one axis.
  summaries$key <- paste(
    summaries$region, summaries$feature, sep = "."
  )
  summaries$scaled <- unsplit(
    value = lapply(
      X = split(x = summaries$value, f = summaries$key),
      FUN = normaliseWithinSet
    ),
    f = summaries$key
  )
  summaries$region <- factor(
    x = summaries$region, levels = regionLevels()
  )
  summaries$is_highlighted <- summaries$construct_id %in% highlight

  background <- summaries[!summaries$is_highlighted, , drop = FALSE]
  foreground <- summaries[summaries$is_highlighted, , drop = FALSE]

  plot <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$feature, y = .data$scaled, group = .data$construct_id
    )
  ) +
    # The coding sequence panel is shaded as well as labelled, so
    # that the part of the figure which takes no part in scoring is
    # obvious without reading the strip.
    ggplot2::geom_rect(
      data = data.frame(
        region = factor(x = "cds", levels = regionLevels())
      ),
      mapping = ggplot2::aes(
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ),
      fill = palette[["context_fill"]],
      alpha = 0.4,
      inherit.aes = FALSE
    )

  if (nrow(x = background) > 0L) {
    plot <- plot +
      ggplot2::geom_line(
        data = background,
        colour = palette[["muted"]],
        linewidth = 0.4,
        alpha = 0.6
      ) +
      ggplot2::geom_point(
        data = background,
        colour = palette[["muted"]],
        size = 1.1,
        alpha = 0.6
      )
  }
  if (nrow(x = foreground) > 0L) {
    plot <- plot +
      ggplot2::geom_line(
        data = foreground,
        mapping = ggplot2::aes(colour = .data$construct_id),
        linewidth = 1,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = foreground,
        mapping = ggplot2::aes(colour = .data$construct_id),
        size = 2
      ) +
      ggplot2::scale_colour_manual(
        values = stats::setNames(
          object = candidateHighlightColours()[
            seq_along(along.with = highlight)
          ],
          nm = highlight
        )
      )
  }

  plot +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$region),
      scales = "free_x",
      nrow = 1L,
      labeller = ggplot2::labeller(region = regionFacetLabels())
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = c(0, 0.5, 1)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Scaled within this set of candidates",
      title = "Every candidate across every feature",
      subtitle = paste(
        "0 is the lowest value among these candidates and 1 the",
        "highest. High is not the same as good:",
        "\ndirection belongs to scoring, not to this figure."
      ),
      caption = contextCaption()
    ) +
    utrunwayTheme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, size = 7
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = palette[["grid"]], linewidth = 0.3
      ),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

#' Check the candidates asked to be highlighted
#'
#' @param highlight Value supplied by the caller.
#' @param construct_ids The candidates present in the table.
#' @return A character vector, possibly empty.
#' @noRd
validateHighlight <- function(highlight, construct_ids) {
  if (is.null(x = highlight)) {
    return(character(length = 0L))
  }
  if (!is.character(x = highlight)) {
    stop(
      "`highlight` must be a character vector of construct ",
      "identifiers, not an object of class \"",
      class(x = highlight)[[1]], "\".",
      call. = FALSE
    )
  }
  if (length(x = highlight) > 3L) {
    stop(
      "`highlight` takes at most three candidates, because only ",
      "three colours of this palette separate reliably under ",
      "simulated colour vision deficiency when all pairs are in ",
      "play. ", length(x = highlight), " were given.",
      call. = FALSE
    )
  }
  unknown <- setdiff(x = highlight, y = construct_ids)
  if (length(x = unknown) > 0L) {
    stop(
      "`highlight` names ",
      if (length(x = unknown) == 1L) "a construct" else "constructs",
      " not in `feature_table`: ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(x = highlight) != 0L) {
    stop("`highlight` must not repeat a construct.", call. = FALSE)
  }

  highlight
}
