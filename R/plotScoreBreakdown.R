#' Show What Each Candidate's Score Is Made Of
#'
#' A function that draws the composite score of every candidate as a
#' stack of the per-feature contributions behind it, split by region,
#' so a ranking can be read as a set of reasons rather than as one
#' number. The height of a stack is the candidate's score and each
#' segment is one feature's contribution, so a short segment is a
#' feature on which the candidate did poorly relative to the others.
#' The dashed rule in each panel is the most that region could
#' contribute, so the gap between a stack and the rule is what the
#' candidate gave away. Only the untranslated regions appear, because
#' only they are scored.
#'
#' @param ranking An object of class "utrunway_ranking" produced by
#'    rankCandidates().
#'
#' @return Returns a ggplot2 object of horizontal stacked bars
#'    faceted by region, ordered best candidate first, which can be
#'    customised further.
#'
#' @examples
#' features <- UTRunway::profileConstruct(
#'    constructs = UTRunway::simulatedConstructs[1:6, ])
#' ranking <- UTRunway::rankCandidates(feature_table = features)
#'
#' UTRunway::plotScoreBreakdown(ranking = ranking)
#'
#' # A different weighting gives a different breakdown.
#' UTRunway::plotScoreBreakdown(
#'    ranking = UTRunway::rankCandidates(
#'       feature_table = features,
#'       weights = c(kozak_score_scaled = 3, n_upstream_orf = 1)))
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data
#' Analysis}. Springer-Verlag New York.
#' \href{https://ggplot2.tidyverse.org}{Link}
#'
#' @export
#' @importFrom ggplot2 .data
plotScoreBreakdown <- function(ranking) {
  validateRankingObject(ranking = ranking)
  palette <- utrunwayPalette()

  contributions <- ranking$contributions
  contributions <- contributions[
    contributions$weight > 0, , drop = FALSE
  ]
  if (nrow(x = contributions) == 0L) {
    stop(
      "No feature carries any weight in this ranking, so there is ",
      "no breakdown to draw.",
      call. = FALSE
    )
  }

  # Best candidate at the top, so the figure reads in rank order.
  order_by_rank <- ranking$ranking$construct_id[
    order(ranking$ranking$rank, decreasing = TRUE)
  ]
  contributions$construct_id <- factor(
    x = contributions$construct_id, levels = order_by_rank
  )
  contributions$region <- factor(
    x = contributions$region, levels = regionLevels()
  )
  contributions$region <- droplevels(x = contributions$region)

  # Segments are ordered parts of one whole, so they take one hue
  # light to dark rather than a set of categorical colours.
  feature_levels <- unique(x = contributions$feature[
    order(contributions$region, contributions$feature)
  ])
  contributions$feature <- factor(
    x = contributions$feature, levels = feature_levels
  )

  # One weight per feature, then summed per region. Taking unique
  # weights directly would collapse equally weighted features into
  # one, which is exactly what equal weighting makes them.
  per_feature <- unique(
    x = contributions[, c("region", "feature", "weight")]
  )
  ceilings <- stats::aggregate(
    x = per_feature["weight"],
    by = list(region = per_feature$region),
    FUN = sum
  )
  names(x = ceilings)[[2]] <- "ceiling"

  ggplot2::ggplot(
    data = contributions,
    mapping = ggplot2::aes(
      x = .data$contribution, y = .data$construct_id,
      fill = .data$feature
    )
  ) +
    ggplot2::geom_col(
      width = 0.68,
      colour = palette[["surface"]],
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      data = ceilings,
      mapping = ggplot2::aes(xintercept = .data$ceiling),
      colour = palette[["ink_secondary"]],
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$region),
      nrow = 1L,
      labeller = ggplot2::labeller(region = regionFacetLabels())
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        object = stackRamp(n = length(x = feature_levels)),
        nm = feature_levels
      )
    ) +
    ggplot2::labs(
      x = "Contribution to the composite score",
      y = NULL,
      title = "What each candidate's score is made of",
      subtitle = paste(
        "Segments sum to the score. The dashed rule is the most the",
        "region could contribute,\nso the gap to it is what the",
        "candidate gave away."
      ),
      caption = paste0(
        describeWeights(weights = ranking$weights),
        "\nA relative comparison within this set under these",
        " weights, not a prediction of experimental performance."
      )
    ) +
    utrunwayTheme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 7)
    )
}
