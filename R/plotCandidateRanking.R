#' Plot The Ranking With How Much It Depends On The Weights
#'
#' A function that draws each candidate's rank as a point at its
#' median across many perturbed weightings, with an interval spanning
#' the ranks it took. A tight interval means the conclusion survives
#' the choice of weights; a wide one is a warning that the order is
#' an artefact of them. The open marker is the rank under the
#' unperturbed weights. Read the intervals before the order: a
#' candidate whose interval spans half the field has no settled
#' position, and two candidates whose intervals overlap heavily are
#' not meaningfully separated. A tight interval establishes
#' robustness to the weights only, not that the features or their
#' directions are right.
#'
#' @param stability An object of class
#'    "utrunway_rank_stability" produced by rankStability().
#'
#' @return Returns a ggplot2 object showing each candidate's median
#'    rank, the full range of ranks taken, and the rank under the
#'    weights as given, which can be customised further.
#'
#' @examples
#' features <- UTRunway::profileConstruct(
#'    constructs = UTRunway::simulatedConstructs[1:6, ])
#'
#' set.seed(seed = 1)
#' UTRunway::plotCandidateRanking(
#'    stability = UTRunway::rankStability(feature_table = features,
#'                                        n = 100))
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data
#' Analysis}. Springer-Verlag New York.
#' \href{https://ggplot2.tidyverse.org}{Link}
#'
#' @export
#' @importFrom ggplot2 .data
plotCandidateRanking <- function(stability) {
  validateStabilityObject(stability = stability)
  palette <- utrunwayPalette()

  summary_table <- stability$summary
  n_candidates <- nrow(x = summary_table)

  # Best at the top. Ties on the median are settled by the
  # unperturbed rank, so the order is deterministic.
  ordered_ids <- summary_table$construct_id[
    order(
      summary_table$median_rank, summary_table$base_rank,
      decreasing = TRUE
    )
  ]
  summary_table$construct_id <- factor(
    x = summary_table$construct_id, levels = ordered_ids
  )
  summary_table$span <- summary_table$max_rank -
    summary_table$min_rank
  # Labelled selectively: a row of "0% first" against every
  # candidate that never led would be noise.
  summary_table$best_label <- ifelse(
    test = summary_table$proportion_best > 0,
    yes = paste0(
      round(x = summary_table$proportion_best * 100), "% first"
    ),
    no = NA_character_
  )

  ggplot2::ggplot(
    data = summary_table,
    mapping = ggplot2::aes(y = .data$construct_id)
  ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$min_rank, xend = .data$max_rank,
        yend = .data$construct_id
      ),
      colour = palette[["muted"]],
      linewidth = 2.4,
      lineend = "round",
      alpha = 0.55
    ) +
    # The unperturbed rank, as an open marker behind the median.
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$base_rank),
      shape = 21L,
      fill = palette[["surface"]],
      colour = palette[["ink_secondary"]],
      size = 3.4,
      stroke = 0.7
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$median_rank),
      shape = 21L,
      fill = palette[["upstream_aug_adequate"]],
      colour = palette[["ink"]],
      size = 3,
      stroke = 0.4
    ) +
    ggplot2::geom_text(
      data = summary_table[!is.na(summary_table$best_label), ],
      mapping = ggplot2::aes(
        x = .data$max_rank, label = .data$best_label
      ),
      hjust = -0.25,
      size = 2.7,
      colour = palette[["ink_secondary"]]
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(length.out = n_candidates),
      limits = c(0.5, n_candidates + 1.4),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::labs(
      x = "Rank (1 is best)",
      y = NULL,
      title = "Ranking, and how far it moves with the weights",
      subtitle = paste0(
        "Filled point: median rank over ", stability$n_replicates,
        " perturbed weightings. Open point: rank under the",
        "\nweights as given. Bar: the full range of ranks taken.",
        " A wide bar means the order is not settled."
      ),
      caption = paste(
        "Sensitivity to the weights only. It does not establish that",
        "the features or their directions are\nright, and it is not",
        "a prediction of experimental performance."
      )
    ) +
    utrunwayTheme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        colour = palette[["ink"]]
      )
    )
}
