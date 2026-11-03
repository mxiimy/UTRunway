#' Plot Composition Along The Whole Construct
#'
#' A function that draws the sliding-window GC content and the
#' structure proxy along the whole construct, with the region
#' boundaries marked, so a compositional problem can be seen where it
#' sits rather than only as a number for the region it falls in. The
#' two measures are drawn in separate panels rather than on two y
#' axes, because GC content is a proportion and the complementary run
#' is a length in nucleotides. The coding sequence GC3 is drawn as a
#' level line spanning the coding sequence because that is what it
#' is, one value for the region rather than a sliding window.
#'
#' @param feature_table A data frame produced by profileConstruct().
#'
#' @param construct_id A character vector of length one naming the
#'    construct to draw. It must appear in feature_table.
#'
#' @return Returns a ggplot2 object of two faceted panels sharing the
#'    construct's coordinate axis, which can be customised further.
#'
#' @examples
#' # CONSTRUCT_05 carries a planted G-quadruplex-prone run, which
#' # shows as a GC-rich stretch in its 5' UTR.
#' features <- UTRunway::profileConstruct(
#'    constructs = UTRunway::simulatedConstructs[c(5, 10), ])
#' UTRunway::plotRegionProfile(feature_table = features,
#'                             construct_id = "CONSTRUCT_05")
#'
#' @references
#' Babendure, J. R., Babendure, J. L., Ding, J.-H. and Tsien, R. Y.
#' (2006). Control of mammalian translation by mRNA structure near
#' caps. \emph{RNA}, 12(5), 851-861.
#' \href{https://doi.org/10.1261/rna.2309906}{Link}
#'
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data
#' Analysis}. Springer-Verlag New York.
#' \href{https://ggplot2.tidyverse.org}{Link}
#'
#' @export
#' @importFrom ggplot2 .data
plotRegionProfile <- function(feature_table, construct_id) {
  validateFeatureTable(feature_table = feature_table)
  construct_id <- validateConstructChoice(
    feature_table = feature_table, construct_id = construct_id
  )
  palette <- utrunwayPalette()

  rows <- feature_table[
    feature_table$construct_id == construct_id, , drop = FALSE
  ]
  offsets <- regionOffsets(rows = rows)

  windows <- rows[
    !is.na(rows$position) &
      rows$feature %in% c("gc_window", "stem_window"),
    ,
    drop = FALSE
  ]
  if (nrow(x = windows) == 0L) {
    stop(
      "`feature_table` holds no sliding window rows for \"",
      construct_id, "\", so there is no profile to draw.",
      call. = FALSE
    )
  }
  windows$x <- absolutePosition(rows = windows, offsets = offsets)
  windows$panel <- ifelse(
    test = windows$feature == "gc_window",
    yes = "GC content of the window",
    no = "Longest complementary run (nt)"
  )

  cds_span <- offsets[offsets$region == "cds", ]
  gc3 <- rows$value[
    rows$region == "cds" & rows$feature == "gc3_content"
  ]
  gc3_line <- data.frame(
    x = cds_span$start,
    xend = cds_span$end,
    y = gc3,
    panel = "GC content of the window",
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot() +
    # The coding sequence is shaded, so the part of the profile that
    # is context rather than score is obvious at a glance.
    ggplot2::geom_rect(
      data = cds_span,
      mapping = ggplot2::aes(
        xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf
      ),
      fill = palette[["context_fill"]],
      alpha = 0.45
    ) +
    ggplot2::geom_vline(
      xintercept = offsets$start[-1],
      colour = palette[["muted"]],
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = windows,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$value, group = .data$region
      ),
      colour = palette[["ink_secondary"]],
      linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      data = gc3_line,
      mapping = ggplot2::aes(
        x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y
      ),
      colour = palette[["context"]],
      linetype = "dotted",
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = gc3_line,
      mapping = ggplot2::aes(
        x = (.data$x + .data$xend) / 2, y = .data$y,
        label = paste0("CDS GC3 ", round(x = .data$y * 100), "%")
      ),
      vjust = -0.7,
      size = 2.8,
      colour = palette[["ink_secondary"]]
    ) +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$panel), ncol = 1L, scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::labs(
      x = "Position in construct (nt)",
      y = NULL,
      title = paste("Composition along", construct_id),
      subtitle = paste(
        "Dashed lines mark the region boundaries; the shaded span is",
        "the coding sequence."
      ),
      caption = paste(
        contextCaption(),
        "\nThe lower panel is a sequence-only proxy, not a folding",
        "prediction, and misses anything spanning more than one",
        "window."
      )
    ) +
    utrunwayTheme()
}
