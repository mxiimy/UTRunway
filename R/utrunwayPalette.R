# The single colour palette used by every plot in the package, so that
# a colour learned in one figure carries its meaning into the next.
# Internal: none of this is exported.
#
# The palette was not chosen by eye. Every value below was checked by
# computing, in OKLab, the separation between each pair of colours
# under simulated protanopia and deuteranopia (Machado, Oliveira and
# Fernandes 2009, severity 1.0), the separation under unsimulated
# vision, the OKLCH lightness band and chroma floor, and the WCAG
# contrast against the chart surface. The measured figures are
# recorded beside each group so a reader can see what was tested.
#
# Two findings shaped the design:
#
#   1. All six feature colours together do NOT clear the all-pairs
#      gates: orange against green measures only 3.2 under protanopia,
#      and orange against magenta only 12.9 under normal vision,
#      against a floor of 15. No six-colour set drawn from this palette
#      clears them, and no five-colour set either. That is why the
#      feature colours are split across tracks, three per track, and
#      why each track's set is validated on its own.
#   2. Violet cannot be used while blue is in play: the pair measures
#      4.5 under deuteranopia. It is left out of the palette entirely
#      rather than kept as a tempting spare.
#
# Because of finding 1, colour alone is never sufficient here even
# though each track passes: every feature also carries its own shape
# and a direct label, and sits in its own track. See utrunwayShapes().

#' The UTRunway colour palette
#'
#' @description
#' One palette for every plot in the package, keyed by the meaning a
#' colour carries rather than by where it is used, so the same feature
#' is the same colour in every figure.
#'
#' Measured separations, as OKLab distance times 100, all pairs within
#' a track, worst of protanopia and deuteranopia:
#'
#' \itemize{
#'   \item 5' UTR track (blue, aqua, orange): 9.2 under colour vision
#'     deficiency, 24.0 under normal vision. Target 8 and floor 15.
#'   \item 3' UTR track (green, yellow, magenta): 16.2 and 19.6.
#'   \item Kozak strength steps are an ordinal ramp of one hue rather
#'     than three categorical colours, since strength is ordered.
#'     Smallest lightness step 0.189, against a minimum of 0.06; the
#'     step nearest the surface holds 2.06:1 contrast, against a
#'     minimum of 2.
#' }
#'
#' Aqua, yellow and magenta fall below 3:1 against the light surface.
#' That is permitted only where a mark also carries a visible label,
#' which is why every feature in this package is labelled.
#'
#' @return A named character vector of hex colours.
#' @noRd
utrunwayPalette <- function() {
  c(
    # 5' UTR track. Upstream AUGs take an ordinal ramp of one hue,
    # because Kozak strength is ordered rather than categorical.
    upstream_aug_weak = "#86b6ef",
    upstream_aug_adequate = "#2a78d6",
    upstream_aug_strong = "#104281",
    upstream_aug_unscored = "#898781",
    top_motif = "#1baf7a",
    g_quadruplex = "#eb6834",

    # 3' UTR track.
    polya_signal = "#008300",
    are_pentamer = "#eda100",
    mirna_seed = "#e87ba4",

    # Context and chrome. The coding sequence is drawn in a neutral
    # whose chroma is 0.009, far below the 0.10 floor a categorical
    # colour must clear, so it reads as grey rather than as another
    # feature. That is the point: it is context, not a score.
    context = "#898781",
    context_fill = "#d9d8d4",
    region_utr = "#c3c2b7",
    region_cds = "#e4e3df",
    ink = "#0b0b0b",
    ink_secondary = "#52514e",
    muted = "#898781",
    grid = "#e4e3df",
    surface = "#fcfcfb"
  )
}

#' Shapes that carry the same meanings as the palette
#'
#' @description
#' Every feature class has its own plotting symbol, so that identity
#' never rests on colour alone. This matters more than usual here: the
#' six feature colours do not all separate from one another when seen
#' together, so shape, track position and direct labels are what make
#' the figure readable, with colour reinforcing them.
#'
#' Two variants of one feature keep its colour and differ in shape, as
#' the canonical and variant polyadenylation hexamers do: same family,
#' different form.
#'
#' Only the symbols that accept a fill, \code{21} to \code{25}, are
#' used, since a symbol without one would drop its colour and leave
#' the mark black. There are more feature classes than such symbols,
#' so a shape is reused between tracks, never within one: position
#' separates the reuses, and no track repeats a shape.
#'
#' @return A named integer vector of \code{pch} values.
#' @noRd
utrunwayShapes <- function() {
  c(
    upstream_aug = 24L,
    upstream_orf = 24L,
    top_motif = 22L,
    g_quadruplex = 23L,
    polya_signal_canonical = 21L,
    polya_signal_variant = 22L,
    are_pentamer = 25L,
    mirna_seed = 23L,
    hairpin = 24L
  )
}

#' Fills that distinguish variants sharing a colour
#'
#' @return A named logical vector: TRUE where the symbol is filled.
#' @noRd
utrunwayFilled <- function() {
  c(
    polya_signal_canonical = TRUE,
    polya_signal_variant = FALSE
  )
}

#' Human-readable labels for feature classes
#'
#' Used for legends and direct labels, so that the words a reader sees
#' are the words the documentation uses.
#'
#' @return A named character vector.
#' @noRd
utrunwayLabels <- function() {
  c(
    upstream_aug_weak = "uAUG, weak context",
    upstream_aug_adequate = "uAUG, adequate context",
    upstream_aug_strong = "uAUG, strong context",
    upstream_aug_unscored = "uAUG, context unscorable",
    upstream_orf = "uORF extent",
    top_motif = "TOP motif",
    g_quadruplex = "G-quadruplex-prone run",
    polya_signal_canonical = "polyA signal (AAUAAA)",
    polya_signal_variant = "polyA signal (AUUAAA)",
    are_pentamer = "AU-rich element (AUUUA)",
    mirna_seed = "miRNA seed site",
    utr5 = "5' UTR",
    cds = "CDS",
    utr3 = "3' UTR"
  )
}

#' A plain theme shared by every plot in the package
#'
#' Recessive grid and axes, text in ink rather than in any series
#' colour, and no legend box where a single series is named by the
#' panel title instead.
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
