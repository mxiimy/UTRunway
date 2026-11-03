#' Draw One Construct As An Annotated Map
#'
#' A function that draws a single candidate construct on one shared
#' coordinate axis as three stacked tracks: the three regions and
#' their boundaries, every annotated feature, and the sequence-only
#' structure profile. Positions in the feature table are relative to
#' each region, so the regions are laid end to end using the region
#' lengths that are themselves rows of the table. Within the
#' annotation track each class of feature occupies its own row, named
#' on the axis, so identity rests on position and label rather than
#' on colour: the six feature colours were measured pairwise under
#' simulated colour vision deficiency and do not all separate when
#' drawn together. The structure track is a proxy, not a folding
#' prediction.
#'
#' @param feature_table A data frame produced by profileConstruct().
#'
#' @param construct_id A character vector of length one naming the
#'    construct to draw. It must appear in feature_table.
#'
#' @param show_cds_context A logical of length one. Whether to give
#'    coding sequence features their own rows in the annotation
#'    track. The map is a figure about untranslated region design and
#'    coding sequence features are never scored, so they stay out of
#'    it unless asked for; the region bar still shows where the
#'    coding sequence begins and ends. Setting it to TRUE adds rows
#'    for premature in-frame stop codons and out-of-frame ATGs, drawn
#'    in grey and named as context. Default is FALSE.
#'
#' @return Returns a patchwork object, which prints as the assembled
#'    figure and can be customised further like any ggplot2 object.
#'    The individual tracks are accessible by position, so the second
#'    element is the annotation track.
#'
#' @examples
#' # CONSTRUCT_03 carries three upstream AUGs, two of which open
#' # uORFs, and a premature stop codon in its coding sequence.
#' features <- UTRunway::profileConstruct(
#'    constructs = UTRunway::simulatedConstructs[c(3, 12), ])
#' UTRunway::plotConstructMap(feature_table = features,
#'                            construct_id = "CONSTRUCT_03")
#'
#' # The premature stop is not drawn unless asked for, since nothing
#' # scores it.
#' UTRunway::plotConstructMap(feature_table = features,
#'                            construct_id = "CONSTRUCT_03",
#'                            show_cds_context = TRUE)
#'
#' @references
#' Hinnebusch, A. G., Ivanov, I. P. and Sonenberg, N. (2016).
#' Translational control by 5'-untranslated regions of eukaryotic
#' mRNAs. \emph{Science}, 352(6292), 1413-1416.
#' \href{https://doi.org/10.1126/science.aad9868}{Link}
#'
#' Mayr, C. (2019). What are 3' UTRs doing? \emph{Cold Spring Harbor
#' Perspectives in Biology}, 11(10), a034728.
#' \href{https://doi.org/10.1101/cshperspect.a034728}{Link}
#'
#' Pedersen, T. L. (2025). \emph{patchwork: The Composer of Plots}.
#' R package version 1.3.2.
#' \href{https://CRAN.R-project.org/package=patchwork}{Link}
#'
#' @export
#' @importFrom ggplot2 .data
plotConstructMap <- function(feature_table, construct_id,
                             show_cds_context = FALSE) {
  validateFeatureTable(feature_table = feature_table)
  construct_id <- validateConstructChoice(
    feature_table = feature_table, construct_id = construct_id
  )
  if (!is.logical(x = show_cds_context) ||
        length(x = show_cds_context) != 1L ||
        is.na(x = show_cds_context)) {
    stop(
      "`show_cds_context` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  rows <- feature_table[
    feature_table$construct_id == construct_id, , drop = FALSE
  ]
  offsets <- regionOffsets(rows = rows)

  region_track <- buildRegionTrack(offsets = offsets)
  annotation_track <- buildAnnotationTrack(
    rows = rows, offsets = offsets,
    show_cds_context = show_cds_context
  )
  structure_track <- buildStructureTrack(
    rows = rows, offsets = offsets
  )

  # Only the bottom track carries the axis: the tracks share one
  # coordinate system, so repeating its labels would add ink without
  # adding information.
  assembled <- patchwork::wrap_plots(
    region_track + hideXAxis(),
    annotation_track + hideXAxis(),
    structure_track,
    ncol = 1L,
    heights = c(0.32, 1.5, 1)
  )

  assembled + patchwork::plot_annotation(
    title = paste("Construct map:", construct_id),
    subtitle = paste(
      "All three regions on one axis, shown at the foot.",
      "Each feature class has its own labelled row."
    ),
    caption = paste(
      if (show_cds_context) {
        paste(
          "Rows marked as context are coding sequence features. They",
          "take no part in any score or ranking:\n"
        )
      } else {
        paste(
          "Coding sequence features are omitted. They take no part",
          "in any score or ranking:\n"
        )
      },
      "UTRunway compares untranslated region designs around a coding",
      "sequence already chosen.\nThe structure track is a",
      "sequence-only proxy, not a folding prediction."
    ),
    theme = utrunwayTheme()
  )
}

#' Check a requested construct is in the table
#'
#' @param feature_table A validated feature table.
#' @param construct_id Value supplied by the caller.
#' @return The construct identifier, as a character string.
#' @noRd
validateConstructChoice <- function(feature_table, construct_id) {
  available <- unique(x = feature_table$construct_id)

  if (is.null(x = construct_id)) {
    stop(
      "`construct_id` must name one construct, not NULL. The table ",
      "holds: ", paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.character(x = construct_id) ||
        length(x = construct_id) != 1L) {
    stop(
      "`construct_id` must be a single character string naming one ",
      "construct, but it has class \"",
      class(x = construct_id)[[1]], "\" and length ",
      length(x = construct_id), ".",
      call. = FALSE
    )
  }
  if (is.na(x = construct_id) || !construct_id %in% available) {
    stop(
      "`construct_id` \"", construct_id, "\" is not in ",
      "`feature_table`, which holds: ",
      paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }

  construct_id
}

#' Where each region starts in whole-construct coordinates
#'
#' Positions in the feature table are relative to their own region,
#' so drawing them on one axis means offsetting each region by the
#' lengths of those before it.
#'
#' @param rows The feature rows for one construct.
#' @return A data frame of region, length, start and end.
#' @noRd
regionOffsets <- function(rows) {
  regions <- regionLevels()
  lengths <- vapply(
    X = regions,
    FUN = function(region) {
      value <- rows$value[
        rows$region == region & rows$feature == "length"
      ]
      if (length(x = value) != 1L) {
        stop(
          "`feature_table` has no single length for the ", region,
          " of this construct, so the regions cannot be laid out. ",
          "Was it produced by profileConstruct()?",
          call. = FALSE
        )
      }
      value
    },
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  )

  starts <- cumsum(c(0, lengths[-length(x = lengths)]))

  data.frame(
    region = regions,
    length = lengths,
    start = starts,
    end = starts + lengths,
    stringsAsFactors = FALSE
  )
}

#' Convert region-relative positions to construct coordinates
#'
#' @param rows Feature rows carrying region and position.
#' @param offsets The region offset table.
#' @return A numeric vector of absolute positions.
#' @noRd
absolutePosition <- function(rows, offsets) {
  offsets$start[match(x = rows$region, table = offsets$region)] +
    rows$position
}

#' A shared x scale, so every track lines up
#'
#' @param offsets The region offset table.
#' @return A list of ggplot2 scale and coordinate components.
#' @noRd
sharedXScale <- function(offsets) {
  total <- max(offsets$end)
  list(
    ggplot2::scale_x_continuous(
      limits = c(0, total),
      expand = ggplot2::expansion(mult = 0.01)
    ),
    ggplot2::geom_vline(
      xintercept = offsets$start[-1],
      colour = utrunwayPalette()[["muted"]],
      linetype = "dashed",
      linewidth = 0.3
    )
  )
}

#' The region bar: where each region begins and ends
#'
#' @param offsets The region offset table.
#' @return A ggplot object.
#' @noRd
buildRegionTrack <- function(offsets) {
  palette <- utrunwayPalette()
  labels <- utrunwayLabels()

  offsets$fill <- ifelse(
    test = offsets$region == "cds",
    yes = palette[["region_cds"]],
    no = palette[["region_utr"]]
  )
  offsets$label <- labels[offsets$region]
  offsets$midpoint <- (offsets$start + offsets$end) / 2

  ggplot2::ggplot(data = offsets) +
    # A small gap between the fills, so adjacent regions read as
    # separate blocks rather than one continuous bar.
    ggplot2::geom_rect(
      mapping = ggplot2::aes(
        xmin = .data$start + 1, xmax = .data$end - 1,
        ymin = 0, ymax = 1, fill = .data$fill
      )
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = .data$midpoint, y = 0.5,
        label = paste0(.data$label, " (", .data$length, " nt)")
      ),
      colour = palette[["ink"]],
      size = 3
    ) +
    ggplot2::scale_fill_identity() +
    sharedXScale(offsets = offsets) +
    ggplot2::labs(x = NULL, y = NULL) +
    utrunwayTheme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

#' The rows the annotation track can hold, in 5' to 3' order
#'
#' Keyed by the feature name in the tidy table, or by a prefix where
#' a family of features shares a row. Order here is the order the
#' rows appear, top to bottom, so a figure of one construct can be
#' compared with a figure of another.
#'
#' @return A data frame of feature, row label, palette key and
#'   whether the row is coding sequence context.
#' @noRd
annotationRowSpec <- function() {
  labels <- utrunwayLabels()

  data.frame(
    feature = c(
      "upstream_aug_strong", "upstream_aug_adequate",
      "upstream_aug_weak", "upstream_aug_unscored",
      "upstream_orf", "motif_TOP", "motif_G_quadruplex",
      "premature_stop", "out_of_frame_atg",
      "polya_signal_canonical", "polya_signal_variant",
      "are_pentamer", "mirna_seed"
    ),
    label = c(
      labels[["upstream_aug_strong"]],
      labels[["upstream_aug_adequate"]],
      labels[["upstream_aug_weak"]],
      labels[["upstream_aug_unscored"]],
      labels[["upstream_orf"]],
      labels[["top_motif"]],
      labels[["g_quadruplex"]],
      "CDS context: premature stop",
      "CDS context: out-of-frame ATG",
      labels[["polya_signal_canonical"]],
      labels[["polya_signal_variant"]],
      labels[["are_pentamer"]],
      labels[["mirna_seed"]]
    ),
    colour_key = c(
      "upstream_aug_strong", "upstream_aug_adequate",
      "upstream_aug_weak", "upstream_aug_unscored",
      "upstream_aug_adequate", "top_motif", "g_quadruplex",
      "context", "context",
      "polya_signal", "polya_signal",
      "are_pentamer", "mirna_seed"
    ),
    shape_key = c(
      "upstream_aug", "upstream_aug", "upstream_aug", "upstream_aug",
      "upstream_orf", "top_motif", "g_quadruplex",
      "premature_stop", "out_of_frame_atg",
      "polya_signal_canonical", "polya_signal_variant",
      "are_pentamer", "mirna_seed"
    ),
    is_context = c(
      rep(x = FALSE, times = 7L), TRUE, TRUE,
      rep(x = FALSE, times = 4L)
    ),
    stringsAsFactors = FALSE
  )
}

#' Shapes for the rows of the annotation track
#'
#' @return A named integer vector of pch values.
#' @noRd
annotationShapes <- function() {
  shapes <- utrunwayShapes()
  c(
    shapes,
    premature_stop = 4L,
    out_of_frame_atg = 3L
  )
}

#' Gather every positional feature into one table of marks
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @param show_cds_context Whether to keep coding sequence rows.
#' @return A data frame of marks, possibly with zero rows.
#' @noRd
annotationMarks <- function(rows, offsets, show_cds_context) {
  spec <- annotationRowSpec()
  if (!show_cds_context) {
    spec <- spec[!spec$is_context, , drop = FALSE]
  }

  positional <- rows[!is.na(rows$position), , drop = FALSE]
  # miRNA seed rows carry the family in the name, so they are matched
  # by prefix and all share one row.
  positional$row_key <- ifelse(
    test = grepl(pattern = "^mirna_seed_", x = positional$feature),
    yes = "mirna_seed",
    no = positional$feature
  )
  marks <- positional[positional$row_key %in% spec$feature, ,
                      drop = FALSE]
  if (nrow(x = marks) == 0L) {
    return(marks)
  }

  matched <- match(x = marks$row_key, table = spec$feature)
  marks$label <- spec$label[matched]
  marks$colour <- utrunwayPalette()[spec$colour_key[matched]]
  marks$shape <- annotationShapes()[spec$shape_key[matched]]
  marks$is_context <- spec$is_context[matched]
  # The outline carries the colour for shapes that take no fill,
  # which the two coding sequence markers do. Without this they
  # would draw in ink and read as prominently as a scored feature,
  # whatever their fill said.
  marks$outline <- ifelse(
    test = marks$is_context,
    yes = utrunwayPalette()[["context"]],
    no = utrunwayPalette()[["ink"]]
  )
  marks$x <- absolutePosition(rows = marks, offsets = offsets)
  marks$xend <- marks$x + marks$value - 1
  marks$midpoint <- (marks$x + marks$xend) / 2

  # Rows appear top to bottom in 5' to 3' order, so the factor is
  # reversed: ggplot2 draws the first level at the bottom.
  present <- spec$label[spec$label %in% marks$label]
  marks$label <- factor(x = marks$label, levels = rev(x = present))
  marks
}

#' The annotation track: every feature, one row per class
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @param show_cds_context Whether to draw coding sequence rows.
#' @return A ggplot object.
#' @noRd
buildAnnotationTrack <- function(rows, offsets, show_cds_context) {
  palette <- utrunwayPalette()
  marks <- annotationMarks(
    rows = rows, offsets = offsets,
    show_cds_context = show_cds_context
  )

  if (nrow(x = marks) == 0L) {
    return(
      emptyTrack(
        offsets = offsets,
        title = "Annotated features: none found in this construct"
      )
    )
  }

  ggplot2::ggplot() +
    sharedXScale(offsets = offsets) +
    # The extent of each feature, which matters for a uORF or a
    # motif and is invisibly short for a hexamer or a start codon.
    ggplot2::geom_segment(
      data = marks,
      mapping = ggplot2::aes(
        x = .data$x, xend = .data$xend,
        y = .data$label, yend = .data$label,
        colour = .data$colour
      ),
      linewidth = 1.5,
      lineend = "butt"
    ) +
    ggplot2::geom_point(
      data = marks,
      mapping = ggplot2::aes(
        x = .data$midpoint, y = .data$label,
        fill = .data$colour, shape = .data$shape,
        colour = .data$outline
      ),
      size = 2.4,
      stroke = 0.4
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Annotated features across the whole construct"
    ) +
    utrunwayTheme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        colour = palette[["ink_secondary"]], size = 7.5
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = palette[["grid"]], linewidth = 0.25
      )
    )
}

#' The structure track: the sequence-only GC and stem proxy
#'
#' Drawn without a categorical colour: it is one series, named by the
#' track title, with the windows holding a possible stem marked by
#' shape. Nothing here depends on telling colours apart.
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @return A ggplot object.
#' @noRd
buildStructureTrack <- function(rows, offsets) {
  palette <- utrunwayPalette()

  windows <- rows[
    !is.na(rows$position) & rows$feature == "gc_window", , drop = FALSE
  ]
  stems <- rows[
    !is.na(rows$position) & rows$feature == "stem_window", ,
    drop = FALSE
  ]

  if (nrow(x = windows) == 0L) {
    return(emptyTrack(
      offsets = offsets, title = "Structure proxy: no windows"
    ))
  }

  windows$x <- absolutePosition(rows = windows, offsets = offsets)
  plot <- ggplot2::ggplot() +
    sharedXScale(offsets = offsets) +
    ggplot2::geom_hline(
      yintercept = 0.5, colour = palette[["grid"]], linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = windows,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$value, group = .data$region
      ),
      colour = palette[["ink_secondary"]],
      linewidth = 0.6
    )

  # The longest complementary run per window, marked only where it is
  # long enough to be worth a look, since short runs are everywhere
  # in random sequence.
  if (nrow(x = stems) > 0L) {
    stems$x <- absolutePosition(rows = stems, offsets = offsets)
    notable <- stems[stems$value >= 8, , drop = FALSE]
    if (nrow(x = notable) > 0L) {
      notable$gc <- windows$value[
        match(x = notable$x, table = windows$x)
      ]
      plot <- plot +
        ggplot2::geom_point(
          data = notable,
          mapping = ggplot2::aes(x = .data$x, y = .data$gc),
          shape = utrunwayShapes()[["hairpin"]],
          fill = palette[["surface"]],
          colour = palette[["ink"]],
          size = 2.4,
          stroke = 0.5
        )
    }
  }

  plot +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = c(0, 0.5, 1)
    ) +
    ggplot2::labs(
      x = "Position in construct (nt)", y = "GC",
      title = "Structure proxy: sliding-window GC content",
      subtitle = paste(
        "Open triangles mark windows holding a complementary run of",
        "8 nt or more. A proxy, not a folding prediction."
      )
    ) +
    utrunwayTheme()
}

#' A placeholder track, for a region with nothing to draw
#'
#' @param offsets The region offset table.
#' @param title Title explaining what is absent.
#' @return A ggplot object.
#' @noRd
emptyTrack <- function(offsets, title) {
  ggplot2::ggplot() +
    sharedXScale(offsets = offsets) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    utrunwayTheme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
}
