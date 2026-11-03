#' Draw One Construct as an Annotated Multi-Track Map
#'
#' @description
#' Draws a single candidate construct as stacked tracks on one shared
#' coordinate axis: the three regions and their boundaries, the 5' UTR
#' features that govern initiation, the sequence-only structure
#' profile, the coding sequence as context, and the 3' UTR elements
#' that govern stability. It is the package's signature figure, and
#' the one place a whole design can be read at once.
#'
#' @details
#' Positions in the feature table are relative to each region, so the
#' regions are laid end to end here, using the region lengths that are
#' themselves rows of the table. The x axis is therefore the position
#' in the whole construct, and every track shares it.
#'
#' @section The coding sequence track is context, not score:
#' The coding sequence track is deliberately drawn in grey while the
#' untranslated regions carry colour, and the figure's caption says
#' so. UTRunway scores and ranks on untranslated region features
#' alone, so nothing in the CDS track contributes to a score. Its
#' GC3 profile and its cryptic codons are there to let a design be
#' read as a whole, not to be weighed. See \code{\link{profileCDS}}
#' and \code{\link{rankCandidates}}.
#'
#' @section Colour never carries meaning alone:
#' Each feature class has its own plotting symbol and its own track as
#' well as its own colour, and the legend names every class. This is
#' not decoration. The six feature colours were checked pairwise under
#' simulated protanopia and deuteranopia, and while the three used in
#' any one track separate well, the full six do not separate from one
#' another when seen together. Shape, track and label are what make
#' the figure readable for every reader; colour reinforces them.
#'
#' @section The structure track is a proxy:
#' The structure track shows sliding-window GC content with the
#' windows holding a possible hairpin stem marked. This is a
#' sequence-only proxy, not a folding prediction: see
#' \code{\link{estimateStructure}} for what it does and does not
#' establish.
#'
#' @param feature_table A tidy feature table from
#'   \code{\link{profileConstruct}}.
#' @param construct_id Character string of length one naming the
#'   construct to draw. It must appear in \code{feature_table}.
#'
#' @return A \pkg{patchwork} object, which prints as the assembled
#'   figure and can be customised further like any
#'   \pkg{ggplot2} object, for example by adding a theme or a title.
#'   The individual tracks are accessible by position, so
#'   \code{map[[2]]} is the 5' UTR track.
#'
#' @references
#' Hinnebusch, A. G., Ivanov, I. P. and Sonenberg, N. (2016).
#' Translational control by 5'-untranslated regions of eukaryotic
#' mRNAs. \emph{Science} 352(6292), 1413-1416.
#' \doi{10.1126/science.aad9868}
#'
#' Mayr, C. (2019). What are 3' UTRs doing? \emph{Cold Spring Harbor
#' Perspectives in Biology} 11(10), a034728.
#' \doi{10.1101/cshperspect.a034728}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # CONSTRUCT_03 carries three upstream AUGs, two of which open
#' # uORFs, and a premature stop codon in its coding sequence.
#' features <- profileConstruct(
#'   constructs = simulatedConstructs[c(3, 12), ]
#' )
#' map <- plotConstructMap(
#'   feature_table = features, construct_id = "CONSTRUCT_03"
#' )
#' map
#'
#' # The result is a patchwork object, so it can be customised.
#' map + patchwork::plot_annotation(
#'   title = "Candidate CONSTRUCT_03"
#' )
#'
#' # CONSTRUCT_12 carries nine AU-rich element pentamers instead.
#' plotConstructMap(
#'   feature_table = features, construct_id = "CONSTRUCT_12"
#' )
#'
#' @importFrom ggplot2 .data
#' @export
plotConstructMap <- function(feature_table, construct_id) {
  validateFeatureTable(feature_table = feature_table)
  construct_id <- validateConstructChoice(
    feature_table = feature_table, construct_id = construct_id
  )

  rows <- feature_table[
    feature_table$construct_id == construct_id, , drop = FALSE
  ]
  offsets <- regionOffsets(rows = rows)

  region_track <- buildRegionTrack(offsets = offsets)
  utr5_track <- buildUtr5Track(rows = rows, offsets = offsets)
  structure_track <- buildStructureTrack(
    rows = rows, offsets = offsets
  )
  cds_track <- buildCdsTrack(rows = rows, offsets = offsets)
  utr3_track <- buildUtr3Track(rows = rows, offsets = offsets)

  # Only the bottom track carries the axis: the tracks share one
  # coordinate system, so repeating its labels five times would add
  # ink without adding information.
  upper_tracks <- lapply(
    X = list(region_track, utr5_track, structure_track, cds_track),
    FUN = function(track) track + hideXAxis()
  )

  assembled <- patchwork::wrap_plots(
    c(upper_tracks, list(utr3_track)),
    ncol = 1L,
    heights = c(0.4, 1.1, 1.0, 0.8, 1.2)
  ) +
    patchwork::plot_layout(guides = "collect")

  assembled + patchwork::plot_annotation(
    title = paste("Construct map:", construct_id),
    subtitle = paste(
      "Tracks share one coordinate axis, shown at the foot.",
      "Coloured tracks are scored; the CDS track is context."
    ),
    caption = paste(
      "The CDS track is grey because coding sequence features are",
      "reported as context and take no part in any score or\n",
      "ranking. The structure track is a sequence-only proxy, not a",
      "folding prediction. Every feature carries a shape and a\n",
      "label as well as a colour, so no reading depends on telling",
      "colours apart."
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
  if (is.na(x = construct_id) ||
        !construct_id %in% available) {
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
#' Positions in the feature table are relative to their own region, so
#' drawing them on one axis means offsetting each region by the lengths
#' of those before it.
#'
#' @param rows The feature rows for one construct.
#' @return A data frame of region, length, start and end.
#' @noRd
regionOffsets <- function(rows) {
  regions <- c("utr5", "cds", "utr3")
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

#' Select the positional rows for one feature family
#'
#' @param rows The feature rows for one construct.
#' @param region Region to take rows from.
#' @param pattern Regular expression matching the feature names.
#' @param offsets The region offset table.
#' @return A data frame with an added `x` column, possibly empty.
#' @noRd
featureMarks <- function(rows, region, pattern, offsets) {
  selected <- rows[
    rows$region == region &
      !is.na(rows$position) &
      grepl(pattern = pattern, x = rows$feature),
    ,
    drop = FALSE
  ]
  if (nrow(x = selected) == 0L) {
    return(selected)
  }
  selected$x <- absolutePosition(rows = selected, offsets = offsets)
  selected
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

#' The 5' UTR track: upstream starts, uORF extents and motifs
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @return A ggplot object.
#' @noRd
buildUtr5Track <- function(rows, offsets) {
  palette <- utrunwayPalette()
  labels <- utrunwayLabels()
  shapes <- utrunwayShapes()

  augs <- featureMarks(
    rows = rows, region = "utr5", pattern = "^upstream_aug_",
    offsets = offsets
  )
  orfs <- featureMarks(
    rows = rows, region = "utr5", pattern = "^upstream_orf$",
    offsets = offsets
  )
  motifs <- featureMarks(
    rows = rows, region = "utr5", pattern = "^motif_",
    offsets = offsets
  )

  plot <- ggplot2::ggplot() + sharedXScale(offsets = offsets)

  # uORF extents first, as segments beneath the start codons that
  # open them, so a reader sees how far each one reads.
  if (nrow(x = orfs) > 0L) {
    orfs$xend <- orfs$x + orfs$value - 1
    plot <- plot +
      ggplot2::geom_segment(
        data = orfs,
        mapping = ggplot2::aes(
          x = .data$x, xend = .data$xend, y = 0.45, yend = 0.45
        ),
        colour = palette[["upstream_aug_adequate"]],
        linewidth = 1.6,
        lineend = "round"
      )
    # Lengths are labelled only when the uORFs are far enough apart
    # for the labels not to collide; otherwise the segments speak for
    # themselves and the counts are in the feature table.
    orfs$midpoint <- (orfs$x + orfs$xend) / 2
    spacing <- if (nrow(x = orfs) > 1L) {
      min(diff(sort(x = orfs$midpoint)))
    } else {
      Inf
    }
    if (spacing > max(offsets$end) * 0.06) {
      plot <- plot +
        ggplot2::geom_text(
          data = orfs,
          mapping = ggplot2::aes(
            x = .data$midpoint, y = 0.15,
            label = paste0(.data$value, " nt")
          ),
          colour = palette[["ink_secondary"]], size = 2.6
        )
    }
  }

  if (nrow(x = augs) > 0L) {
    augs$strength <- sub(
      pattern = "^upstream_aug_", replacement = "", x = augs$feature
    )
    augs$class <- labels[augs$feature]
    plot <- plot +
      ggplot2::geom_point(
        data = augs,
        mapping = ggplot2::aes(
          x = .data$x, y = 1, fill = .data$class
        ),
        shape = shapes[["upstream_aug"]],
        colour = palette[["ink"]],
        size = 3,
        stroke = 0.4
      )
  }

  if (nrow(x = motifs) > 0L) {
    motifs$class <- sub(
      pattern = "^motif_", replacement = "", x = motifs$feature
    )
    motifs$class <- ifelse(
      test = motifs$class == "TOP",
      yes = labels[["top_motif"]],
      no = labels[["g_quadruplex"]]
    )
    motifs$shape <- ifelse(
      test = grepl(pattern = "TOP", x = motifs$feature),
      yes = shapes[["top_motif"]],
      no = shapes[["g_quadruplex"]]
    )
    motifs$xend <- motifs$x + motifs$value - 1
    plot <- plot +
      ggplot2::geom_segment(
        data = motifs,
        mapping = ggplot2::aes(
          x = .data$x, xend = .data$xend, y = 1.55, yend = 1.55,
          colour = .data$class
        ),
        linewidth = 1.6,
        lineend = "round",
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = motifs,
        mapping = ggplot2::aes(
          x = (.data$x + .data$xend) / 2, y = 1.55,
          fill = .data$class, shape = .data$class
        ),
        colour = palette[["ink"]],
        size = 2.6,
        stroke = 0.4
      )
  }

  # A manual scale is only added when a layer actually maps to it;
  # adding one with nothing mapped warns about unshared levels.
  if (nrow(x = augs) > 0L || nrow(x = motifs) > 0L) {
    plot <- plot +
      ggplot2::scale_fill_manual(
        values = featureFillValues(), drop = TRUE,
        na.translate = FALSE
      )
  }
  if (nrow(x = motifs) > 0L) {
    plot <- plot +
      ggplot2::scale_colour_manual(
        values = featureFillValues(), drop = TRUE, guide = "none"
      ) +
      ggplot2::scale_shape_manual(
        values = featureShapeValues(), drop = TRUE
      )
  }

  plot +
    ggplot2::scale_y_continuous(limits = c(0, 2), breaks = NULL) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "5' UTR: initiation features",
      subtitle = if (nrow(x = augs) == 0L) {
        "No upstream AUG, and no TOP or G-quadruplex motif"
      } else {
        NULL
      }
    ) +
    utrunwayTheme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

#' Colours keyed by the labels the legend shows
#'
#' @return A named character vector of hex colours.
#' @noRd
featureFillValues <- function() {
  palette <- utrunwayPalette()
  labels <- utrunwayLabels()

  stats::setNames(
    object = c(
      palette[["upstream_aug_weak"]],
      palette[["upstream_aug_adequate"]],
      palette[["upstream_aug_strong"]],
      palette[["upstream_aug_unscored"]],
      palette[["top_motif"]],
      palette[["g_quadruplex"]],
      palette[["polya_signal"]],
      palette[["polya_signal"]],
      palette[["are_pentamer"]],
      palette[["mirna_seed"]]
    ),
    nm = c(
      labels[["upstream_aug_weak"]],
      labels[["upstream_aug_adequate"]],
      labels[["upstream_aug_strong"]],
      labels[["upstream_aug_unscored"]],
      labels[["top_motif"]],
      labels[["g_quadruplex"]],
      labels[["polya_signal_canonical"]],
      labels[["polya_signal_variant"]],
      labels[["are_pentamer"]],
      labels[["mirna_seed"]]
    )
  )
}

#' Shapes keyed by the labels the legend shows
#'
#' @return A named integer vector of pch values.
#' @noRd
featureShapeValues <- function() {
  labels <- utrunwayLabels()
  shapes <- utrunwayShapes()

  stats::setNames(
    object = c(
      rep(x = shapes[["upstream_aug"]], times = 4L),
      shapes[["top_motif"]],
      shapes[["g_quadruplex"]],
      shapes[["polya_signal_canonical"]],
      shapes[["polya_signal_variant"]],
      shapes[["are_pentamer"]],
      shapes[["mirna_seed"]]
    ),
    nm = c(
      labels[["upstream_aug_weak"]],
      labels[["upstream_aug_adequate"]],
      labels[["upstream_aug_strong"]],
      labels[["upstream_aug_unscored"]],
      labels[["top_motif"]],
      labels[["g_quadruplex"]],
      labels[["polya_signal_canonical"]],
      labels[["polya_signal_variant"]],
      labels[["are_pentamer"]],
      labels[["mirna_seed"]]
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
      x = NULL, y = "GC",
      title = "Structure proxy: sliding-window GC content",
      subtitle = paste(
        "Open triangles mark windows holding a complementary run of",
        "8 nt or more. A proxy, not a folding prediction."
      )
    ) +
    utrunwayTheme()
}

#' The coding sequence track, drawn as context
#'
#' Everything here is grey, because none of it is scored.
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @return A ggplot object.
#' @noRd
buildCdsTrack <- function(rows, offsets) {
  palette <- utrunwayPalette()

  gc3 <- rows$value[rows$region == "cds" & rows$feature == "gc3_content"]
  gc_overall <- rows$value[
    rows$region == "cds" & rows$feature == "gc_content"
  ]
  stops <- featureMarks(
    rows = rows, region = "cds", pattern = "^premature_stop$",
    offsets = offsets
  )
  cryptic <- featureMarks(
    rows = rows, region = "cds", pattern = "^out_of_frame_atg$",
    offsets = offsets
  )
  cds_row <- offsets[offsets$region == "cds", ]

  plot <- ggplot2::ggplot() +
    sharedXScale(offsets = offsets) +
    ggplot2::geom_rect(
      data = cds_row,
      mapping = ggplot2::aes(
        xmin = .data$start + 1, xmax = .data$end - 1,
        ymin = 0.35, ymax = 0.65
      ),
      fill = palette[["context_fill"]]
    )

  if (nrow(x = cryptic) > 0L) {
    plot <- plot +
      ggplot2::geom_point(
        data = cryptic,
        mapping = ggplot2::aes(x = .data$x, y = 0.5),
        shape = 124L,
        colour = palette[["context"]],
        size = 2.2
      )
  }
  if (nrow(x = stops) > 0L) {
    plot <- plot +
      ggplot2::geom_point(
        data = stops,
        mapping = ggplot2::aes(x = .data$x, y = 0.5),
        shape = 4L,
        colour = palette[["ink"]],
        size = 2.6,
        stroke = 0.7
      )
  }

  plot +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "CDS: context only, never scored",
      subtitle = paste0(
        "GC ", round(x = gc_overall * 100), "%, GC3 ",
        round(x = gc3 * 100), "%. Vertical ticks mark out-of-frame ",
        "ATGs; crosses mark premature in-frame stops."
      )
    ) +
    utrunwayTheme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

#' The 3' UTR track: polyA signals, AU-rich elements and seed sites
#'
#' @param rows The feature rows for one construct.
#' @param offsets The region offset table.
#' @return A ggplot object.
#' @noRd
buildUtr3Track <- function(rows, offsets) {
  palette <- utrunwayPalette()
  labels <- utrunwayLabels()

  elements <- featureMarks(
    rows = rows,
    region = "utr3",
    pattern = "^(polya_signal_|are_pentamer$|mirna_seed_)",
    offsets = offsets
  )

  if (nrow(x = elements) == 0L) {
    return(emptyTrack(
      offsets = offsets,
      title = "3' UTR: no decay or repression elements found"
    ))
  }

  elements$class <- ifelse(
    test = elements$feature == "polya_signal_canonical",
    yes = labels[["polya_signal_canonical"]],
    no = ifelse(
      test = elements$feature == "polya_signal_variant",
      yes = labels[["polya_signal_variant"]],
      no = ifelse(
        test = elements$feature == "are_pentamer",
        yes = labels[["are_pentamer"]],
        no = labels[["mirna_seed"]]
      )
    )
  )
  # One row per class, so overlapping pentamers do not pile up on a
  # single line and become unreadable.
  class_levels <- c(
    labels[["polya_signal_canonical"]],
    labels[["polya_signal_variant"]],
    labels[["are_pentamer"]],
    labels[["mirna_seed"]]
  )
  class_levels <- class_levels[class_levels %in% elements$class]
  elements$y <- match(x = elements$class, table = class_levels)

  ggplot2::ggplot() +
    sharedXScale(offsets = offsets) +
    # Each class has its own row, named on the axis. That direct
    # label is the identity channel here, so the legend would only
    # repeat it.
    ggplot2::geom_point(
      data = elements,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y,
        fill = .data$class, shape = .data$class
      ),
      colour = palette[["ink"]],
      size = 2.8,
      stroke = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = featureFillValues(), drop = TRUE, na.translate = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = featureShapeValues(), drop = TRUE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, length(x = class_levels) + 0.5),
      breaks = seq_along(along.with = class_levels),
      labels = class_levels
    ) +
    ggplot2::labs(
      x = "Position in construct (nt)", y = NULL,
      title = "3' UTR: decay and repression elements"
    ) +
    utrunwayTheme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        colour = palette[["ink_secondary"]], size = 7
      )
    )
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
