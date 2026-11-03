# Build the package overview image, inst/extdata/MO_M_A1.png.
#
# The diagram shows what goes into UTRunway and what comes out: the
# three regions of a candidate construct entering readConstruct(),
# flowing through profileConstruct() into the tidy feature table, and
# branching from there to the ranking and to the figures.
#
# The coding sequence path is drawn differently on purpose. It is
# grey, its border and its arrow are dashed, and it is labelled as
# context, because coding sequence features are profiled and plotted
# but never scored. The package's central scope rule is therefore
# visible in the diagram itself rather than only in its caption.
#
# Colours are the package's own palette, so a reader who has seen a
# UTRunway figure recognises them. Nothing here is carried by colour
# alone: every box is labelled and the two paths differ in line type
# as well as in shade.
#
# Run with:  source("data-raw/makeOverviewImage.R")

library(ggplot2)

# The package palette, repeated here so the script stands alone.
ink <- "#0b0b0b"
ink_soft <- "#52514e"
muted <- "#898781"
surface <- "#fcfcfb"
rule <- "#c3c2b7"
step_fill <- "#cde2fb"
table_fill <- "#e8f1fd"
utr_fill <- "#d7e9d7"
context_fill <- "#e4e3df"
output_fill <- "#fbe6cd"

boxes <- data.frame(
  key = c(
    "utr5", "cds", "utr3", "read", "profile", "table", "rank",
    "plots", "out_rank", "out_map"
  ),
  x = c(11, 11, 11, 33, 55, 79, 105, 105, 129, 129),
  y = c(40, 29, 18, 29, 29, 29, 42, 15, 42, 15),
  width = c(16, 16, 16, 18, 18, 22, 22, 24, 20, 20),
  height = c(7.5, 7.5, 7.5, 10, 10, 14, 10, 16, 10, 10),
  fill = c(
    utr_fill, context_fill, utr_fill, step_fill, step_fill,
    table_fill, step_fill, step_fill, output_fill, output_fill
  ),
  border = c(rule, muted, rule, muted, muted, muted, muted, muted,
             muted, muted),
  linetype = c(
    "solid", "dashed", "solid", "solid", "solid", "solid", "solid",
    "solid", "solid", "solid"
  ),
  label = c(
    "5' UTR",
    "CDS",
    "3' UTR",
    "readConstruct()",
    "profileConstruct()",
    "Tidy feature table",
    "rankCandidates()\nrankStability()",
    paste(
      "plotConstructMap()",
      "plotRegionProfile()",
      "plotCandidateComparison()",
      "plotScoreBreakdown()",
      "plotCandidateRanking()",
      sep = "\n"
    ),
    "Ranked comparison",
    "Construct map"
  ),
  stringsAsFactors = FALSE
)

# Second lines, set smaller and softer than the box titles.
subtitles <- data.frame(
  x = c(11, 33, 55, 79, 129, 129),
  y = c(26.4, 25.6, 25.6, 25.4, 38.6, 11.6),
  label = c(
    "context only",
    "read and validate",
    "measure features",
    paste(
      "construct_id | region |",
      "feature | value | position",
      sep = "\n"
    ),
    "scores, contributions,\nweights, rank intervals",
    "every feature on one axis"
  ),
  stringsAsFactors = FALSE
)

# Solid arrows carry sequence that will be scored; the dashed one
# carries the coding sequence, which will not be.
arrows <- data.frame(
  x = c(19, 19, 19, 42, 64, 90, 90, 115, 115),
  y = c(40, 29, 18, 29, 29, 31, 27, 42, 15),
  xend = c(23.6, 23.6, 23.6, 45.6, 67.6, 93.4, 92.6, 118.6, 118.6),
  yend = c(32.5, 29, 25.5, 29, 29, 39, 20, 42, 15),
  linetype = c(
    "solid", "dashed", "solid", "solid", "solid", "solid", "solid",
    "solid", "solid"
  ),
  colour = c(ink_soft, muted, ink_soft, ink_soft, ink_soft, ink_soft,
             ink_soft, ink_soft, ink_soft),
  stringsAsFactors = FALSE
)

overview <- ggplot(boxes) +
  geom_segment(
    data = arrows,
    mapping = aes(
      x = x, xend = xend, y = y, yend = yend,
      linetype = I(linetype), colour = I(colour)
    ),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.5
  ) +
  geom_rect(
    mapping = aes(
      xmin = x - width / 2, xmax = x + width / 2,
      ymin = y - height / 2, ymax = y + height / 2,
      fill = I(fill), colour = I(border), linetype = I(linetype)
    ),
    linewidth = 0.4
  ) +
  geom_text(
    mapping = aes(x = x, y = y, label = label),
    size = 3.1, colour = ink, lineheight = 1.25, fontface = "plain"
  ) +
  geom_text(
    data = subtitles,
    mapping = aes(x = x, y = y, label = label),
    size = 2.5, colour = ink_soft, lineheight = 1.2
  ) +
  # Column headings, so the direction of travel is named.
  annotate(
    geom = "text",
    x = c(11, 44, 79, 105, 129),
    y = rep(x = 50.5, times = 5L),
    label = c(
      "Input: one candidate", "Read and profile", "One table",
      "Rank and draw", "Output"
    ),
    size = 2.9, colour = ink, fontface = "bold"
  ) +
  annotate(
    geom = "segment",
    x = 3, xend = 139, y = 48.2, yend = 48.2,
    colour = rule, linewidth = 0.3
  ) +
  # The scope rule, stated where the diagram can be read without it.
  annotate(
    geom = "segment",
    x = 46, xend = 52, y = 5, yend = 5,
    colour = muted, linewidth = 0.5, linetype = "dashed"
  ) +
  annotate(
    geom = "text", x = 53.5, y = 5, hjust = 0,
    label = paste(
      "Dashed and grey: the coding sequence is profiled and plotted",
      "as context, and is never scored."
    ),
    size = 2.7, colour = ink_soft
  ) +
  annotate(
    geom = "text", x = 3, y = 1.2, hjust = 0,
    label = paste(
      "UTRunway: profile and rank candidate mRNA construct designs",
      "from sequence alone."
    ),
    size = 2.7, colour = muted
  ) +
  scale_x_continuous(limits = c(0, 142)) +
  scale_y_continuous(limits = c(0, 53)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = surface, colour = NA),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
  )

destination <- file.path("inst", "extdata", "MO_M_A1.png")
ggsave(
  filename = destination,
  plot = overview,
  width = 12.5,
  height = 4.9,
  dpi = 200,
  bg = surface
)

message("Wrote ", destination, " (",
        round(file.size(destination) / 1024), " KB)")
