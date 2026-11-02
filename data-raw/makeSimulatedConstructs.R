# Build the bundled `simulatedConstructs` dataset.
#
# Twelve synthetic mRNA construct designs with deliberately planted
# sequence features. The planted features are recorded alongside the
# sequences so that tests can check detection against a known ground
# truth. These are simulated sequences only: they carry no measured
# biological outcomes and are not derived from any real transcript.
#
# Run with:  source("data-raw/makeSimulatedConstructs.R")

set.seed(seed = 23)

bases <- c("A", "C", "G", "T")

stop_codons <- c("TAA", "TAG", "TGA")

allCodons <- function() {
  as.vector(outer(
    X = outer(X = bases, Y = bases, FUN = paste0),
    Y = bases,
    FUN = paste0
  ))
}

sense_codons <- setdiff(x = allCodons(), y = stop_codons)

# Random sequence of a given length, uniform over the four bases.
randomSequence <- function(len) {
  paste(sample(x = bases, size = len, replace = TRUE), collapse = "")
}

# Count occurrences of `motif`, counting overlapping hits separately.
# AUUUA pentamers in an AU-rich element overlap ("ATTTATTTA" holds
# two), so gregexpr(), which skips past each match, would undercount.
countMotif <- function(sequence, motif) {
  n_starts <- nchar(x = sequence) - nchar(x = motif) + 1L
  if (n_starts < 1L) {
    return(0L)
  }
  starts <- seq_len(length.out = n_starts)
  sum(
    substring(
      text = sequence,
      first = starts,
      last = starts + nchar(x = motif) - 1L
    ) == motif
  )
}

# Remove every occurrence of `motif` from `sequence` by mutating one
# base inside each hit, so that unplanted background cannot imitate a
# planted feature. Loops because a mutation can reveal a new overlap.
maskMotif <- function(sequence, motif, at, replacement) {
  repeat {
    hit <- regexpr(pattern = motif, text = sequence, fixed = TRUE)
    if (hit[[1]] == -1L) {
      return(sequence)
    }
    pos <- hit[[1]] + at - 1L
    substr(x = sequence, start = pos, stop = pos) <- replacement
  }
}

# Count complete upstream ORFs in a 5' UTR: AUGs whose in-frame stop
# codon also lies inside the UTR. An AUG with no in-frame stop before
# the UTR ends overlaps the CDS instead and is not counted here.
countUpstreamOrfs <- function(utr5) {
  len <- nchar(x = utr5)
  n_orf <- 0L
  aug_starts <- which(
    substring(
      text = utr5,
      first = seq_len(length.out = max(len - 2L, 0L)),
      last = seq_len(length.out = max(len - 2L, 0L)) + 2L
    ) == "ATG"
  )
  for (start_i in aug_starts) {
    codon_starts <- seq(from = start_i, to = len - 2L, by = 3L)
    codons <- substring(
      text = utr5, first = codon_starts, last = codon_starts + 2L
    )
    if (any(codons %in% stop_codons)) {
      n_orf <- n_orf + 1L
    }
  }
  n_orf
}

# The G-quadruplex-prone arrangement: four G-triplets separated by
# short loops. This is the same pattern the package's tests use, and
# the two must stay in step: screening the background for a bare
# "GGGG" run would not catch an arrangement such as "GGGAGGGTGGGCGGG",
# which contains four G-triplets but no four G's in a row.
g_quad_pattern <- paste0(
  "G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"
)

# Remove every match of a regular expression from `sequence` by
# mutating the first base of each match. Loops because shortening one
# G-run can leave a shorter arrangement that still matches.
maskPattern <- function(sequence, pattern, replacement = "A") {
  repeat {
    hit <- regexpr(pattern = pattern, text = sequence)
    if (hit[[1]] == -1L) {
      return(sequence)
    }
    pos <- hit[[1]]
    substr(x = sequence, start = pos, stop = pos) <- replacement
  }
}

# Count matches of a regular expression in `sequence`.
countPattern <- function(sequence, pattern) {
  hits <- gregexpr(pattern = pattern, text = sequence)[[1]]
  if (hits[[1]] == -1L) 0L else length(x = hits)
}

# Background for a 5' UTR: no ATG (so upstream AUG counts are exact),
# no G-quadruplex-prone arrangement and no four-G run (so planted
# structure is the only structure), and a purine first base (so no
# accidental TOP motif at the 5' end).
makeUtr5Background <- function(len) {
  seq_out <- randomSequence(len = len)
  # Screen G-runs before ATG: these masks turn a G into an A, which
  # can spell a new ATG out of a GTG, whereas the ATG mask writes a C
  # and so cannot create a G-run.
  seq_out <- maskPattern(
    sequence = seq_out, pattern = g_quad_pattern
  )
  seq_out <- maskMotif(
    sequence = seq_out, motif = "GGGG", at = 2L, replacement = "C"
  )
  seq_out <- maskMotif(
    sequence = seq_out, motif = "ATG", at = 2L, replacement = "C"
  )
  substr(x = seq_out, start = 1L, stop = 1L) <-
    sample(x = c("A", "G"), size = 1L)
  seq_out
}

# Background for a 3' UTR: no AUUUA pentamer and no AATAAA hexamer, so
# that planted AU-rich elements and polyA signals are the only ones.
makeUtr3Background <- function(len) {
  seq_out <- randomSequence(len = len)
  seq_out <- maskMotif(
    sequence = seq_out, motif = "ATTTA", at = 3L, replacement = "C"
  )
  seq_out <- maskMotif(
    sequence = seq_out, motif = "AATAAA", at = 3L, replacement = "C"
  )
  seq_out
}

# A coding sequence in frame: ATG, then sense codons only, then a stop.
# `second_codon_start` fixes the +4 base, which is part of the Kozak
# context. `premature_stop_codon_index` plants one in-frame stop codon.
makeCds <- function(n_codons, second_codon_start,
                    premature_stop_codon_index = NA_integer_) {
  body <- sample(x = sense_codons, size = n_codons - 2L, replace = TRUE)
  body[[1]] <- sample(
    x = grep(
      pattern = paste0("^", second_codon_start),
      x = sense_codons,
      value = TRUE
    ),
    size = 1L
  )
  if (!is.na(premature_stop_codon_index)) {
    body[[premature_stop_codon_index]] <- sample(x = stop_codons, size = 1L)
  }
  paste(c("ATG", body, "TAA"), collapse = "")
}

# An in-frame open reading frame to plant inside a 5' UTR: AUG, a few
# sense codons, then a stop codon, all contained within the UTR.
makeUpstreamOrf <- function(n_codons) {
  paste(
    c(
      "ATG",
      sample(x = sense_codons, size = n_codons, replace = TRUE),
      sample(x = stop_codons, size = 1L)
    ),
    collapse = ""
  )
}

# Splice `insert` into `sequence` at `pos`, keeping the total length.
spliceIn <- function(sequence, insert, pos) {
  paste0(
    substr(x = sequence, start = 1L, stop = pos - 1L),
    insert,
    substr(
      x = sequence,
      start = pos + nchar(x = insert),
      stop = nchar(x = sequence)
    )
  )
}

kozak_strong <- "GCCACC"
kozak_weak <- "TTTTTT"
# Four G-triplets, the arrangement scored as G-quadruplex-prone.
g_quad_run <- "GGGAGGGAGGGAGGGA"
# Pyrimidine-rich 5' terminal oligopyrimidine motif.
top_motif <- "CTTTTCTTTCT"
# Overlapping AUUUA pentamers in an AU context.
are_run <- "TATTTATTTATTTATT"
are_run_pentamers <- 3L
polya_canonical <- "AATAAA"
polya_disrupted <- "AAGAAA"

# Assemble one 5' UTR carrying its planted features. Returns NULL if
# splicing happened to create an extra ATG or G-run at a junction, in
# which case the caller retries with a fresh background.
buildUtr5 <- function(len, n_upstream_aug, n_upstream_orf, kozak,
                      has_g_quadruplex, has_top_motif) {
  seq_out <- makeUtr5Background(len = len)

  if (has_top_motif) {
    seq_out <- spliceIn(sequence = seq_out, insert = top_motif, pos = 1L)
  }
  if (has_g_quadruplex) {
    seq_out <- spliceIn(
      sequence = seq_out, insert = g_quad_run, pos = round(x = len * 0.3)
    )
  }

  # Planted uORFs sit in the middle of the UTR, each wholly inside it.
  if (n_upstream_orf > 0L) {
    for (i in seq_len(length.out = n_upstream_orf)) {
      orf <- makeUpstreamOrf(n_codons = 4L + i)
      seq_out <- spliceIn(
        sequence = seq_out,
        insert = orf,
        pos = round(x = len * (0.45 + 0.12 * (i - 1L)))
      )
    }
  }

  # Remaining upstream AUGs are bare starts with no in-frame stop
  # before the UTR ends, placed just before the Kozak block.
  n_bare_aug <- n_upstream_aug - n_upstream_orf
  if (n_bare_aug > 0L) {
    for (i in seq_len(length.out = n_bare_aug)) {
      seq_out <- spliceIn(
        sequence = seq_out,
        insert = "ATG",
        pos = len - nchar(x = kozak_strong) - 4L * i
      )
    }
  }

  # The Kozak context occupies the bases immediately before the CDS.
  if (!is.na(kozak)) {
    seq_out <- spliceIn(
      sequence = seq_out,
      insert = kozak,
      pos = len - nchar(x = kozak) + 1L
    )
  }

  if (countMotif(sequence = seq_out, motif = "ATG") != n_upstream_aug) {
    return(NULL)
  }

  # A bare upstream AUG must not pick up a chance in-frame stop inside
  # the UTR, which would quietly make it a complete uORF.
  if (countUpstreamOrfs(utr5 = seq_out) != n_upstream_orf) {
    return(NULL)
  }

  # A four-G run outside the planted G-quadruplex block would blur the
  # structure ground truth.
  if (countMotif(sequence = seq_out, motif = "GGGG") > 0L) {
    return(NULL)
  }

  # The quadruplex-prone arrangement must appear exactly where it was
  # planted and nowhere else, including across splice junctions.
  n_quad <- countPattern(sequence = seq_out, pattern = g_quad_pattern)
  if (n_quad != (if (has_g_quadruplex) 1L else 0L)) {
    return(NULL)
  }

  seq_out
}

# Assemble one 3' UTR carrying its planted features. Like buildUtr5(),
# returns NULL when a splice junction created an unplanned motif.
buildUtr3 <- function(len, n_are_pentamer, is_polya_disrupted) {
  seq_out <- makeUtr3Background(len = len)

  # Each run contributes three overlapping AUUUA pentamers.
  if (n_are_pentamer > 0L) {
    n_runs <- n_are_pentamer / are_run_pentamers
    for (i in seq_len(length.out = n_runs)) {
      seq_out <- spliceIn(
        sequence = seq_out,
        insert = are_run,
        pos = round(x = len * (0.25 + 0.2 * (i - 1L)))
      )
    }
  }

  signal <- if (is_polya_disrupted) polya_disrupted else polya_canonical
  seq_out <- spliceIn(
    sequence = seq_out,
    insert = signal,
    pos = len - nchar(x = signal) - 19L
  )

  # A run abutting the background can spell one extra pentamer, and a
  # junction can spell an unplanned polyA hexamer.
  if (countMotif(sequence = seq_out, motif = "ATTTA") != n_are_pentamer) {
    return(NULL)
  }
  n_polya <- countMotif(sequence = seq_out, motif = polya_canonical)
  if (n_polya != (if (is_polya_disrupted) 0L else 1L)) {
    return(NULL)
  }

  seq_out
}

# Retry wrapper: backgrounds are random, so a junction artefact just
# means drawing another one.
retryBuild <- function(builder, ..., max_tries = 500L) {
  for (try_i in seq_len(length.out = max_tries)) {
    candidate <- builder(...)
    if (!is.null(x = candidate)) {
      return(candidate)
    }
  }
  stop("Could not build a region with the requested planted features.")
}

# ---- Planted-feature plan -------------------------------------------
#
# One row per construct. Features are spread so that most constructs
# carry more than one, and so that each feature has both positive and
# negative examples.

plan <- data.frame(
  construct_id = sprintf(fmt = "CONSTRUCT_%02d", 1:12),
  utr5_length = c(
    120L, 160L, 200L, 95L, 140L, 180L, 110L, 150L, 190L, 85L, 130L, 170L
  ),
  cds_codons = c(
    120L, 150L, 180L, 110L, 140L, 170L, 105L, 135L, 165L, 100L, 125L, 195L
  ),
  utr3_length = c(
    260L, 320L, 400L, 210L, 300L, 380L, 240L, 340L, 460L, 200L, 280L, 500L
  ),
  n_upstream_aug = c(
    1L, 2L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L
  ),
  n_upstream_orf = c(
    1L, 1L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L
  ),
  kozak_class = c(
    "moderate", "moderate", "moderate", "strong", "strong", "strong",
    "weak", "weak", "weak", "moderate", "moderate", "moderate"
  ),
  has_g_quadruplex = c(
    FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
    TRUE, FALSE, FALSE
  ),
  has_top_motif = c(
    FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE
  ),
  n_are_pentamer = c(
    3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 9L
  ),
  is_polya_disrupted = c(
    FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE,
    FALSE, FALSE, FALSE
  ),
  has_premature_stop = c(
    FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
    FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

kozak_insert <- c(strong = kozak_strong, weak = kozak_weak)
second_base <- c(strong = "G", weak = "C", moderate = "A")

utr5_vec <- character(length = nrow(x = plan))
cds_vec <- character(length = nrow(x = plan))
utr3_vec <- character(length = nrow(x = plan))
premature_stop_codon <- rep(x = NA_integer_, times = nrow(x = plan))

for (i in seq_len(length.out = nrow(x = plan))) {
  row_i <- plan[i, ]

  kozak_i <- if (row_i$kozak_class == "moderate") {
    NA_character_
  } else {
    kozak_insert[[row_i$kozak_class]]
  }

  utr5_vec[[i]] <- retryBuild(
    builder = buildUtr5,
    len = row_i$utr5_length,
    n_upstream_aug = row_i$n_upstream_aug,
    n_upstream_orf = row_i$n_upstream_orf,
    kozak = kozak_i,
    has_g_quadruplex = row_i$has_g_quadruplex,
    has_top_motif = row_i$has_top_motif
  )

  if (row_i$has_premature_stop) {
    premature_stop_codon[[i]] <- round(x = row_i$cds_codons * 0.4)
  }

  cds_vec[[i]] <- makeCds(
    n_codons = row_i$cds_codons,
    second_codon_start = second_base[[row_i$kozak_class]],
    premature_stop_codon_index = premature_stop_codon[[i]]
  )

  utr3_vec[[i]] <- retryBuild(
    builder = buildUtr3,
    len = row_i$utr3_length,
    n_are_pentamer = row_i$n_are_pentamer,
    is_polya_disrupted = row_i$is_polya_disrupted
  )
}

simulatedConstructs <- data.frame(
  construct_id = plan$construct_id,
  utr5 = utr5_vec,
  cds = cds_vec,
  utr3 = utr3_vec,
  n_upstream_aug = plan$n_upstream_aug,
  n_upstream_orf = plan$n_upstream_orf,
  has_strong_kozak = plan$kozak_class == "strong",
  has_weak_kozak = plan$kozak_class == "weak",
  has_g_quadruplex = plan$has_g_quadruplex,
  has_top_motif = plan$has_top_motif,
  n_are_pentamer = plan$n_are_pentamer,
  has_canonical_polya = !plan$is_polya_disrupted,
  has_premature_stop = plan$has_premature_stop,
  stringsAsFactors = FALSE
)

# ---- Check the ground truth before saving ---------------------------
#
# Every planted feature is counted back out of the finished sequences,
# so the stored truth columns cannot drift from the sequences.

inFrameStopCount <- function(cds) {
  codons <- substring(
    text = cds,
    first = seq(from = 1L, to = nchar(x = cds) - 2L, by = 3L),
    last = seq(from = 3L, to = nchar(x = cds), by = 3L)
  )
  sum(codons[-length(x = codons)] %in% stop_codons)
}

for (i in seq_len(length.out = nrow(x = simulatedConstructs))) {
  row_i <- simulatedConstructs[i, ]

  stopifnot(
    "5' UTR length is outside 80-200 nt" =
      nchar(x = row_i$utr5) >= 80L && nchar(x = row_i$utr5) <= 200L,
    "3' UTR length is outside 200-500 nt" =
      nchar(x = row_i$utr3) >= 200L && nchar(x = row_i$utr3) <= 500L,
    "CDS length is outside 300-600 nt" =
      nchar(x = row_i$cds) >= 300L && nchar(x = row_i$cds) <= 600L,
    "CDS length is not a multiple of three" =
      nchar(x = row_i$cds) %% 3L == 0L,
    "CDS does not start with ATG" =
      substr(x = row_i$cds, start = 1L, stop = 3L) == "ATG",
    "CDS does not end with a stop codon" =
      substr(
        x = row_i$cds,
        start = nchar(x = row_i$cds) - 2L,
        stop = nchar(x = row_i$cds)
      ) %in% stop_codons,
    "upstream AUG count does not match the planted truth" =
      countMotif(sequence = row_i$utr5, motif = "ATG") ==
        row_i$n_upstream_aug,
    "upstream ORF count does not match the planted truth" =
      countUpstreamOrfs(utr5 = row_i$utr5) == row_i$n_upstream_orf,
    "AUUUA pentamer count does not match the planted truth" =
      countMotif(sequence = row_i$utr3, motif = "ATTTA") ==
        row_i$n_are_pentamer,
    "G-quadruplex-prone block does not match the planted truth" =
      (countMotif(sequence = row_i$utr5, motif = g_quad_run) > 0L) ==
        row_i$has_g_quadruplex,
    "G-quadruplex arrangement does not match the planted truth" =
      (countPattern(
        sequence = row_i$utr5, pattern = g_quad_pattern
      ) > 0L) == row_i$has_g_quadruplex,
    "polyA signal does not match the planted truth" =
      (countMotif(sequence = row_i$utr3, motif = polya_canonical) > 0L) ==
        row_i$has_canonical_polya,
    "premature in-frame stop does not match the planted truth" =
      (inFrameStopCount(cds = row_i$cds) > 0L) == row_i$has_premature_stop,
    "strong Kozak context is not present as planted" =
      !row_i$has_strong_kozak ||
        grepl(
          pattern = paste0(kozak_strong, "ATGG"),
          x = paste0(row_i$utr5, row_i$cds),
          fixed = TRUE
        ),
    "TOP motif is not at the 5' end as planted" =
      (substr(x = row_i$utr5, start = 1L, stop = nchar(x = top_motif)) ==
        top_motif) == row_i$has_top_motif
  )
}

# Planted-feature totals required by the dataset design.
stopifnot(
  "expected 3 constructs with upstream AUGs" =
    sum(simulatedConstructs$n_upstream_aug > 0L) == 3L,
  "expected 3 constructs with a strong Kozak context" =
    sum(simulatedConstructs$has_strong_kozak) == 3L,
  "expected 3 constructs with a weak Kozak context" =
    sum(simulatedConstructs$has_weak_kozak) == 3L,
  "expected 2 constructs with a G-quadruplex-prone run" =
    sum(simulatedConstructs$has_g_quadruplex) == 2L,
  "expected 2 constructs with a TOP motif" =
    sum(simulatedConstructs$has_top_motif) == 2L,
  "expected 3 constructs with AU-rich elements" =
    sum(simulatedConstructs$n_are_pentamer > 0L) == 3L,
  "expected 2 constructs with a disrupted polyA signal" =
    sum(!simulatedConstructs$has_canonical_polya) == 2L,
  "expected 2 constructs with a premature in-frame stop" =
    sum(simulatedConstructs$has_premature_stop) == 2L
)

usethis::use_data(simulatedConstructs, overwrite = TRUE)

# The shipped FASTA example holds three constructs chosen to span the
# planted features: uORFs, a strong Kozak context with a TOP motif, and
# AU-rich elements.
example_ids <- c("CONSTRUCT_01", "CONSTRUCT_06", "CONSTRUCT_12")
example_rows <- simulatedConstructs[
  match(x = example_ids, table = simulatedConstructs$construct_id),
]

wrapSequence <- function(sequence, width = 60L) {
  starts <- seq(from = 1L, to = nchar(x = sequence), by = width)
  substring(
    text = sequence,
    first = starts,
    last = pmin(starts + width - 1L, nchar(x = sequence))
  )
}

fasta_lines <- character(length = 0L)
for (i in seq_len(length.out = nrow(x = example_rows))) {
  row_i <- example_rows[i, ]
  for (region in c("utr5", "cds", "utr3")) {
    fasta_lines <- c(
      fasta_lines,
      paste0(">", row_i$construct_id, "|region=", region),
      wrapSequence(sequence = row_i[[region]])
    )
  }
}

writeLines(
  text = fasta_lines,
  con = file.path("inst", "extdata", "exampleConstructs.fasta")
)
