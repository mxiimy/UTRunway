#' Estimate Local RNA Structure With a Sequence-Only Proxy
#'
#' @description
#' Scans a sequence in overlapping windows and reports, for each
#' window, its GC content and whether it contains a pair of
#' reverse-complementary runs that could pair into a hairpin stem.
#' Together these give a cheap, sequence-only indication of where a
#' transcript is likely to be structured.
#'
#' @section This is a proxy, not a folding prediction:
#' \strong{This function does not fold RNA.} It performs no
#' thermodynamic calculation: there is no nearest-neighbour energy
#' model, no partition function, no minimum free energy and no base
#' pair probability anywhere in it. GC content is counted and
#' complementary runs are matched, and that is all. The output is a
#' heuristic flag, not a structure.
#'
#' What follows from that, and matters when reading the results:
#'
#' \itemize{
#'   \item \strong{Long-range interactions are missed.} Only pairings
#'     that fall inside a single window are considered, so a stem
#'     formed between two distant parts of a transcript is invisible
#'     however stable it would be. Widening \code{window} widens the
#'     reach, at the cost of blurring where the signal is.
#'   \item \strong{No stability is estimated.} A stem of eight A-U
#'     pairs counts the same as eight G-C pairs, though the latter is
#'     far more stable. GC content is reported alongside partly to
#'     offset this.
#'   \item \strong{Competing structures are ignored.} A real sequence
#'     folds into one structure out of many possible ones; this scan
#'     reports that a pairing is possible, not that it forms.
#' }
#'
#' @section Why the default stem length is eight:
#' Short complementary runs are everywhere in random sequence, so a
#' short stem requirement reports nothing but base composition. The
#' rate at which windows of the bundled simulated constructs are
#' flagged, against the rate for the same sequences shuffled, which
#' destroys any structure while preserving composition:
#'
#' \tabular{rrr}{
#'   \strong{min_stem} \tab \strong{flagged} \tab \strong{shuffled} \cr
#'   4 \tab 79\% \tab - \cr
#'   5 \tab 28\% \tab 37\% \cr
#'   7 \tab 1.3\% \tab - \cr
#'   8 \tab 0.5\% \tab 0\%
#' }
#'
#' At a five nucleotide stem, shuffled sequence is flagged at least as
#' often as real sequence: the flag carries no information about
#' structure at that setting. A chance five nucleotide pairing is
#' expected about once per forty nucleotide window, against about once
#' per hundred windows at eight. The default is therefore eight, where
#' a flag is rarer than chance and so worth looking at.
#'
#' Because any threshold is a blunt instrument, \code{max_stem_nt}
#' reports the longest complementary run actually found in each
#' window, whatever \code{min_stem} is set to. Prefer it to
#' \code{hairpin_flag} when comparing candidates: it is a measurement
#' rather than the result of a cutoff.
#'
#' Use it to compare candidate designs against one another on the same
#' settings, which is what UTRunway does with it. Do not read it as a
#' structure prediction. If a real structure prediction is needed, fold
#' the sequence with a thermodynamic folding program, which is outside
#' this package by design, because UTRunway depends on nothing beyond
#' R.
#'
#' @details
#' Windows start at position one and advance by \code{step} for as long
#' as a whole window fits, so trailing nucleotides that cannot fill a
#' window are not covered. When the sequence is shorter than
#' \code{window}, a single window spanning the whole sequence is used
#' and a warning is raised.
#'
#' Within each window, the scan finds the longest pair of runs that
#' are the reverse complement of one another and are separated by at
#' least \code{min_loop} nucleotides, which is the smallest
#' arrangement that could close a hairpin. That length is reported as
#' \code{max_stem_nt}, and \code{hairpin_flag} is simply whether it
#' reaches \code{min_stem}.
#'
#' @param sequence Character string of length one, the sequence to
#'   scan. Upper or lower case, DNA or RNA alphabet.
#' @param window Single whole number, the window width in nucleotides.
#' @param step Single whole number, how far each window advances.
#' @param min_stem Single whole number, the shortest complementary run
#'   counted as a possible stem for \code{hairpin_flag}. Defaults to
#'   8; see the section on why. It does not affect
#'   \code{max_stem_nt}.
#' @param min_loop Single whole number, the fewest nucleotides that
#'   must separate the two runs, standing in for the turn a hairpin
#'   loop has to make.
#'
#' @return A data frame with one row per window and these columns:
#'   \describe{
#'     \item{window_start}{Integer. 1-based position of the window's
#'       first nucleotide.}
#'     \item{window_end}{Integer. 1-based position of the window's last
#'       nucleotide.}
#'     \item{gc_content}{Numeric. Proportion of G and C in the window,
#'       between 0 and 1.}
#'     \item{max_stem_nt}{Integer. Length of the longest pair of
#'       reverse-complementary runs found in the window that are far
#'       enough apart to leave a loop, or \code{0} if none reach two
#'       nucleotides. Independent of \code{min_stem}.}
#'     \item{hairpin_flag}{Logical. \code{TRUE} when
#'       \code{max_stem_nt} is at least \code{min_stem}.}
#'   }
#'
#' @references
#' Babendure, J. R., Babendure, J. L., Ding, J.-H. and Tsien, R. Y.
#' (2006). Control of mammalian translation by mRNA structure near
#' caps. \emph{RNA} 12(5), 851-861. \doi{10.1261/rna.2309906}
#'
#' Hinnebusch, A. G., Ivanov, I. P. and Sonenberg, N. (2016).
#' Translational control by 5'-untranslated regions of eukaryotic
#' mRNAs. \emph{Science} 352(6292), 1413-1416.
#' \doi{10.1126/science.aad9868}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # CONSTRUCT_05 carries a planted G-quadruplex-prone run in its
#' # 5' UTR, which shows up as a GC-rich window.
#' structured <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_05",
#' ]
#' profile <- estimateStructure(sequence = structured$utr5)
#' profile
#'
#' # The most GC-rich window, as a rough guide to where structure is
#' # most likely.
#' profile[which.max(x = profile$gc_content), ]
#'
#' # max_stem_nt is a measurement, so it can be compared directly
#' # across windows without choosing a cutoff.
#' profile[which.max(x = profile$max_stem_nt), ]
#'
#' # Finer windows resolve more detail but reach shorter distances.
#' estimateStructure(sequence = structured$utr5, window = 20, step = 20)
#'
#' # A clear hairpin: two eight nucleotide complementary runs around a
#' # three nucleotide loop.
#' estimateStructure(
#'   sequence = "GGGGCCCCAAAGGGGCCCCTTTTTTTTTTTTTTTTTTTT",
#'   window = 39
#' )
#'
#' @importFrom Biostrings DNAString letterFrequency reverseComplement
#' @export
estimateStructure <- function(sequence, window = 40, step = 10,
                              min_stem = 8, min_loop = 3) {
  sequence <- normaliseSequence(sequence = sequence, arg_name = "sequence")
  window <- validateWholeNumber(
    value = window, arg_name = "window", min_value = 1L
  )
  step <- validateWholeNumber(
    value = step, arg_name = "step", min_value = 1L
  )
  min_stem <- validateWholeNumber(
    value = min_stem, arg_name = "min_stem", min_value = 2L
  )
  min_loop <- validateWholeNumber(
    value = min_loop, arg_name = "min_loop", min_value = 1L
  )

  sequence_length <- nchar(x = sequence)

  if (window > sequence_length) {
    warning(
      "`window` is ", window, " but `sequence` is only ",
      sequence_length, " nucleotides long; using a single window ",
      "spanning the whole sequence.",
      call. = FALSE
    )
    window <- sequence_length
  }

  shortest_hairpin <- 2L * min_stem + min_loop
  if (window < shortest_hairpin) {
    warning(
      "`window` is ", window, ", which is narrower than the ",
      shortest_hairpin, " nucleotides needed for two stems of ",
      min_stem, " and a loop of ", min_loop,
      "; `hairpin_flag` will be FALSE throughout.",
      call. = FALSE
    )
  }

  window_start <- seq(
    from = 1L, to = sequence_length - window + 1L, by = step
  )
  window_end <- window_start + window - 1L

  gc_content <- numeric(length = length(x = window_start))
  max_stem_nt <- integer(length = length(x = window_start))

  for (i in seq_along(along.with = window_start)) {
    window_sequence <- substr(
      x = sequence, start = window_start[[i]], stop = window_end[[i]]
    )
    gc_content[[i]] <- as.numeric(x = Biostrings::letterFrequency(
      x = Biostrings::DNAString(x = window_sequence),
      letters = "GC",
      as.prob = TRUE
    ))
    max_stem_nt[[i]] <- longestStemLength(
      sequence = window_sequence, min_loop = min_loop
    )
  }

  data.frame(
    window_start = as.integer(x = window_start),
    window_end = as.integer(x = window_end),
    gc_content = gc_content,
    max_stem_nt = max_stem_nt,
    hairpin_flag = max_stem_nt >= min_stem,
    stringsAsFactors = FALSE
  )
}

#' Longest possible hairpin stem in a window
#'
#' Finds the longest pair of runs within `sequence` that are the
#' reverse complement of one another and far enough apart to leave a
#' loop. This is a complementarity scan only: it says a pairing is
#' geometrically possible, not that it is stable or that it forms.
#'
#' Stem lengths are tried from the longest that could fit downwards,
#' so the first length that matches is the longest. The window's
#' reverse complement is taken once and the partner of each candidate
#' run is read out of it, since the reverse complement of
#' `sequence[i:(i + k - 1)]` is the slice of that reverse complement
#' running from `n - i - k + 2` to `n - i + 1`.
#'
#' @param sequence An upper-case DNA character string of length one.
#' @param min_loop Fewest nucleotides between the two runs.
#' @return A single integer: the longest stem length found, or 0.
#' @noRd
longestStemLength <- function(sequence, min_loop) {
  sequence_length <- nchar(x = sequence)

  # Two stems and a loop have to fit inside the window.
  longest_possible <- (sequence_length - min_loop) %/% 2L
  if (longest_possible < 2L) {
    return(0L)
  }

  reverse_complement <- as.character(x = Biostrings::reverseComplement(
    x = Biostrings::DNAString(x = sequence)
  ))

  for (stem_length in seq(from = longest_possible, to = 2L, by = -1L)) {
    last_stem_start <- sequence_length - stem_length + 1L

    for (stem_start in seq_len(length.out = last_stem_start)) {
      # The partner run must begin far enough along to leave a loop.
      search_start <- stem_start + stem_length + min_loop
      if (search_start + stem_length - 1L > sequence_length) {
        next
      }

      partner <- substr(
        x = reverse_complement,
        start = sequence_length - stem_start - stem_length + 2L,
        stop = sequence_length - stem_start + 1L
      )
      downstream <- substr(
        x = sequence, start = search_start, stop = sequence_length
      )
      if (grepl(pattern = partner, x = downstream, fixed = TRUE)) {
        return(as.integer(x = stem_length))
      }
    }
  }

  0L
}
