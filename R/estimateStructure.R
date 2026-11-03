#' Estimate Local RNA Structure With A Sequence-Only Proxy
#'
#' A function that scans a sequence in overlapping windows and
#' reports, for each window, its GC content and the longest pair of
#' reverse-complementary runs that could pair into a hairpin stem.
#' Together these give a cheap indication of where a transcript is
#' likely to be structured. This function does not fold RNA: there is
#' no nearest-neighbour energy model, no partition function and no
#' minimum free energy anywhere in it. Pairings spanning more than
#' one window are invisible, a stem of A-U pairs counts the same as
#' one of G-C pairs, and a possible pairing is not one that forms.
#' Use it to compare candidate designs on the same settings, and
#' fold the sequence with a thermodynamic program if a structure
#' prediction is needed.
#'
#' @param sequence A character vector of length one giving the
#'    sequence to scan. Upper or lower case is accepted, as is the
#'    DNA or the RNA alphabet.
#'
#' @param window A positive integer of length one giving the window
#'    width in nucleotides. Default is 40.
#'
#' @param step A positive integer of length one giving how far each
#'    window advances. Default is 10.
#'
#' @param min_stem A positive integer of length one giving the
#'    shortest complementary run counted as a possible stem for
#'    hairpin_flag. At a five nucleotide stem, shuffled sequence is
#'    flagged as often as real sequence, so the default is stricter.
#'    Default is 8.
#'
#' @param min_loop A positive integer of length one giving the fewest
#'    nucleotides that must separate the two runs. Default is 3.
#'
#' @return Returns a data frame of class "data.frame" with one row
#'    per window containing:
#' \itemize{
#'   \item window_start - An integer giving the 1-based position of
#'      the window's first nucleotide.
#'   \item window_end - An integer giving the position of its last
#'      nucleotide.
#'   \item gc_content - A double giving the proportion of G and C in
#'      the window, between 0 and 1.
#'   \item max_stem_nt - An integer giving the longest pair of
#'      reverse-complementary runs found, or 0. Independent of
#'      min_stem, and preferable to hairpin_flag when comparing
#'      candidates because it is a measurement rather than a cutoff.
#'   \item hairpin_flag - A logical, TRUE when max_stem_nt is at
#'      least min_stem.
#' }
#'
#' @examples
#' # CONSTRUCT_05 carries a planted G-quadruplex-prone run in its
#' # 5' UTR, which shows as a GC-rich stretch.
#' structured <- UTRunway::simulatedConstructs[5, ]
#' structureProfile <- UTRunway::estimateStructure(
#'    sequence = structured$utr5)
#' structureProfile
#'
#' # The most GC-rich window, as a rough guide to where structure is
#' # most likely.
#' structureProfile[which.max(x = structureProfile$gc_content), ]
#'
#' # Finer windows resolve more detail but reach shorter distances.
#' UTRunway::estimateStructure(sequence = structured$utr5,
#'                             window = 20, step = 20)
#'
#' @references
#' Babendure, J. R., Babendure, J. L., Ding, J.-H. and Tsien, R. Y.
#' (2006). Control of mammalian translation by mRNA structure near
#' caps. \emph{RNA}, 12(5), 851-861.
#' \href{https://doi.org/10.1261/rna.2309906}{Link}
#'
#' Hinnebusch, A. G., Ivanov, I. P. and Sonenberg, N. (2016).
#' Translational control by 5'-untranslated regions of eukaryotic
#' mRNAs. \emph{Science}, 352(6292), 1413-1416.
#' \href{https://doi.org/10.1126/science.aad9868}{Link}
#'
#' @export
#' @importFrom Biostrings DNAString letterFrequency reverseComplement
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
