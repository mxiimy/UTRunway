#' Locate Upstream AUGs and Upstream Open Reading Frames in a 5' UTR
#'
#' @description
#' Finds every AUG in a 5' UTR and describes what each one does: the
#' frame it sits in relative to the main start codon, whether it opens
#' an upstream open reading frame (uORF) that terminates before the
#' coding sequence begins, and how strong its own initiation context
#' is.
#'
#' @details
#' Every AUG is reported in its own row, whether or not other AUGs lie
#' nearby, so overlapping and nested uORFs are all described rather
#' than collapsed into one call. A uORF nested inside a longer one
#' appears as a separate row, as does each AUG of a pair that share a
#' stop codon.
#'
#' An upstream AUG falls into one of three cases, distinguished by
#' \code{opens_uorf} and \code{frame}:
#'
#' \itemize{
#'   \item It reaches an in-frame stop codon inside the 5' UTR, so it
#'     opens a self-contained uORF. \code{opens_uorf} is \code{TRUE}.
#'   \item It has no in-frame stop codon before the UTR ends and is out
#'     of frame with the main start codon, so its reading frame runs
#'     through the start codon and into the coding sequence in the
#'     wrong frame. \code{out_of_frame_overlap} is \code{TRUE}.
#'   \item It has no in-frame stop codon before the UTR ends and is in
#'     frame with the main start codon, so it acts as an N-terminal
#'     extension of the main reading frame. \code{extends_into_cds} is
#'     \code{TRUE} and \code{frame} is \code{0}.
#' }
#'
#' Positions are 1-based and are given relative to the 5' UTR, not to
#' the whole construct.
#'
#' @param utr5 Character string of length one, the 5' UTR sequence.
#'   Upper or lower case, DNA or RNA alphabet.
#' @param cds Character string of length one, the coding sequence, or
#'   \code{NULL}. The coding sequence is used only to complete the
#'   Kozak context window of an upstream AUG lying within four
#'   nucleotides of the end of the 5' UTR. Supplying it does not change
#'   which AUGs or uORFs are found. When it is \code{NULL}, such an AUG
#'   gets \code{NA} for its Kozak columns.
#'
#' @return A data frame with one row per upstream AUG, ordered by
#'   position, and these columns:
#'   \describe{
#'     \item{aug_position}{Integer. 1-based position in \code{utr5} of
#'       the A of the upstream AUG.}
#'     \item{distance_to_cds}{Integer. Nucleotides from the A of the
#'       upstream AUG to the A of the main start codon.}
#'     \item{frame}{Integer. Frame relative to the main start codon:
#'       \code{0} if in frame with it, otherwise \code{1} or \code{2}.}
#'     \item{in_frame_with_cds}{Logical. \code{TRUE} when \code{frame}
#'       is \code{0}.}
#'     \item{opens_uorf}{Logical. \code{TRUE} if an in-frame stop codon
#'       lies wholly inside the 5' UTR, so the uORF terminates before
#'       the coding sequence.}
#'     \item{stop_position}{Integer. 1-based position in \code{utr5} of
#'       the first base of that stop codon, or \code{NA_integer_}.}
#'     \item{uorf_length_nt}{Integer. Length of the uORF in
#'       nucleotides from the A of the AUG through the last base of the
#'       stop codon inclusive, or \code{NA_integer_}.}
#'     \item{extends_into_cds}{Logical. \code{TRUE} if the reading
#'       frame from this AUG reaches the end of the 5' UTR without
#'       terminating, so it continues past the main start codon.}
#'     \item{out_of_frame_overlap}{Logical. \code{TRUE} if the AUG both
#'       extends into the coding sequence and is out of frame with the
#'       main start codon.}
#'     \item{kozak_score}{Numeric. Log-odds Kozak score of this
#'       upstream start codon, from \code{\link{scoreKozak}}, or
#'       \code{NA_real_} where the context window does not fit.}
#'     \item{kozak_strength}{Character. \code{"strong"},
#'       \code{"adequate"}, \code{"weak"} or \code{NA_character_}.}
#'   }
#'   When the 5' UTR contains no AUG, the same data frame is returned
#'   with zero rows and these columns and types, never \code{NULL}.
#'
#' @section Kozak context near the 5' end:
#' An upstream AUG within six nucleotides of the 5' end has no complete
#' context window, so its Kozak columns are \code{NA}. This is expected
#' rather than an error, and the warning \code{\link{scoreKozak}} would
#' raise is suppressed here.
#'
#' @references
#' Calvo, S. E., Pervouchine, D. D. and Mootha, V. K. (2009). Upstream
#' open reading frames cause widespread reduction of protein expression
#' and are polymorphic among humans. \emph{Proceedings of the National
#' Academy of Sciences} 106(18), 7507-7512.
#' \doi{10.1073/pnas.0810916106}
#'
#' Hinnebusch, A. G., Ivanov, I. P. and Sonenberg, N. (2016).
#' Translational control by 5'-untranslated regions of eukaryotic
#' mRNAs. \emph{Science} 352(6292), 1413-1416.
#' \doi{10.1126/science.aad9868}
#'
#' Kozak, M. (1987). An analysis of 5'-noncoding sequences from 699
#' vertebrate messenger RNAs. \emph{Nucleic Acids Research} 15(20),
#' 8125-8148. \doi{10.1093/nar/15.20.8125}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # CONSTRUCT_03 has three planted upstream AUGs, two of which open
#' # uORFs that terminate inside the 5' UTR.
#' with_uorfs <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_03",
#' ]
#' findUpstreamORFs(utr5 = with_uorfs$utr5, cds = with_uorfs$cds)
#'
#' # CONSTRUCT_05 has none, so the result has zero rows but the same
#' # columns.
#' without_uorfs <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_05",
#' ]
#' findUpstreamORFs(utr5 = without_uorfs$utr5)
#'
#' # Nested and overlapping starts are each reported in their own row.
#' findUpstreamORFs(utr5 = "ATGATGAAATAACCGGTTACGTACGTACG")
#'
#' @export
findUpstreamORFs <- function(utr5, cds = NULL) {
  utr5 <- normaliseSequence(sequence = utr5, arg_name = "utr5")
  if (!is.null(x = cds)) {
    cds <- normaliseSequence(sequence = cds, arg_name = "cds")
  }

  utr5_length <- nchar(x = utr5)
  aug_positions <- findAugPositions(sequence = utr5)

  if (length(x = aug_positions) == 0L) {
    return(emptyUpstreamOrfTable())
  }

  # The Kozak context of an upstream AUG can reach past the end of the
  # 5' UTR, so score against the UTR and CDS joined when the coding
  # sequence is available.
  context_sequence <- if (is.null(x = cds)) utr5 else paste0(utr5, cds)

  n_aug <- length(x = aug_positions)
  distance_to_cds <- utr5_length + 1L - aug_positions
  frame <- as.integer(x = distance_to_cds %% 3L)

  stop_position <- integer(length = n_aug)
  kozak_score <- numeric(length = n_aug)
  kozak_strength <- character(length = n_aug)

  for (i in seq_len(length.out = n_aug)) {
    stop_position[[i]] <- firstInFrameStop(
      sequence = utr5, from = aug_positions[[i]]
    )

    # NA for an AUG too close to the 5' end is expected here, so the
    # warning scoreKozak() raises is not passed on to the caller.
    kozak <- suppressWarnings(scoreKozak(
      sequence = context_sequence, position = aug_positions[[i]]
    ))
    kozak_score[[i]] <- kozak$score
    kozak_strength[[i]] <- kozak$strength
  }

  opens_uorf <- !is.na(x = stop_position)
  uorf_length_nt <- ifelse(
    test = opens_uorf,
    yes = stop_position + 3L - aug_positions,
    no = NA_integer_
  )

  data.frame(
    aug_position = as.integer(x = aug_positions),
    distance_to_cds = as.integer(x = distance_to_cds),
    frame = frame,
    in_frame_with_cds = frame == 0L,
    opens_uorf = opens_uorf,
    stop_position = as.integer(x = stop_position),
    uorf_length_nt = as.integer(x = uorf_length_nt),
    extends_into_cds = !opens_uorf,
    out_of_frame_overlap = !opens_uorf & frame != 0L,
    kozak_score = kozak_score,
    kozak_strength = kozak_strength,
    stringsAsFactors = FALSE
  )
}

#' The upstream ORF table with zero rows
#'
#' Keeps the no-AUG case returning the same columns and types as a
#' populated result, so that callers never have to special-case it.
#'
#' @return A data frame with zero rows and eleven columns.
#' @noRd
emptyUpstreamOrfTable <- function() {
  data.frame(
    aug_position = integer(length = 0L),
    distance_to_cds = integer(length = 0L),
    frame = integer(length = 0L),
    in_frame_with_cds = logical(length = 0L),
    opens_uorf = logical(length = 0L),
    stop_position = integer(length = 0L),
    uorf_length_nt = integer(length = 0L),
    extends_into_cds = logical(length = 0L),
    out_of_frame_overlap = logical(length = 0L),
    kozak_score = numeric(length = 0L),
    kozak_strength = character(length = 0L),
    stringsAsFactors = FALSE
  )
}
