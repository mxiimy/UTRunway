#' Profile a Coding Sequence for Context
#'
#' @description
#' Describes a construct's coding sequence: base composition overall
#' and at third codon positions, premature in-frame stop codons, and
#' AUGs that sit out of frame. These are reported so that a construct
#' map or a candidate comparison can show the coding sequence
#' alongside its untranslated regions.
#'
#' @section Reported for context, never scored:
#' \strong{Nothing this function returns contributes to any UTRunway
#' score or ranking.} UTRunway ranks candidate designs on 5' UTR and
#' 3' UTR features alone, because its scope is untranslated region
#' design: the question it answers is which UTRs to put around a
#' coding sequence that has already been chosen. Coding sequence
#' features are shown as context, so that a design can be read as a
#' whole, and are deliberately excluded from the candidate ranking
#' that \code{rankCandidates()} performs.
#'
#' This also sets the limit of what the function computes. There is no
#' codon adaptation index here, no tRNA adaptation index, no effective
#' number of codons and no codon usage table of any kind, neither
#' bundled nor required. Those measures belong to codon optimisation,
#' which is a different problem from UTR design and is already served
#' well elsewhere.
#'
#' For codon usage analysis and codon optimisation, use the
#' \strong{cubar} package on CRAN:
#' \url{https://CRAN.R-project.org/package=cubar}. It computes codon
#' adaptation indices, tRNA adaptation weights and related statistics
#' from a codon usage table, which is exactly the ground UTRunway does
#' not cover. See Liu et al. (2025) in the references below.
#'
#' @details
#' GC3 is the proportion of G or C at the third position of each
#' complete codon. It is reported separately from overall GC content
#' because the third position is the most free to vary without
#' changing the protein, so the two figures can differ considerably in
#' a sequence that has been recoded.
#'
#' A premature stop codon is an in-frame stop at any codon before the
#' last one. The final codon is treated as the natural terminator
#' whether or not it is a stop codon, which \code{ends_with_stop}
#' reports separately.
#'
#' An out-of-frame AUG is an ATG whose first base does not fall on a
#' codon boundary. Such an AUG can start a competing reading frame if
#' it is reached by a scanning ribosome, which is why it is worth
#' seeing next to the 5' UTR's upstream AUGs; see
#' \code{\link{findUpstreamORFs}}.
#'
#' @param cds Character string of length one, the coding sequence.
#'   Upper or lower case, DNA or RNA alphabet. A sequence whose length
#'   is not a multiple of three is accepted with a warning, and the
#'   trailing incomplete codon is ignored.
#'
#' @return A list of thirteen elements:
#'   \describe{
#'     \item{cds_length}{Integer. Length in nucleotides.}
#'     \item{n_codons}{Integer. Number of complete codons.}
#'     \item{gc_content}{Numeric. Proportion of G or C over the whole
#'       sequence, between 0 and 1.}
#'     \item{gc3_content}{Numeric. Proportion of G or C at third codon
#'       positions, between 0 and 1.}
#'     \item{starts_with_atg}{Logical. Whether the first codon is
#'       \code{ATG}.}
#'     \item{ends_with_stop}{Logical. Whether the last complete codon
#'       is a stop codon.}
#'     \item{is_multiple_of_three}{Logical. Whether the length divides
#'       by three.}
#'     \item{has_premature_stop}{Logical. Whether any in-frame stop
#'       codon occurs before the last codon.}
#'     \item{premature_stop_positions}{Integer vector. 1-based
#'       positions of the first base of each premature stop codon, in
#'       order. Length zero when there are none.}
#'     \item{premature_stop_codons}{Character vector. Which stop codon
#'       was found at each of those positions.}
#'     \item{n_out_of_frame_atg}{Integer. How many ATGs do not sit on
#'       a codon boundary.}
#'     \item{out_of_frame_atg_positions}{Integer vector. 1-based
#'       positions of those ATGs. Length zero when there are none.}
#'     \item{out_of_frame_atg_frames}{Integer vector. For each, its
#'       offset from the coding frame, \code{1} or \code{2}.}
#'   }
#'
#' @references
#' Ikemura, T. (1985). Codon usage and tRNA content in unicellular and
#' multicellular organisms. \emph{Molecular Biology and Evolution}
#' 2(1), 13-34. \doi{10.1093/oxfordjournals.molbev.a040335}
#'
#' Kudla, G., Murray, A. W., Tollervey, D. and Plotkin, J. B. (2009).
#' Coding-sequence determinants of gene expression in Escherichia
#' coli. \emph{Science} 324(5924), 255-258.
#' \doi{10.1126/science.1170160}
#'
#' Liu, M., Zi, B., Zhang, H. and Zhang, H. (2025). cubar: a versatile
#' package for codon usage bias analysis in R. \emph{Genetics},
#' iyaf191. \doi{10.1093/genetics/iyaf191}
#'
#' For codon usage analysis and codon optimisation, which this package
#' deliberately leaves out, see the cubar package:
#' \url{https://CRAN.R-project.org/package=cubar}.
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # A construct with no planted premature stop codon.
#' clean <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_01",
#' ]
#' profileCDS(cds = clean$cds)
#'
#' # CONSTRUCT_03 has one planted about 40 percent of the way in.
#' truncated <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_03",
#' ]
#' profile <- profileCDS(cds = truncated$cds)
#' profile$has_premature_stop
#' profile$premature_stop_positions
#' profile$premature_stop_codons
#'
#' # GC content and GC3 compared across every candidate.
#' data.frame(
#'   construct_id = simulatedConstructs$construct_id,
#'   gc = vapply(
#'     X = simulatedConstructs$cds,
#'     FUN = function(cds) profileCDS(cds = cds)$gc_content,
#'     FUN.VALUE = numeric(1),
#'     USE.NAMES = FALSE
#'   ),
#'   gc3 = vapply(
#'     X = simulatedConstructs$cds,
#'     FUN = function(cds) profileCDS(cds = cds)$gc3_content,
#'     FUN.VALUE = numeric(1),
#'     USE.NAMES = FALSE
#'   )
#' )
#'
#' @importFrom Biostrings DNAString letterFrequency
#' @export
profileCDS <- function(cds) {
  cds <- normaliseSequence(sequence = cds, arg_name = "cds")

  cds_length <- nchar(x = cds)
  if (cds_length < 3L) {
    stop(
      "`cds` must contain at least one complete codon, but it is ",
      cds_length, " nucleotide", if (cds_length == 1L) "" else "s",
      " long.",
      call. = FALSE
    )
  }

  is_multiple_of_three <- cds_length %% 3L == 0L
  if (!is_multiple_of_three) {
    warning(
      "`cds` is ", cds_length, " nucleotides long, which is not a ",
      "multiple of three; the trailing ", cds_length %% 3L,
      " nucleotide", if (cds_length %% 3L == 1L) "" else "s",
      " will be ignored.",
      call. = FALSE
    )
  }

  n_codons <- cds_length %/% 3L
  codon_starts <- seq(
    from = 1L, by = 3L, length.out = n_codons
  )
  codons <- substring(
    text = cds, first = codon_starts, last = codon_starts + 2L
  )

  gc_content <- as.numeric(x = Biostrings::letterFrequency(
    x = Biostrings::DNAString(x = cds), letters = "GC", as.prob = TRUE
  ))

  third_positions <- codon_starts + 2L
  third_bases <- paste(
    substring(text = cds, first = third_positions, last = third_positions),
    collapse = ""
  )
  gc3_content <- as.numeric(x = Biostrings::letterFrequency(
    x = Biostrings::DNAString(x = third_bases),
    letters = "GC",
    as.prob = TRUE
  ))

  # The last complete codon is the natural terminator, so only stops
  # before it count as premature.
  internal_is_stop <- codons[-n_codons] %in% stopCodons()
  premature_index <- which(internal_is_stop)

  atg_positions <- findAugPositions(sequence = cds)
  atg_offset <- (atg_positions - 1L) %% 3L
  out_of_frame <- atg_offset != 0L

  list(
    cds_length = as.integer(x = cds_length),
    n_codons = as.integer(x = n_codons),
    gc_content = gc_content,
    gc3_content = gc3_content,
    starts_with_atg = codons[[1]] == "ATG",
    ends_with_stop = codons[[n_codons]] %in% stopCodons(),
    is_multiple_of_three = is_multiple_of_three,
    has_premature_stop = length(x = premature_index) > 0L,
    premature_stop_positions = as.integer(
      x = codon_starts[premature_index]
    ),
    premature_stop_codons = codons[premature_index],
    n_out_of_frame_atg = sum(out_of_frame),
    out_of_frame_atg_positions = as.integer(
      x = atg_positions[out_of_frame]
    ),
    out_of_frame_atg_frames = as.integer(x = atg_offset[out_of_frame])
  )
}
