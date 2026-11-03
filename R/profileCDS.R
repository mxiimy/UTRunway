#' Profile A Coding Sequence For Context
#'
#' A function that describes a construct's coding sequence: base
#' composition overall and at third codon positions, premature
#' in-frame stop codons, and AUGs sitting out of frame. Nothing this
#' function returns contributes to any UTRunway score or ranking.
#' The package ranks candidate designs on 5' and 3' untranslated
#' region features alone, because its scope is untranslated region
#' design: which UTRs to put around a coding sequence already
#' chosen. Coding sequence features are shown as context so a design
#' can be read as a whole. There is deliberately no codon adaptation
#' index here and no codon usage table; for codon usage analysis and
#' codon optimisation use the cubar package on CRAN.
#'
#' @param cds A character vector of length one giving the coding
#'    sequence. Upper or lower case is accepted, as is the DNA or the
#'    RNA alphabet. A length that is not a multiple of three is
#'    accepted with a warning, and the trailing incomplete codon is
#'    ignored.
#'
#' @return Returns a list of length thirteen containing:
#' \itemize{
#'   \item cds_length - An integer giving the length in nucleotides.
#'   \item n_codons - An integer giving the number of complete
#'      codons.
#'   \item gc_content - A double giving the proportion of G or C over
#'      the whole sequence.
#'   \item gc3_content - A double giving the proportion of G or C at
#'      third codon positions, reported separately because that
#'      position is freest to vary without changing the protein.
#'   \item starts_with_atg - A logical, TRUE if the first codon is
#'      ATG.
#'   \item ends_with_stop - A logical, TRUE if the last complete
#'      codon is a stop codon.
#'   \item is_multiple_of_three - A logical, TRUE if the length
#'      divides by three.
#'   \item has_premature_stop - A logical, TRUE if an in-frame stop
#'      codon occurs before the last codon.
#'   \item premature_stop_positions - An integer vector of the
#'      1-based positions of those stop codons, length zero if none.
#'   \item premature_stop_codons - A character vector naming which
#'      stop codon was found at each position.
#'   \item n_out_of_frame_atg - An integer giving how many ATGs do
#'      not sit on a codon boundary.
#'   \item out_of_frame_atg_positions - An integer vector of their
#'      positions, length zero if none.
#'   \item out_of_frame_atg_frames - An integer vector giving each
#'      one's offset from the coding frame, 1 or 2.
#' }
#'
#' @examples
#' # A construct with no planted premature stop codon.
#' clean <- UTRunway::simulatedConstructs[1, ]
#' UTRunway::profileCDS(cds = clean$cds)
#'
#' # CONSTRUCT_03 has one planted about 40 percent of the way in.
#' truncated <- UTRunway::simulatedConstructs[3, ]
#' cdsProfile <- UTRunway::profileCDS(cds = truncated$cds)
#' cdsProfile$has_premature_stop
#' cdsProfile$premature_stop_positions
#' cdsProfile$premature_stop_codons
#'
#' @references
#' Ikemura, T. (1985). Codon usage and tRNA content in unicellular
#' and multicellular organisms. \emph{Molecular Biology and
#' Evolution}, 2(1), 13-34.
#' \href{https://doi.org/10.1093/oxfordjournals.molbev.a040335}{Link}
#'
#' Kudla, G., Murray, A. W., Tollervey, D. and Plotkin, J. B. (2009).
#' Coding-sequence determinants of gene expression in Escherichia
#' coli. \emph{Science}, 324(5924), 255-258.
#' \href{https://doi.org/10.1126/science.1170160}{Link}
#'
#' Liu, M., Zi, B., Zhang, H. and Zhang, H. (2025). cubar: a
#' versatile package for codon usage bias analysis in R.
#' \emph{Genetics}, iyaf191.
#' \href{https://doi.org/10.1093/genetics/iyaf191}{Link}
#'
#' @export
#' @importFrom Biostrings DNAString letterFrequency
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
