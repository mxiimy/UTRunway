#' Default miRNA Seed Match Sites
#'
#' @description
#' A small, illustrative set of 7mer-m8 seed match sites for five
#' well-characterised, broadly conserved miRNA families, used as the
#' default for \code{\link{scanDecayElements}}.
#'
#' @details
#' Each site is the reverse complement of nucleotides 2 to 8 of the
#' mature miRNA, which is the 7mer-m8 site of the TargetScan
#' convention. The derivation is given so that each entry can be
#' checked against the mature sequence in miRBase, and so that sites
#' for other miRNAs can be built the same way.
#'
#' \strong{This set is illustrative, not comprehensive.} It covers
#' five families out of hundreds, and it is not organism-specific.
#' Supply your own seeds for real work.
#'
#' @return A named character vector of seed match sites, as DNA.
#' @noRd
defaultMirnaSeeds <- function() {
  c(
    # Reverse complement of nucleotides 2-8 of each mature miRNA.
    "miR-1/miR-206" = "ACATTCC",
    "let-7" = "CTACCTC",
    "miR-21" = "ATAAGCT",
    "miR-122" = "ACACTCC",
    "miR-155" = "AGCATTA"
  )
}

#' Scan A 3' UTR For Decay And Repression Elements
#'
#' A function that finds the 3' untranslated region (UTR) elements
#' governing transcript stability and repression: polyadenylation
#' signals, AU-rich elements and miRNA seed match sites. Every
#' element is reported in its own row, and overlapping matches are
#' kept, because overlapping AUUUA pentamers are the characteristic
#' form of a strong AU-rich element rather than a redundant double
#' count. These are sequence matches only: an AAUAAA hexamer far from
#' the 3' end may not be used, and a miRNA seed match is necessary
#' but far from sufficient for repression.
#'
#' @param utr3 A character vector of length one giving the 3' UTR
#'    sequence. Upper or lower case is accepted, as is the DNA or the
#'    RNA alphabet.
#'
#' @param seeds A named character vector of miRNA seed match sites
#'    written as IUPAC nucleotide patterns, or NULL to use the
#'    default set of five well-characterised families, or
#'    character(0) to skip the miRNA scan. Default is NULL.
#'
#' @param polya_end_window A non-negative integer of length one. A
#'    polyA signal ending within this many nucleotides of the 3' end
#'    is flagged as near_three_prime_end. The window is absolute
#'    rather than a fraction of region length, because a functional
#'    signal sits 10 to 30 nucleotides upstream of the cleavage site
#'    whatever the length of the UTR. Default is 30.
#'
#' @param are_cluster_fraction A double of length one between 0 and 1
#'    giving the fraction of the region's length within which two
#'    AUUUA pentamers count as clustered. Default is 0.05.
#'
#' @return Returns a data frame of class "data.frame" with one row
#'    per element, ordered by position, containing:
#' \itemize{
#'   \item element_type - A character value, one of
#'      "polyA_signal", "ARE_pentamer" or "miRNA_seed".
#'   \item element_name - A character value giving the hexamer, the
#'      pentamer, or the miRNA family name.
#'   \item start, end - Integers giving the 1-based positions of the
#'      element's first and last nucleotides.
#'   \item width - An integer giving its length in nucleotides.
#'   \item sequence - A character value giving the matched sequence.
#'   \item distance_to_end - An integer giving the nucleotides from
#'      the element's last base to the end of the UTR.
#'   \item is_canonical - A logical, TRUE for AAUAAA and FALSE for
#'      AUUAAA, NA for other element types.
#'   \item near_three_prime_end - A logical, TRUE if the element ends
#'      within polya_end_window of the 3' end.
#'   \item in_cluster - A logical, TRUE for a pentamer with another
#'      nearby, NA for other element types.
#'   \item distance_to_nearest_are - An integer giving the distance
#'      to the nearest other pentamer, or NA.
#' }
#'    When nothing is found, the same data frame is returned with
#'    zero rows, never NULL.
#'
#' @examples
#' # CONSTRUCT_12 carries nine planted AUUUA pentamers and a
#' # canonical polyA signal.
#' rich <- UTRunway::simulatedConstructs[12, ]
#' elements <- UTRunway::scanDecayElements(utr3 = rich$utr3)
#' head(x = elements, n = 10)
#'
#' # The summaries, derived from the result.
#' sum(elements$element_type == "ARE_pentamer")
#' any(elements$is_canonical & elements$near_three_prime_end,
#'     na.rm = TRUE)
#'
#' # CONSTRUCT_08 has its polyA signal disrupted, so no canonical
#' # hexamer is found.
#' disrupted <- UTRunway::simulatedConstructs[8, ]
#' UTRunway::scanDecayElements(utr3 = disrupted$utr3)
#'
#' @references
#' Bartel, D. P. (2009). MicroRNAs: target recognition and regulatory
#' functions. \emph{Cell}, 136(2), 215-233.
#' \href{https://doi.org/10.1016/j.cell.2009.01.002}{Link}
#'
#' Lewis, B. P., Burge, C. B. and Bartel, D. P. (2005). Conserved
#' seed pairing, often flanked by adenosines, indicates that
#' thousands of human genes are microRNA targets. \emph{Cell},
#' 120(1), 15-20.
#' \href{https://doi.org/10.1016/j.cell.2004.12.035}{Link}
#'
#' Proudfoot, N. J. (2011). Ending the message: poly(A) signals then
#' and now. \emph{Genes & Development}, 25(17), 1770-1782.
#' \href{https://doi.org/10.1101/gad.17268411}{Link}
#'
#' Shaw, G. and Kamen, R. (1986). A conserved AU sequence from the 3'
#' untranslated region of GM-CSF mRNA mediates selective mRNA
#' degradation. \emph{Cell}, 46(5), 659-667.
#' \href{https://doi.org/10.1016/0092-8674(86)90341-7}{Link}
#'
#' @export
scanDecayElements <- function(utr3, seeds = NULL,
                              polya_end_window = 30,
                              are_cluster_fraction = 0.05) {
  utr3 <- normaliseSequence(sequence = utr3, arg_name = "utr3")
  polya_end_window <- validateWholeNumber(
    value = polya_end_window, arg_name = "polya_end_window",
    min_value = 0L
  )
  are_cluster_fraction <- validateProportion(
    value = are_cluster_fraction, arg_name = "are_cluster_fraction"
  )

  if (is.null(x = seeds)) {
    seeds <- defaultMirnaSeeds()
  } else if (length(x = seeds) > 0L) {
    seeds <- validateMotifVector(motifs = seeds, arg_name = "seeds")
  } else if (!is.character(x = seeds)) {
    stop(
      "`seeds` must be a named character vector, NULL, or ",
      "character(0) to skip the miRNA scan.",
      call. = FALSE
    )
  }

  utr3_length <- nchar(x = utr3)
  rows <- list()

  # PolyA signals: the canonical hexamer and the common variant.
  for (hexamer in c("AATAAA", "ATTAAA")) {
    found <- matchIupacMotif(sequence = utr3, motif = hexamer)
    rows[[hexamer]] <- decayElementFrame(
      element_type = "polyA_signal",
      element_name = hexamer,
      found = found,
      utr3_length = utr3_length,
      is_canonical = hexamer == "AATAAA"
    )
  }

  # AU-rich element pentamers.
  rows[["ARE"]] <- decayElementFrame(
    element_type = "ARE_pentamer",
    element_name = "AUUUA",
    found = matchIupacMotif(sequence = utr3, motif = "ATTTA"),
    utr3_length = utr3_length,
    is_canonical = NA
  )

  # miRNA seed match sites.
  for (seed_name in names(x = seeds)) {
    rows[[paste0("seed_", seed_name)]] <- decayElementFrame(
      element_type = "miRNA_seed",
      element_name = seed_name,
      found = matchIupacMotif(
        sequence = utr3, motif = seeds[[seed_name]]
      ),
      utr3_length = utr3_length,
      is_canonical = NA
    )
  }

  elements <- do.call(what = rbind, args = rows)
  if (is.null(x = elements) || nrow(x = elements) == 0L) {
    return(emptyDecayElementTable())
  }

  # The polyA window is absolute, because the spacing it stands for
  # does not scale with region length. Clustering is relative, so that
  # designs of different lengths compare on the same footing.
  are_cluster_window <- round(x = are_cluster_fraction * utr3_length)

  elements$near_three_prime_end <-
    elements$distance_to_end <= polya_end_window

  # Clustering is defined for AU-rich pentamers only: a pentamer
  # clusters when another one starts close by on either side.
  is_are <- elements$element_type == "ARE_pentamer"
  are_starts <- elements$start[is_are]
  elements$in_cluster <- NA
  elements$distance_to_nearest_are <- NA_integer_
  if (length(x = are_starts) > 0L) {
    nearest_distance <- vapply(
      X = are_starts,
      FUN = function(this_start) {
        others <- are_starts[are_starts != this_start]
        if (length(x = others) == 0L) {
          return(NA_integer_)
        }
        as.integer(x = min(abs(x = others - this_start)))
      },
      FUN.VALUE = integer(length = 1L)
    )
    elements$distance_to_nearest_are[is_are] <- nearest_distance
    elements$in_cluster[is_are] <- !is.na(x = nearest_distance) &
      nearest_distance <= are_cluster_window
  }

  ordered <- elements[
    order(elements$start, elements$element_type, elements$element_name), ,
    drop = FALSE
  ]
  rownames(x = ordered) <- NULL
  ordered[, names(x = emptyDecayElementTable())]
}

#' Assemble the rows for one kind of element
#'
#' @param element_type Element type label.
#' @param element_name Specific element label.
#' @param found A list of `start`, `end` and `match` vectors.
#' @param utr3_length Length of the 3' UTR, for end distances.
#' @param is_canonical Canonical flag, or NA where not applicable.
#' @return A data frame of element rows, possibly with zero rows.
#' @noRd
decayElementFrame <- function(element_type, element_name, found,
                              utr3_length, is_canonical) {
  n_hits <- length(x = found$start)
  if (n_hits == 0L) {
    return(NULL)
  }

  data.frame(
    element_type = rep(x = element_type, times = n_hits),
    element_name = rep(x = element_name, times = n_hits),
    start = found$start,
    end = found$end,
    width = found$end - found$start + 1L,
    sequence = found$match,
    distance_to_end = utr3_length - found$end,
    is_canonical = rep(x = is_canonical, times = n_hits),
    stringsAsFactors = FALSE
  )
}

#' The decay element table with zero rows
#'
#' @return A data frame with zero rows and ten columns.
#' @noRd
emptyDecayElementTable <- function() {
  data.frame(
    element_type = character(length = 0L),
    element_name = character(length = 0L),
    start = integer(length = 0L),
    end = integer(length = 0L),
    width = integer(length = 0L),
    sequence = character(length = 0L),
    distance_to_end = integer(length = 0L),
    is_canonical = logical(length = 0L),
    near_three_prime_end = logical(length = 0L),
    in_cluster = logical(length = 0L),
    distance_to_nearest_are = integer(length = 0L),
    stringsAsFactors = FALSE
  )
}
