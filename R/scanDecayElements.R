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

#' Scan a 3' UTR for Decay and Repression Elements
#'
#' @description
#' Finds the 3' UTR elements that govern transcript stability and
#' repression: polyadenylation signals, AU-rich elements, and miRNA
#' seed match sites.
#'
#' @details
#' Every element found is reported in its own row. Overlapping matches
#' are all kept, because overlapping AUUUA pentamers are the
#' characteristic form of a strong AU-rich element rather than a
#' redundant double count.
#'
#' Three element types are reported in \code{element_type}:
#'
#' \describe{
#'   \item{\code{"polyA_signal"}}{The canonical hexamer \code{AATAAA}
#'     and the common variant \code{ATTAAA}. \code{is_canonical}
#'     distinguishes them. A functional signal normally lies a short
#'     distance upstream of the cleavage site, so
#'     \code{near_three_prime_end} flags signals ending within
#'     \code{polya_end_window} nucleotides of the 3' end.}
#'   \item{\code{"ARE_pentamer"}}{The \code{AUUUA} pentamer of an
#'     AU-rich element. \code{in_cluster} flags a pentamer with
#'     another pentamer starting within \code{are_cluster_fraction} of
#'     the region's length, since clustered and overlapping repeats
#'     drive decay far more strongly than an isolated pentamer does.}
#'   \item{\code{"miRNA_seed"}}{A match to one of the seed sites in
#'     \code{seeds}. \code{element_name} gives the miRNA family.}
#' }
#'
#' The two summaries the caller usually wants come straight off the
#' result. Whether a canonical polyA signal is present near the 3' end:
#'
#' \preformatted{any(elements$is_canonical &
#'     elements$near_three_prime_end, na.rm = TRUE)}
#'
#' and the AU-rich element count, with whether any of them cluster:
#'
#' \preformatted{sum(elements$element_type == "ARE_pentamer")
#' any(elements$in_cluster, na.rm = TRUE)}
#'
#' @section A match is not a functional site:
#' These are sequence matches. An \code{AATAAA} hexamer far from the 3'
#' end may not be used, and a miRNA seed match is necessary but nowhere
#' near sufficient for repression, which also depends on site context,
#' accessibility and whether the miRNA is expressed at all. Treat a
#' seed match as a candidate site to consider, never as a validated
#' target.
#'
#' @section On the two thresholds:
#' The two thresholds are deliberately expressed in different units,
#' because the two pieces of biology behave differently.
#'
#' \code{polya_end_window} is an absolute number of nucleotides. A
#' functional polyadenylation signal sits roughly 10 to 30 nucleotides
#' upstream of the cleavage site whether the 3' UTR is 200 or 2000
#' nucleotides long, so a window expressed as a fraction of region
#' length would be too tight on a short UTR and far too generous on a
#' long one.
#'
#' \code{are_cluster_fraction} is a fraction of the region's length,
#' so that candidate designs of different lengths are compared on the
#' same relative footing. Overlapping pentamers lie four nucleotides
#' apart and so cluster under any sensible threshold; the fraction
#' matters only for deciding how far apart two separate runs may be
#' and still count together.
#'
#' Both raw distances are reported as columns, so any other rule can
#' be applied to the result without rescanning. For a strict 20
#' nucleotide polyA window, or clustering within an absolute 40
#' nucleotides:
#'
#' \preformatted{is_signal <- elements$element_type == "polyA_signal"
#' any(is_signal & elements$is_canonical &
#'       elements$distance_to_end <= 20)
#' elements$distance_to_nearest_are <= 40}
#'
#' @param utr3 Character string of length one, the 3' UTR sequence.
#'   Upper or lower case, DNA or RNA alphabet.
#' @param seeds Named character vector of miRNA seed match sites, as
#'   IUPAC nucleotide patterns, or \code{NULL} to use the default set
#'   of five well-characterised families. Pass
#'   \code{character(0)} to skip the miRNA scan entirely.
#' @param polya_end_window Single whole number of nucleotides. A
#'   polyA signal ending within this many nucleotides of the 3' end is
#'   flagged as \code{near_three_prime_end}. The default of \code{30}
#'   covers the 10 to 30 nucleotide spacing at which a functional
#'   signal sits upstream of the cleavage site. This window is
#'   absolute rather than a fraction of the region's length, because
#'   that spacing does not scale with the length of the 3' UTR.
#' @param are_cluster_fraction Single number between 0 and 1. Two
#'   AUUUA pentamers starting within this fraction of the region's
#'   length count as clustered. Defaults to \code{0.05}.
#'
#' @return A data frame with one row per element, ordered by position,
#'   and these columns:
#'   \describe{
#'     \item{element_type}{Character. \code{"polyA_signal"},
#'       \code{"ARE_pentamer"} or \code{"miRNA_seed"}.}
#'     \item{element_name}{Character. The specific element: the hexamer
#'       for a polyA signal, \code{"AUUUA"} for an AU-rich pentamer, or
#'       the miRNA family name for a seed match.}
#'     \item{start}{Integer. 1-based position of the first nucleotide.}
#'     \item{end}{Integer. 1-based position of the last nucleotide.}
#'     \item{width}{Integer. Element length in nucleotides.}
#'     \item{sequence}{Character. The matched sequence as it occurs in
#'       \code{utr3}.}
#'     \item{distance_to_end}{Integer. Nucleotides from the element's
#'       last base to the last base of the 3' UTR.}
#'     \item{is_canonical}{Logical. For a polyA signal, \code{TRUE} for
#'       \code{AATAAA} and \code{FALSE} for \code{ATTAAA};
#'       \code{NA} for other element types.}
#'     \item{near_three_prime_end}{Logical. \code{TRUE} if the element
#'       ends within \code{polya_end_window} nucleotides of the 3'
#'       end. Reported for every element type, since it is a plain
#'       distance test, but it is interpreted chiefly for polyA
#'       signals.}
#'     \item{in_cluster}{Logical. For an AU-rich pentamer, \code{TRUE}
#'       if another pentamer starts within
#'       \code{are_cluster_fraction} of the region's length;
#'       \code{NA} for other element types, and \code{FALSE} for a
#'       lone pentamer.}
#'     \item{distance_to_nearest_are}{Integer. For an AU-rich
#'       pentamer, the distance in nucleotides to the nearest other
#'       pentamer's start. \code{NA} for other element types and for a
#'       lone pentamer. Reported so that a different clustering rule
#'       can be applied without rescanning.}
#'   }
#'   When nothing is found, the same data frame is returned with zero
#'   rows and these columns and types, never \code{NULL}.
#'
#' @references
#' Proudfoot, N. J. (2011). Ending the message: poly(A) signals then
#' and now. \emph{Genes & Development} 25(17), 1770-1782.
#' \doi{10.1101/gad.17268411}
#'
#' Beaudoing, E., Freier, S., Wyatt, J. R., Claverie, J. M. and
#' Gautheret, D. (2000). Patterns of variant polyadenylation signal
#' usage in human genes. \emph{Genome Research} 10(7), 1001-1010.
#' \doi{10.1101/gr.10.7.1001}
#'
#' Shaw, G. and Kamen, R. (1986). A conserved AU sequence from the 3'
#' untranslated region of GM-CSF mRNA mediates selective mRNA
#' degradation. \emph{Cell} 46(5), 659-667.
#' \doi{10.1016/0092-8674(86)90341-7}
#'
#' Bakheet, T., Hitti, E., Al-Saif, M., Moghrabi, W. N. and Khabar, K.
#' S. A. (2018). The AU-rich element landscape across human
#' transcriptome reveals a large proportion in introns and regulation
#' by ELAVL1/HuR. \emph{Biochimica et Biophysica Acta} 1861(2),
#' 167-177. \doi{10.1016/j.bbagrm.2017.12.006}
#'
#' Lewis, B. P., Burge, C. B. and Bartel, D. P. (2005). Conserved seed
#' pairing, often flanked by adenosines, indicates that thousands of
#' human genes are microRNA targets. \emph{Cell} 120(1), 15-20.
#' \doi{10.1016/j.cell.2004.12.035}
#'
#' Bartel, D. P. (2009). MicroRNAs: target recognition and regulatory
#' functions. \emph{Cell} 136(2), 215-233.
#' \doi{10.1016/j.cell.2009.01.002}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # CONSTRUCT_12 has nine planted AUUUA pentamers and a canonical
#' # polyA signal.
#' rich <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_12",
#' ]
#' elements <- scanDecayElements(utr3 = rich$utr3)
#' head(x = elements, n = 10)
#'
#' # The summaries, derived from the result.
#' sum(elements$element_type == "ARE_pentamer")
#' any(elements$in_cluster, na.rm = TRUE)
#' any(elements$is_canonical & elements$near_three_prime_end,
#'     na.rm = TRUE)
#'
#' # CONSTRUCT_08 has its polyA signal disrupted, so no canonical
#' # hexamer is found.
#' disrupted <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_08",
#' ]
#' scanDecayElements(utr3 = disrupted$utr3)
#'
#' # A supplied seed set, with an IUPAC ambiguity code.
#' scanDecayElements(
#'   utr3 = rich$utr3,
#'   seeds = c(my_mirna = "ACATTCC", ambiguous = "ACTAAY")
#' )
#'
#' # Skip the miRNA scan.
#' nrow(x = scanDecayElements(utr3 = rich$utr3, seeds = character(0)))
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
