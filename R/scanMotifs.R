#' Scan A Sequence For Regulatory Motifs
#'
#' A function that finds motif occurrences in one region of a
#' construct. Two motifs are built in: the 5' terminal
#' oligopyrimidine (TOP) motif, which is anchored to the first
#' position and scanned only for a 5' UTR, and the
#' G-quadruplex-prone run of four G-triplets separated by short
#' loops. Any set of nucleotide motifs may be supplied instead,
#' including RNA-binding protein motifs written with IUPAC ambiguity
#' codes. Overlapping matches are all reported. Built-in motifs are
#' regular expressions because they need anchoring and quantifiers;
#' supplied motifs are nucleotide patterns and are rejected if they
#' contain anything outside the IUPAC alphabet.
#'
#' @param sequence A character vector of length one giving the
#'    sequence to scan. Upper or lower case is accepted, as is the
#'    DNA or the RNA alphabet.
#'
#' @param region A character vector of length one naming the region
#'    the sequence comes from, one of "utr5", "cds" or "utr3".
#'    It is matched exactly, with no partial matching, and selects
#'    which built-in motifs apply.
#'
#' @param motifs A named character vector of IUPAC nucleotide motifs,
#'    or NULL to scan the built-in set. Supplying motifs replaces the
#'    built-in set rather than adding to it. Default is NULL.
#'
#' @return Returns a data frame of class "data.frame" with one row
#'    per match, ordered by position, containing:
#' \itemize{
#'   \item motif - A character value giving the name of the motif
#'      that matched.
#'   \item region - A character value carrying the region argument
#'      through, so results from several regions can be combined.
#'   \item start, end - Integers giving the 1-based positions of the
#'      match's first and last nucleotides.
#'   \item width - An integer giving the match length.
#'   \item match - A character value giving the matched sequence as
#'      it occurs, which for an ambiguous motif is the concrete
#'      sequence found.
#'   \item pattern - A character value giving the pattern searched.
#' }
#'    When nothing matches, the same data frame is returned with zero
#'    rows, never NULL.
#'
#' @examples
#' # CONSTRUCT_06 has a planted TOP motif at the 5' end.
#' topRow <- UTRunway::simulatedConstructs[6, ]
#' UTRunway::scanMotifs(sequence = topRow$utr5, region = "utr5")
#'
#' # The same sequence scanned as a 3' UTR: the anchored TOP motif
#' # no longer applies, so it is not reported.
#' UTRunway::scanMotifs(sequence = topRow$utr5, region = "utr3")
#'
#' # Supplied motifs, with IUPAC ambiguity codes. Y is C or T, so
#' # this pattern matches both forms of the site.
#' UTRunway::scanMotifs(sequence = topRow$utr3, region = "utr3",
#'                      motifs = c(AU_rich = "ATTTA",
#'                                 ambiguous = "ACTAAY"))
#'
#' @references
#' Huppert, J. L. and Balasubramanian, S. (2005). Prevalence of
#' quadruplexes in the human genome. \emph{Nucleic Acids Research},
#' 33(9), 2908-2916.
#' \href{https://doi.org/10.1093/nar/gki609}{Link}
#'
#' Meyuhas, O. and Kahan, T. (2015). The race to decipher the top
#' secrets of TOP mRNAs. \emph{Biochimica et Biophysica Acta},
#' 1849(7), 801-811.
#' \href{https://doi.org/10.1016/j.bbagrm.2014.08.015}{Link}
#'
#' Ray, D., Kazan, H., Cook, K. B., Weirauch, M. T., Najafabadi, H.
#' S., Li, X., et al. (2013). A compendium of RNA-binding motifs for
#' decoding gene regulation. \emph{Nature}, 499(7457), 172-177.
#' \href{https://doi.org/10.1038/nature12311}{Link}
#'
#' @export
#' @importFrom Biostrings DNAString matchPattern
#' @importFrom IRanges end start
scanMotifs <- function(sequence, region, motifs = NULL) {
  sequence <- normaliseSequence(sequence = sequence, arg_name = "sequence")
  region <- validateRegion(region = region, arg_name = "region")

  if (is.null(x = motifs)) {
    return(scanBuiltInMotifs(sequence = sequence, region = region))
  }

  motifs <- validateMotifVector(motifs = motifs, arg_name = "motifs")

  hits <- lapply(
    X = names(x = motifs),
    FUN = function(motif_name) {
      found <- matchIupacMotif(
        sequence = sequence, motif = motifs[[motif_name]]
      )
      motifHitFrame(
        motif = motif_name,
        region = region,
        found = found,
        pattern = motifs[[motif_name]]
      )
    }
  )

  orderMotifHits(hits = do.call(what = rbind, args = hits))
}

#' Scan the built-in motif set
#'
#' @param sequence An upper-case DNA character string of length one.
#' @param region One of "utr5", "cds" or "utr3".
#' @return A motif hit data frame.
#' @noRd
scanBuiltInMotifs <- function(sequence, region) {
  hits <- list()

  # The TOP motif is defined by sitting at the transcription start
  # site, so it only means anything at the 5' end of a 5' UTR.
  if (region == "utr5") {
    hits[["TOP"]] <- motifHitFrame(
      motif = "TOP",
      region = region,
      found = matchRegexPattern(
        sequence = sequence, pattern = topMotifPattern()
      ),
      pattern = topMotifPattern()
    )
  }

  hits[["G_quadruplex"]] <- motifHitFrame(
    motif = "G_quadruplex",
    region = region,
    found = matchRegexPattern(
      sequence = sequence, pattern = gQuadruplexPattern()
    ),
    pattern = gQuadruplexPattern()
  )

  orderMotifHits(hits = do.call(what = rbind, args = hits))
}

#' Assemble the hit rows for one motif
#'
#' @param motif Name of the motif.
#' @param region Region the sequence came from.
#' @param found A list of `start`, `end` and `match` vectors.
#' @param pattern The pattern that was searched.
#' @return A motif hit data frame, possibly with zero rows.
#' @noRd
motifHitFrame <- function(motif, region, found, pattern) {
  n_hits <- length(x = found$start)
  if (n_hits == 0L) {
    return(emptyMotifTable())
  }

  data.frame(
    motif = rep(x = motif, times = n_hits),
    region = rep(x = region, times = n_hits),
    start = found$start,
    end = found$end,
    width = found$end - found$start + 1L,
    match = found$match,
    pattern = rep(x = pattern, times = n_hits),
    stringsAsFactors = FALSE
  )
}

#' Put motif hits in positional order and reset row names
#'
#' @param hits A motif hit data frame, or NULL.
#' @return A motif hit data frame.
#' @noRd
orderMotifHits <- function(hits) {
  if (is.null(x = hits) || nrow(x = hits) == 0L) {
    return(emptyMotifTable())
  }

  ordered <- hits[order(hits$start, hits$motif), , drop = FALSE]
  rownames(x = ordered) <- NULL
  ordered
}

#' The motif table with zero rows
#'
#' @return A data frame with zero rows and seven columns.
#' @noRd
emptyMotifTable <- function() {
  data.frame(
    motif = character(length = 0L),
    region = character(length = 0L),
    start = integer(length = 0L),
    end = integer(length = 0L),
    width = integer(length = 0L),
    match = character(length = 0L),
    pattern = character(length = 0L),
    stringsAsFactors = FALSE
  )
}
