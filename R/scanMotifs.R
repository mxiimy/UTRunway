#' Scan a Sequence for Regulatory Motifs
#'
#' @description
#' Finds motif occurrences in one region of a construct. Two motifs are
#' built in, the 5' terminal oligopyrimidine (TOP) motif and the
#' G-quadruplex-prone run, and any set of nucleotide motifs can be
#' supplied instead, including RNA-binding protein motifs written with
#' IUPAC ambiguity codes.
#'
#' @details
#' Overlapping matches are all reported, since a motif occurring twice
#' in overlapping register is two sites.
#'
#' Two kinds of pattern are matched, and which one is used depends on
#' where the motif came from:
#'
#' \itemize{
#'   \item \strong{Built-in motifs} are regular expressions, because
#'     they need quantifiers and anchoring that nucleotide patterns
#'     cannot express.
#'   \item \strong{Supplied motifs} are nucleotide patterns in the
#'     IUPAC alphabet (\code{ACGTMRWSYKVHDBN}), matched with
#'     \code{\link[Biostrings]{matchPattern}}. So \code{"TGTANATA"}
#'     matches with any base at the \code{N}, and \code{"ACTAAY"}
#'     matches \code{ACTAAC} or \code{ACTAAT}. They are
#'     \strong{not} regular expressions: a motif containing anything
#'     outside the IUPAC alphabet is rejected rather than being matched
#'     as a regular expression.
#' }
#'
#' Supplying \code{motifs} \strong{replaces} the built-in set rather
#' than adding to it, so that the result holds exactly what was asked
#' for. To get both, call the function twice and \code{rbind} the
#' results.
#'
#' @section Built-in motifs:
#' \describe{
#'   \item{\code{TOP}}{A 5' terminal oligopyrimidine motif: a C at the
#'     first position of the region followed by at least four further
#'     pyrimidines. Because it is defined by sitting at the
#'     transcription start site, it is anchored to position 1 and is
#'     only scanned when \code{region} is \code{"utr5"}. It is silently
#'     skipped for the other regions, where it has no meaning.}
#'   \item{\code{G_quadruplex}}{Four runs of at least three G's, each
#'     pair separated by a loop of one to seven nucleotides, which is
#'     the arrangement conventionally treated as quadruplex-prone.
#'     Scanned in every region. A sequence match says a quadruplex is
#'     possible, not that one forms; see
#'     \code{\link{estimateStructure}} for the same caution about
#'     sequence-only structure proxies.}
#' }
#'
#' @param sequence Character string of length one, the sequence to
#'   scan. Upper or lower case, DNA or RNA alphabet.
#' @param region Character string naming the region \code{sequence}
#'   comes from: \code{"utr5"}, \code{"cds"} or \code{"utr3"}. Matched
#'   exactly, with no partial matching. It selects which built-in
#'   motifs apply and is recorded in the result.
#' @param motifs Named character vector of IUPAC nucleotide motifs, or
#'   \code{NULL} to scan the built-in set. Names are used to label the
#'   matches and must be unique and non-empty. An empty vector is an
#'   error rather than a silent scan of nothing: pass \code{NULL} to
#'   ask for the built-in motifs.
#'
#' @return A data frame with one row per match, ordered by position,
#'   and these columns:
#'   \describe{
#'     \item{motif}{Character. Name of the motif that matched.}
#'     \item{region}{Character. The \code{region} argument, carried
#'       through so results from several regions can be combined.}
#'     \item{start}{Integer. 1-based position of the match's first
#'       nucleotide.}
#'     \item{end}{Integer. 1-based position of its last nucleotide.}
#'     \item{width}{Integer. Match length in nucleotides.}
#'     \item{match}{Character. The matched sequence as it occurs,
#'       which for an ambiguous motif is the concrete sequence found.}
#'     \item{pattern}{Character. The pattern that was searched, so the
#'       result records what was actually looked for.}
#'   }
#'   When nothing matches, the same data frame is returned with zero
#'   rows and these columns and types, never \code{NULL}.
#'
#' @references
#' Meyuhas, O. and Kahan, T. (2015). The race to decipher the top
#' secrets of TOP mRNAs. \emph{Biochimica et Biophysica Acta} 1849(7),
#' 801-811. \doi{10.1016/j.bbagrm.2014.08.015}
#'
#' Huppert, J. L. and Balasubramanian, S. (2005). Prevalence of
#' quadruplexes in the human genome. \emph{Nucleic Acids Research}
#' 33(9), 2908-2916. \doi{10.1093/nar/gki609}
#'
#' Ray, D., Kazan, H., Cook, K. B., Weirauch, M. T., Najafabadi, H. S.,
#' Li, X., et al. (2013). A compendium of RNA-binding motifs for
#' decoding gene regulation. \emph{Nature} 499(7457), 172-177.
#' \doi{10.1038/nature12311}
#'
#' @examples
#' data(simulatedConstructs, package = "UTRunway")
#'
#' # CONSTRUCT_06 has a planted TOP motif at the 5' end.
#' top_row <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_06",
#' ]
#' scanMotifs(sequence = top_row$utr5, region = "utr5")
#'
#' # CONSTRUCT_05 has a planted G-quadruplex-prone run instead.
#' quad_row <- simulatedConstructs[
#'   simulatedConstructs$construct_id == "CONSTRUCT_05",
#' ]
#' scanMotifs(sequence = quad_row$utr5, region = "utr5")
#'
#' # The same sequence scanned as a 3' UTR: the anchored TOP motif no
#' # longer applies, so it is not reported.
#' scanMotifs(sequence = top_row$utr5, region = "utr3")
#'
#' # Supplied motifs, with IUPAC ambiguity codes. Y is C or T, so this
#' # pattern matches both forms of the site.
#' scanMotifs(
#'   sequence = top_row$utr3,
#'   region = "utr3",
#'   motifs = c(PUM_like = "TGTANATA", AU_rich = "ATTTA", test = "ACTAAY")
#' )
#'
#' @importFrom Biostrings DNAString matchPattern
#' @importFrom IRanges end start
#' @export
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
