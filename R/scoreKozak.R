#' Position weight matrix for the Kozak context
#'
#' Columns are the scored context positions relative to the A of the
#' AUG, which is position +1. The three bases of the start codon are
#' not scored because they are fixed. Probabilities are derived from
#' the \code{gccRccAUGG} consensus: the two positions known to dominate
#' initiation efficiency, the purine at -3 and the G at +4, are given
#' sharper distributions than the flanking positions.
#'
#' @return A 4 by 7 numeric matrix whose columns each sum to one, with
#'   rows named for the four bases.
#' @noRd
kozakPwm <- function() {
  pwm <- cbind(
    "-6" = c(A = 0.20, C = 0.20, G = 0.40, T = 0.20),
    "-5" = c(A = 0.20, C = 0.40, G = 0.20, T = 0.20),
    "-4" = c(A = 0.20, C = 0.40, G = 0.20, T = 0.20),
    "-3" = c(A = 0.400, C = 0.125, G = 0.350, T = 0.125),
    "-2" = c(A = 0.20, C = 0.40, G = 0.20, T = 0.20),
    "-1" = c(A = 0.20, C = 0.40, G = 0.20, T = 0.20),
    "+4" = c(A = 0.15, C = 0.15, G = 0.60, T = 0.10)
  )
  pwm
}

#' Offsets of the scored context positions, relative to the A of the AUG
#'
#' @return An integer vector of length seven.
#' @noRd
kozakOffsets <- function() {
  c(-6L, -5L, -4L, -3L, -2L, -1L, 3L)
}

#' Score The Kozak Context Of A Start Codon
#'
#' A function that scores the sequence context around a start codon
#' against a position weight matrix built from the Kozak consensus
#' gccRccAUGG, and classifies it as strong, adequate or weak. The
#' scored window runs from position -6 to -1 and then to position +4,
#' counting the A of the AUG as +1; the three bases of the start
#' codon are not scored because they are fixed. The categorical
#' strength follows the rule that a purine at -3 and a G at +4 are
#' the two positions that dominate initiation efficiency: both gives
#' "strong", one gives "adequate", neither gives "weak".
#'
#' @param sequence A character vector of length one giving the
#'    sequence containing the start codon and its context. To score a
#'    start codon in a construct, pass the 5' UTR and the coding
#'    sequence pasted together so the window can reach both sides.
#'
#' @param position A positive integer of length one giving the
#'    1-based position in sequence of the A of the AUG to score.
#'
#' @return Returns a list of length four containing:
#' \itemize{
#'   \item score - A double giving the log-odds score summed over the
#'      seven scored positions, where zero means no better than
#'      chance, or NA_real_ if the context window does not fit.
#'   \item score_scaled - A double giving the same score on a zero to
#'      one scale, or NA_real_.
#'   \item strength - A character value, one of "strong",
#'      "adequate" or "weak", or NA_character_.
#'   \item context - A character value giving the scored context
#'      written out, such as "gccacc ATG g", or NA_character_.
#' }
#'    A full context window needs six bases before the AUG and one
#'    base four positions after it. If either is missing, which is
#'    common for an upstream AUG close to the 5' end, the function
#'    warns and returns NA rather than failing.
#'
#' @examples
#' # CONSTRUCT_04 carries a planted strong Kozak context. The main
#' # start codon sits immediately after the 5' UTR.
#' strongRow <- UTRunway::simulatedConstructs[4, ]
#' UTRunway::scoreKozak(
#'    sequence = paste0(strongRow$utr5, strongRow$cds),
#'    position = nchar(x = strongRow$utr5) + 1)
#'
#' # CONSTRUCT_07 carries a planted weak context, for comparison.
#' weakRow <- UTRunway::simulatedConstructs[7, ]
#' UTRunway::scoreKozak(
#'    sequence = paste0(weakRow$utr5, weakRow$cds),
#'    position = nchar(x = weakRow$utr5) + 1)
#'
#' # A start codon too close to the 5' end for a full context window
#' # returns NA, with a warning.
#' suppressWarnings(
#'    UTRunway::scoreKozak(sequence = "AATGGCCCCC", position = 2))
#'
#' @references
#' Kozak, M. (1986). Point mutations define a sequence flanking the
#' AUG initiator codon that modulates translation by eukaryotic
#' ribosomes. \emph{Cell}, 44(2), 283-292.
#' \href{https://doi.org/10.1016/0092-8674(86)90762-2}{Link}
#'
#' Kozak, M. (1987). An analysis of 5'-noncoding sequences from 699
#' vertebrate messenger RNAs. \emph{Nucleic Acids Research}, 15(20),
#' 8125-8148. \href{https://doi.org/10.1093/nar/15.20.8125}{Link}
#'
#' @export
scoreKozak <- function(sequence, position) {
  sequence <- normaliseSequence(sequence = sequence, arg_name = "sequence")
  position <- validateWholeNumber(
    value = position, arg_name = "position", min_value = 1L
  )

  sequence_length <- nchar(x = sequence)
  if (position + 2L > sequence_length) {
    stop(
      "`position` is ", position, ", which leaves no room for a start ",
      "codon in a sequence of ", sequence_length, " nucleotides.",
      call. = FALSE
    )
  }

  start_codon <- substr(
    x = sequence, start = position, stop = position + 2L
  )
  if (start_codon != "ATG") {
    stop(
      "`position` must point at the A of an ATG start codon, but the ",
      "three nucleotides at position ", position, " are \"",
      start_codon, "\".",
      call. = FALSE
    )
  }

  offsets <- kozakOffsets()
  scored_positions <- position + offsets

  if (
    min(scored_positions) < 1L || max(scored_positions) > sequence_length
  ) {
    warning(
      "The Kozak context window (-6 to +4) does not fit inside ",
      "`sequence` for a start codon at position ", position,
      "; returning NA.",
      call. = FALSE
    )
    return(list(
      score = NA_real_,
      score_scaled = NA_real_,
      strength = NA_character_,
      context = NA_character_
    ))
  }

  pwm <- kozakPwm()
  observed_bases <- substring(
    text = sequence, first = scored_positions, last = scored_positions
  )

  probabilities <- vapply(
    X = seq_along(along.with = observed_bases),
    FUN = function(i) pwm[observed_bases[[i]], i],
    FUN.VALUE = numeric(length = 1L)
  )
  score <- sum(log2(x = probabilities / 0.25))

  # Scale against the worst and best contexts this matrix can award,
  # so that the result is comparable between candidate designs.
  best_score <- sum(log2(x = apply(X = pwm, MARGIN = 2L, FUN = max) / 0.25))
  worst_score <- sum(log2(x = apply(X = pwm, MARGIN = 2L, FUN = min) / 0.25))
  score_scaled <- (score - worst_score) / (best_score - worst_score)

  base_minus_3 <- observed_bases[[which(offsets == -3L)]]
  base_plus_4 <- observed_bases[[which(offsets == 3L)]]
  n_favourable <- sum(base_minus_3 %in% c("A", "G"), base_plus_4 == "G")
  strength <- switch(
    EXPR = as.character(x = n_favourable),
    "2" = "strong",
    "1" = "adequate",
    "0" = "weak"
  )

  context <- paste(
    tolower(x = paste(observed_bases[1:6], collapse = "")),
    "ATG",
    tolower(x = base_plus_4)
  )

  list(
    score = score,
    score_scaled = score_scaled,
    strength = strength,
    context = context
  )
}
