# Internal helpers shared by the sequence feature functions. None of
# these are exported.

#' Validate and normalise a nucleotide sequence argument
#'
#' Accepts a single character string of unambiguous nucleotides, in
#' either case and in either DNA or RNA alphabet, and returns it as
#' upper-case DNA so that the rest of the package only ever handles one
#' representation.
#'
#' @param sequence Value supplied by the caller.
#' @param arg_name Name of the caller's argument, for error messages.
#' @return An upper-case DNA character string of length one.
#' @noRd
normaliseSequence <- function(sequence, arg_name) {
  if (is.null(x = sequence)) {
    stop(
      "`", arg_name, "` must be a character string, not NULL.",
      call. = FALSE
    )
  }
  # Checked before the type test so that a bare NA, which is logical,
  # is reported as a missing value rather than as the wrong class.
  if (is.atomic(x = sequence) && length(x = sequence) == 1L &&
        is.na(x = sequence)) {
    stop("`", arg_name, "` must not be NA.", call. = FALSE)
  }
  if (!is.character(x = sequence)) {
    stop(
      "`", arg_name, "` must be a character string, not an object of ",
      "class \"", class(x = sequence)[[1]], "\". A Biostrings object ",
      "can be converted with as.character().",
      call. = FALSE
    )
  }
  if (length(x = sequence) != 1L) {
    stop(
      "`", arg_name, "` must be a single sequence, but it has length ",
      length(x = sequence), ".",
      call. = FALSE
    )
  }

  normalised <- gsub(
    pattern = "U", replacement = "T", x = toupper(x = sequence)
  )

  if (!nzchar(x = normalised)) {
    stop("`", arg_name, "` must not be an empty string.", call. = FALSE)
  }
  if (grepl(pattern = "[^ACGT]", x = normalised)) {
    bad_letters <- unique(x = strsplit(
      x = gsub(pattern = "[ACGT]", replacement = "", x = normalised),
      split = ""
    )[[1]])
    stop(
      "`", arg_name, "` must contain only unambiguous nucleotides ",
      "(A, C, G, T or U), but it also contains: ",
      paste(bad_letters, collapse = ", "), ".",
      call. = FALSE
    )
  }

  normalised
}

#' Validate a single whole number argument
#'
#' @param value Value supplied by the caller.
#' @param arg_name Name of the caller's argument, for error messages.
#' @param min_value Smallest value the argument may take.
#' @return The value as an integer of length one.
#' @noRd
validateWholeNumber <- function(value, arg_name, min_value = 1L) {
  if (is.null(x = value)) {
    stop(
      "`", arg_name, "` must be a single number, not NULL.",
      call. = FALSE
    )
  }
  if (!is.numeric(x = value) || length(x = value) != 1L) {
    stop(
      "`", arg_name, "` must be a single number, but it has class \"",
      class(x = value)[[1]], "\" and length ", length(x = value), ".",
      call. = FALSE
    )
  }
  if (is.na(x = value) || !is.finite(x = value)) {
    stop(
      "`", arg_name, "` must be a finite number, not ", value, ".",
      call. = FALSE
    )
  }
  if (value != round(x = value)) {
    stop(
      "`", arg_name, "` must be a whole number, not ", value, ".",
      call. = FALSE
    )
  }
  if (value < min_value) {
    stop(
      "`", arg_name, "` must be at least ", min_value, ", not ", value,
      ".",
      call. = FALSE
    )
  }

  as.integer(x = value)
}

#' Start positions of every ATG in a sequence
#'
#' Overlapping and nested matches are all reported, because an upstream
#' AUG is a feature in its own right whether or not another one lies
#' nearby.
#'
#' @param sequence An upper-case DNA character string of length one.
#' @return An integer vector of 1-based start positions, possibly empty.
#' @noRd
findAugPositions <- function(sequence) {
  n_starts <- nchar(x = sequence) - 2L
  if (n_starts < 1L) {
    return(integer(length = 0L))
  }
  starts <- seq_len(length.out = n_starts)
  starts[
    substring(text = sequence, first = starts, last = starts + 2L) == "ATG"
  ]
}

#' The three stop codons, as DNA
#'
#' @return A character vector of length three.
#' @noRd
stopCodons <- function() {
  c("TAA", "TAG", "TGA")
}

#' Position of the first in-frame stop codon at or after a start
#'
#' Reads in steps of three from `from`, and only counts a stop codon
#' that lies wholly inside `sequence`.
#'
#' @param sequence An upper-case DNA character string of length one.
#' @param from 1-based position to start reading from.
#' @return The 1-based position of the stop codon's first base, or
#'   \code{NA_integer_} if the reading frame runs off the end of the
#'   sequence without terminating.
#' @noRd
firstInFrameStop <- function(sequence, from) {
  last_start <- nchar(x = sequence) - 2L
  if (from > last_start) {
    return(NA_integer_)
  }
  codon_starts <- seq(from = from, to = last_start, by = 3L)
  codons <- substring(
    text = sequence, first = codon_starts, last = codon_starts + 2L
  )
  hit <- which(codons %in% stopCodons())
  if (length(x = hit) == 0L) {
    return(NA_integer_)
  }
  as.integer(x = codon_starts[[hit[[1]]]])
}

#' The G-quadruplex-prone sequence pattern
#'
#' Four runs of at least three G's, each separated by a loop of one to
#' seven nucleotides. Written as a regular expression because it needs
#' quantifiers, which IUPAC motif matching cannot express. The same
#' pattern is used to screen the backgrounds of the bundled dataset,
#' so detection and that dataset's ground truth agree by construction.
#'
#' @return A character string holding a regular expression.
#' @noRd
gQuadruplexPattern <- function() {
  paste0(
    "G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"
  )
}

#' The 5' terminal oligopyrimidine (TOP) motif pattern
#'
#' A C at the very 5' end followed by a run of at least four further
#' pyrimidines. Anchored, because a TOP motif is defined by sitting at
#' the transcription start site rather than by occurring anywhere.
#'
#' @return A character string holding a regular expression.
#' @noRd
topMotifPattern <- function() {
  "^C[CT]{4,}"
}

#' The IUPAC nucleotide alphabet
#'
#' @return A character string of the permitted motif letters.
#' @noRd
iupacLetters <- function() {
  "ACGTMRWSYKVHDBN"
}

#' Validate a region argument
#'
#' Checked against the permitted values exactly, with no partial
#' matching, so that a mistyped region fails rather than being guessed.
#'
#' @param region Value supplied by the caller.
#' @param arg_name Name of the caller's argument, for error messages.
#' @return The region as a character string of length one.
#' @noRd
validateRegion <- function(region, arg_name = "region") {
  permitted <- c("utr5", "cds", "utr3")

  if (is.null(x = region)) {
    stop(
      "`", arg_name, "` must be one of ",
      paste0("\"", permitted, "\"", collapse = ", "), ", not NULL.",
      call. = FALSE
    )
  }
  if (!is.character(x = region) || length(x = region) != 1L) {
    stop(
      "`", arg_name, "` must be a single character string, one of ",
      paste0("\"", permitted, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.na(x = region) || !region %in% permitted) {
    stop(
      "`", arg_name, "` must be one of ",
      paste0("\"", permitted, "\"", collapse = ", "), ", not \"",
      region, "\".",
      call. = FALSE
    )
  }

  region
}

#' Validate and normalise a named vector of IUPAC motifs
#'
#' @param motifs Value supplied by the caller.
#' @param arg_name Name of the caller's argument, for error messages.
#' @return A named character vector of upper-case DNA IUPAC motifs.
#' @noRd
validateMotifVector <- function(motifs, arg_name) {
  if (!is.character(x = motifs)) {
    stop(
      "`", arg_name, "` must be a named character vector of motifs, ",
      "not an object of class \"", class(x = motifs)[[1]], "\".",
      call. = FALSE
    )
  }
  if (length(x = motifs) == 0L) {
    stop(
      "`", arg_name, "` must contain at least one motif.",
      call. = FALSE
    )
  }
  if (anyNA(x = motifs)) {
    stop("`", arg_name, "` must not contain NA.", call. = FALSE)
  }

  motif_names <- names(x = motifs)
  if (is.null(x = motif_names) || any(!nzchar(x = motif_names))) {
    stop(
      "`", arg_name, "` must be named, so that each match can be ",
      "reported against the motif it came from.",
      call. = FALSE
    )
  }
  if (anyDuplicated(x = motif_names) != 0L) {
    stop(
      "`", arg_name, "` must have unique names, but these are ",
      "repeated: ",
      paste(unique(x = motif_names[duplicated(x = motif_names)]),
            collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  normalised <- gsub(
    pattern = "U", replacement = "T", x = toupper(x = motifs)
  )
  if (any(!nzchar(x = normalised))) {
    stop(
      "`", arg_name, "` must not contain empty motifs, but these are ",
      "empty: ",
      paste(motif_names[!nzchar(x = normalised)], collapse = ", "), ".",
      call. = FALSE
    )
  }

  permitted <- paste0("[^", iupacLetters(), "]")
  is_bad <- grepl(pattern = permitted, x = normalised)
  if (any(is_bad)) {
    stop(
      "`", arg_name, "` must contain only IUPAC nucleotide codes (",
      iupacLetters(), "), but these motifs do not: ",
      paste(motif_names[is_bad], collapse = ", "),
      ". Motifs are nucleotide patterns, not regular expressions.",
      call. = FALSE
    )
  }

  stats::setNames(object = normalised, nm = motif_names)
}

#' Find every match of an IUPAC motif, overlapping matches included
#'
#' @param sequence An upper-case DNA character string of length one.
#' @param motif An upper-case DNA IUPAC motif.
#' @return A list of equal-length `start`, `end` and `match` vectors.
#' @noRd
matchIupacMotif <- function(sequence, motif) {
  if (nchar(x = motif) > nchar(x = sequence)) {
    return(list(
      start = integer(length = 0L),
      end = integer(length = 0L),
      match = character(length = 0L)
    ))
  }

  hits <- Biostrings::matchPattern(
    pattern = motif,
    subject = Biostrings::DNAString(x = sequence),
    fixed = FALSE
  )

  list(
    start = as.integer(x = IRanges::start(x = hits)),
    end = as.integer(x = IRanges::end(x = hits)),
    match = as.character(x = hits)
  )
}

#' Find every match of a regular expression
#'
#' @param sequence An upper-case DNA character string of length one.
#' @param pattern A regular expression.
#' @return A list of equal-length `start`, `end` and `match` vectors.
#' @noRd
matchRegexPattern <- function(sequence, pattern) {
  hits <- gregexpr(pattern = pattern, text = sequence)[[1]]

  if (hits[[1]] == -1L) {
    return(list(
      start = integer(length = 0L),
      end = integer(length = 0L),
      match = character(length = 0L)
    ))
  }

  start <- as.integer(x = hits)
  end <- start + attr(x = hits, which = "match.length") - 1L

  list(
    start = start,
    end = end,
    match = substring(text = sequence, first = start, last = end)
  )
}

#' Validate a single proportion argument
#'
#' @param value Value supplied by the caller.
#' @param arg_name Name of the caller's argument, for error messages.
#' @return The value as a numeric of length one, between 0 and 1.
#' @noRd
validateProportion <- function(value, arg_name) {
  if (is.null(x = value)) {
    stop(
      "`", arg_name, "` must be a single proportion between 0 and 1, ",
      "not NULL.",
      call. = FALSE
    )
  }
  if (!is.numeric(x = value) || length(x = value) != 1L) {
    stop(
      "`", arg_name, "` must be a single number between 0 and 1, but ",
      "it has class \"", class(x = value)[[1]], "\" and length ",
      length(x = value), ".",
      call. = FALSE
    )
  }
  if (is.na(x = value) || !is.finite(x = value)) {
    stop(
      "`", arg_name, "` must be a finite proportion, not ", value, ".",
      call. = FALSE
    )
  }
  if (value < 0 || value > 1) {
    stop(
      "`", arg_name, "` must be a proportion between 0 and 1, not ",
      value, ". To express it as a percentage of region length, ",
      "divide by 100: 5 percent is 0.05.",
      call. = FALSE
    )
  }

  as.numeric(x = value)
}
