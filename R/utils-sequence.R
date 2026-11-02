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
