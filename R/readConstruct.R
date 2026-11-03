#' Read And Validate Candidate mRNA Constructs
#'
#' A function that reads candidate construct designs from a FASTA
#' file or a data frame, checks that each is a usable three-part
#' construct, and returns them as a validated object the rest of the
#' package relies on without revalidating. Sequences are upper-cased
#' and U is rewritten to T, so downstream functions see only
#' upper-case DNA. Every problem found across every construct is
#' reported together in a single error, one per line, each naming the
#' construct and region at fault, so a file with several mistakes
#' needs one run to diagnose rather than one run per mistake.
#'
#' @param x Either a character vector of length one giving the path
#'    to a FASTA file, or a data frame with the columns
#'    construct_id, utr5, cds and utr3. Further columns of a data
#'    frame are carried through unchanged. An object already returned
#'    by this function is passed back unchanged. In a FASTA file each
#'    region is a separate record whose header names the construct
#'    and the region, as in ">CONSTRUCT_01|region=utr5"; the
#'    shorter form ">CONSTRUCT_01|utr5" is accepted too.
#'
#' @return Returns a data frame of class "utrunway_constructs" with
#'    one row per construct containing at least:
#' \itemize{
#'   \item construct_id - A character value giving the construct's
#'      identifier, which must be non-empty and unique.
#'   \item utr5 - A character value giving the 5' UTR as upper-case
#'      DNA.
#'   \item cds - A character value giving the coding sequence, which
#'      must be a multiple of three, begin with ATG and end with a
#'      stop codon.
#'   \item utr3 - A character value giving the 3' UTR.
#'   \item has_ambiguous - A logical, TRUE if any region contains an
#'      N. Ambiguous bases warn rather than failing, because a design
#'      may legitimately carry an undetermined position, but the
#'      feature functions require unambiguous sequence.
#' }
#'
#' @examples
#' # From the FASTA file shipped with the package.
#' fastaPath <- system.file("extdata", "exampleConstructs.fasta",
#'                          package = "UTRunway")
#' constructs <- UTRunway::readConstruct(x = fastaPath)
#' constructs$construct_id
#' nchar(x = constructs$utr5)
#'
#' # From a data frame, such as the bundled dataset. Its extra
#' # ground-truth columns are carried through.
#' fromFrame <- UTRunway::readConstruct(
#'    x = UTRunway::simulatedConstructs)
#' class(x = fromFrame)
#'
#' # RNA input is accepted and normalised to DNA.
#' UTRunway::readConstruct(x = data.frame(construct_id = "rna",
#'                                        utr5 = "gccaccgccacc",
#'                                        cds = "auggccuaa",
#'                                        utr3 = "aauaaaccc"))$cds
#'
#' @references
#' Pagès, H., Aboyoun, P., Gentleman, R. and DebRoy, S. (2026).
#' \emph{Biostrings: Efficient manipulation of biological strings}.
#' R package version 2.80.2.
#' \href{https://doi.org/10.18129/B9.bioc.Biostrings}{Link}
#'
#' Wickham, H. and Bryan, J. (2023). \emph{R Packages} (2nd ed.).
#' O'Reilly Media. \href{https://r-pkgs.org/}{Link}
#'
#' @export
readConstruct <- function(x) {
  if (inherits(x = x, what = "utrunway_constructs")) {
    return(x)
  }

  if (is.null(x = x)) {
    stop(
      "`x` must be a FASTA file path or a data frame, not NULL.",
      call. = FALSE
    )
  }

  constructs <- if (is.data.frame(x = x)) {
    constructsFromDataFrame(x = x)
  } else if (is.character(x = x)) {
    constructsFromFasta(path = x)
  } else {
    stop(
      "`x` must be a FASTA file path or a data frame, not an object ",
      "of class \"", class(x = x)[[1]], "\".",
      call. = FALSE
    )
  }

  validateConstructs(constructs = constructs)
}

#' Take construct columns from a data frame
#'
#' @param x A data frame supplied by the caller.
#' @return A data frame with the four required columns first.
#' @noRd
constructsFromDataFrame <- function(x) {
  required <- c("construct_id", "utr5", "cds", "utr3")
  missing_columns <- setdiff(x = required, y = names(x = x))
  if (length(x = missing_columns) > 0L) {
    stop(
      "`x` must have the columns ",
      paste0("`", required, "`", collapse = ", "),
      ", but ",
      if (length(x = missing_columns) == 1L) "this one is" else
        "these are",
      " missing: ",
      paste0("`", missing_columns, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(x = x) == 0L) {
    stop("`x` must contain at least one construct.", call. = FALSE)
  }

  for (column in required) {
    if (!is.character(x = x[[column]])) {
      x[[column]] <- as.character(x = x[[column]])
    }
  }

  # Keep the required columns first, then anything else the caller
  # was carrying alongside them.
  x[, c(required, setdiff(x = names(x = x), y = required)), drop = FALSE]
}

#' Read construct regions from a FASTA file
#'
#' @param path Path supplied by the caller.
#' @return A data frame with the four required columns.
#' @noRd
constructsFromFasta <- function(path) {
  if (length(x = path) != 1L) {
    stop(
      "`x` must be a single FASTA file path, but it has length ",
      length(x = path), ".",
      call. = FALSE
    )
  }
  if (is.na(x = path) || !nzchar(x = path)) {
    stop(
      "`x` must be a FASTA file path, not an empty string. ",
      "system.file() returns an empty string when a file is not ",
      "found, which may be what happened here.",
      call. = FALSE
    )
  }
  if (!file.exists(path)) {
    stop("`x` is not an existing file: ", path, call. = FALSE)
  }

  lines <- readLines(con = path, warn = FALSE)
  lines <- trimws(x = lines)
  lines <- lines[nzchar(x = lines)]

  header_at <- grep(pattern = "^>", x = lines)
  if (length(x = header_at) == 0L) {
    stop(
      "No FASTA records found in ", path,
      ". A FASTA file must have at least one header line beginning ",
      "with \">\".",
      call. = FALSE
    )
  }
  if (header_at[[1]] != 1L) {
    stop(
      "The first line of ", path, " is sequence, not a FASTA header ",
      "beginning with \">\".",
      call. = FALSE
    )
  }

  block_end <- c(header_at[-1] - 1L, length(x = lines))

  record_id <- character(length = length(x = header_at))
  record_region <- character(length = length(x = header_at))
  record_sequence <- character(length = length(x = header_at))

  for (i in seq_along(along.with = header_at)) {
    header <- sub(
      pattern = "^>", replacement = "", x = lines[[header_at[[i]]]]
    )
    parsed <- parseFastaHeader(header = header, path = path)
    record_id[[i]] <- parsed$construct_id
    record_region[[i]] <- parsed$region

    if (block_end[[i]] <= header_at[[i]]) {
      record_sequence[[i]] <- ""
    } else {
      record_sequence[[i]] <- paste(
        lines[(header_at[[i]] + 1L):block_end[[i]]],
        collapse = ""
      )
    }
  }

  assembleRegions(
    record_id = record_id,
    record_region = record_region,
    record_sequence = record_sequence,
    path = path
  )
}

#' Split a FASTA header into construct identifier and region
#'
#' @param header One header line, without its leading ">".
#' @param path The file the header came from, for error messages.
#' @return A list with `construct_id` and `region`.
#' @noRd
parseFastaHeader <- function(header, path) {
  permitted <- c("utr5", "cds", "utr3")
  parts <- trimws(
    x = strsplit(x = header, split = "|", fixed = TRUE)[[1]]
  )
  parts <- parts[nzchar(x = parts)]

  if (length(x = parts) < 2L) {
    stop(
      "FASTA header \">", header, "\" in ", path, " does not name a ",
      "region. Expected a header of the form ",
      "\">CONSTRUCT_01|region=utr5\", where the region is one of ",
      paste0("\"", permitted, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }

  region_part <- sub(
    pattern = "^region=", replacement = "", x = parts[-1]
  )
  is_region <- tolower(x = region_part) %in% permitted

  if (!any(is_region)) {
    stop(
      "FASTA header \">", header, "\" in ", path, " names the region ",
      "as \"", paste(parts[-1], collapse = "|"),
      "\", which is not one of ",
      paste0("\"", permitted, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (sum(is_region) > 1L) {
    stop(
      "FASTA header \">", header, "\" in ", path, " names more than ",
      "one region.",
      call. = FALSE
    )
  }

  list(
    construct_id = parts[[1]],
    region = tolower(x = region_part[is_region])
  )
}

#' Turn one record per region into one row per construct
#'
#' @param record_id Construct identifier of each record.
#' @param record_region Region of each record.
#' @param record_sequence Sequence of each record.
#' @param path The file the records came from, for error messages.
#' @return A data frame with the four required columns.
#' @noRd
assembleRegions <- function(record_id, record_region, record_sequence,
                            path) {
  permitted <- c("utr5", "cds", "utr3")
  # Keep the order the constructs first appear in the file.
  construct_ids <- unique(x = record_id)

  problems <- character(length = 0L)
  for (construct in construct_ids) {
    present <- record_region[record_id == construct]
    missing_regions <- setdiff(x = permitted, y = present)
    duplicated_regions <- unique(x = present[duplicated(x = present)])

    if (length(x = missing_regions) > 0L) {
      problems <- c(
        problems,
        paste0(
          "construct \"", construct, "\" has no ",
          paste(missing_regions, collapse = " and no "), " record"
        )
      )
    }
    if (length(x = duplicated_regions) > 0L) {
      problems <- c(
        problems,
        paste0(
          "construct \"", construct, "\" has more than one ",
          paste(duplicated_regions, collapse = " and "), " record"
        )
      )
    }
  }

  if (length(x = problems) > 0L) {
    stop(
      "Every construct needs exactly one utr5, cds and utr3 record. ",
      "In ", path, ":\n  - ",
      paste(problems, collapse = "\n  - "),
      call. = FALSE
    )
  }

  regionSequence <- function(region) {
    vapply(
      X = construct_ids,
      FUN = function(construct) {
        record_sequence[record_id == construct & record_region == region]
      },
      FUN.VALUE = character(length = 1L),
      USE.NAMES = FALSE
    )
  }

  data.frame(
    construct_id = construct_ids,
    utr5 = regionSequence(region = "utr5"),
    cds = regionSequence(region = "cds"),
    utr3 = regionSequence(region = "utr3"),
    stringsAsFactors = FALSE
  )
}

#' Normalise and validate assembled constructs
#'
#' Collects every problem across every construct and reports them in
#' one error, so that a file with several mistakes needs one pass to
#' diagnose rather than one pass per mistake.
#'
#' @param constructs A data frame with the four required columns.
#' @return A validated object of class "utrunway_constructs".
#' @noRd
validateConstructs <- function(constructs) {
  regions <- c("utr5", "cds", "utr3")
  stop_codons <- stopCodons()
  problems <- character(length = 0L)

  ids <- constructs$construct_id
  if (anyNA(x = ids) || any(!nzchar(x = ids))) {
    problems <- c(
      problems,
      paste0(
        "construct ", which(is.na(x = ids) | !nzchar(x = ids)),
        " has a missing or empty construct_id"
      )
    )
  }
  repeated <- unique(x = ids[duplicated(x = ids)])
  if (length(x = repeated) > 0L) {
    problems <- c(
      problems,
      paste0(
        "construct_id must be unique, but \"",
        paste(repeated, collapse = "\", \""), "\" ",
        if (length(x = repeated) == 1L) "is" else "are", " repeated"
      )
    )
  }

  has_ambiguous <- logical(length = nrow(x = constructs))

  for (region in regions) {
    sequences <- toupper(x = constructs[[region]])
    sequences[is.na(x = sequences)] <- ""
    sequences <- gsub(pattern = "U", replacement = "T", x = sequences)
    constructs[[region]] <- sequences

    for (i in seq_along(along.with = sequences)) {
      label <- paste0("construct \"", ids[[i]], "\", ", region)

      if (!nzchar(x = sequences[[i]])) {
        problems <- c(problems, paste0(label, " is empty"))
        next
      }

      bad_letters <- unique(x = strsplit(
        x = gsub(
          pattern = "[ACGTN]", replacement = "", x = sequences[[i]]
        ),
        split = ""
      )[[1]])
      if (length(x = bad_letters) > 0L) {
        problems <- c(
          problems,
          paste0(
            label, " contains ",
            if (length(x = bad_letters) == 1L) "the letter " else
              "the letters ",
            paste(bad_letters, collapse = ", "),
            "; only A, C, G, T, U and N are allowed"
          )
        )
        next
      }

      if (grepl(pattern = "N", x = sequences[[i]], fixed = TRUE)) {
        has_ambiguous[[i]] <- TRUE
      }
    }
  }

  # Coding sequence structure, checked only where the sequence itself
  # is usable, so that one mistake does not produce two messages.
  for (i in seq_len(length.out = nrow(x = constructs))) {
    cds <- constructs$cds[[i]]
    label <- paste0("construct \"", ids[[i]], "\", cds")
    if (!nzchar(x = cds) || grepl(pattern = "[^ACGTN]", x = cds)) {
      next
    }

    frame_is_whole <- nchar(x = cds) %% 3L == 0L
    if (!frame_is_whole) {
      problems <- c(
        problems,
        paste0(
          label, " is ", nchar(x = cds), " nucleotides long, which is ",
          "not a multiple of three"
        )
      )
    }
    if (substr(x = cds, start = 1L, stop = 3L) != "ATG") {
      problems <- c(
        problems,
        paste0(
          label, " begins with \"",
          substr(x = cds, start = 1L, stop = 3L),
          "\" rather than the start codon ATG"
        )
      )
    }

    # Which codon is last only means something once the length
    # divides by three, so this check waits for that to hold rather
    # than reporting a second problem caused by the first.
    if (frame_is_whole) {
      last_codon <- substr(
        x = cds, start = nchar(x = cds) - 2L, stop = nchar(x = cds)
      )
      if (!last_codon %in% stop_codons) {
        problems <- c(
          problems,
          paste0(
            label, " ends with \"", last_codon, "\" rather than a ",
            "stop codon (", paste(stop_codons, collapse = ", "), ")"
          )
        )
      }
    }
  }

  if (length(x = problems) > 0L) {
    stop(
      length(x = problems),
      if (length(x = problems) == 1L) " problem" else " problems",
      " found in the constructs supplied:\n  - ",
      paste(problems, collapse = "\n  - "),
      call. = FALSE
    )
  }

  if (any(has_ambiguous)) {
    warning(
      "Ambiguous base N found in ",
      sum(has_ambiguous),
      if (sum(has_ambiguous) == 1L) " construct: " else " constructs: ",
      paste0("\"", ids[has_ambiguous], "\"", collapse = ", "),
      ". These are flagged in the has_ambiguous column. The feature ",
      "functions require unambiguous sequence and will reject them.",
      call. = FALSE
    )
  }

  constructs$has_ambiguous <- has_ambiguous
  constructs <- constructs[
    ,
    c(
      "construct_id", "utr5", "cds", "utr3", "has_ambiguous",
      setdiff(
        x = names(x = constructs),
        y = c("construct_id", "utr5", "cds", "utr3", "has_ambiguous")
      )
    ),
    drop = FALSE
  ]
  rownames(x = constructs) <- NULL

  class(x = constructs) <- c("utrunway_constructs", "data.frame")
  constructs
}
