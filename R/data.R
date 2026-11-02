#' Simulated mRNA Construct Designs With Known Planted Features
#'
#' @description
#' Twelve synthetic candidate mRNA construct designs, each supplied as a
#' 5' untranslated region (UTR), a coding sequence (CDS) and a 3' UTR,
#' together with columns recording the sequence features deliberately
#' planted in them. Because the planted features are known exactly, this
#' dataset gives a ground truth against which feature detection can be
#' checked, and it provides a worked set of candidates to rank and plot.
#'
#' @section Caution:
#' These are **simulated sequences**. They were drawn from a uniform
#' random base composition and then had motifs spliced into them. They
#' are not derived from any real transcript, they have no measured
#' expression, translation or stability outcomes attached, and nothing
#' about them should be read as biological evidence. They exist to
#' exercise and test the package's feature detection, ranking and
#' plotting on inputs whose answers are known in advance.
#'
#' @format A data frame with 12 rows and 13 columns:
#' \describe{
#'   \item{construct_id}{Character. Construct identifier, of the form
#'     \code{"CONSTRUCT_01"} to \code{"CONSTRUCT_12"}.}
#'   \item{utr5}{Character. The 5' UTR sequence, 85 to 200 nucleotides,
#'     written as upper-case DNA (\code{A}, \code{C}, \code{G},
#'     \code{T}).}
#'   \item{cds}{Character. The coding sequence, 300 to 585 nucleotides.
#'     Always a multiple of three, always beginning with \code{ATG} and
#'     ending in a stop codon, and free of in-frame stop codons except
#'     where one was planted (see \code{has_premature_stop}).}
#'   \item{utr3}{Character. The 3' UTR sequence, 200 to 500
#'     nucleotides.}
#'   \item{n_upstream_aug}{Integer. Number of upstream \code{AUG}
#'     codons planted in the 5' UTR, from 0 to 3.}
#'   \item{n_upstream_orf}{Integer. Number of those upstream AUGs that
#'     begin a complete upstream open reading frame, meaning the
#'     in-frame stop codon also lies inside the 5' UTR. Always less
#'     than or equal to \code{n_upstream_aug}; the remaining upstream
#'     AUGs have no in-frame stop before the UTR ends and so overlap the
#'     CDS.}
#'   \item{has_strong_kozak}{Logical. \code{TRUE} if a strong Kozak
#'     context was planted around the start codon, giving the purine at
#'     position -3 and the \code{G} at position +4 that mark efficient
#'     initiation (the \code{GCCACCAUGG} arrangement).}
#'   \item{has_weak_kozak}{Logical. \code{TRUE} if a deliberately weak
#'     Kozak context was planted, with a pyrimidine at position -3 and
#'     no \code{G} at position +4. Constructs where both
#'     \code{has_strong_kozak} and \code{has_weak_kozak} are
#'     \code{FALSE} carry an unmodified, intermediate context.}
#'   \item{has_g_quadruplex}{Logical. \code{TRUE} if a
#'     G-quadruplex-prone run of four G-triplets was planted in the
#'     5' UTR.}
#'   \item{has_top_motif}{Logical. \code{TRUE} if a 5' terminal
#'     oligopyrimidine (TOP) motif was planted at the very 5' end of
#'     the 5' UTR.}
#'   \item{n_are_pentamer}{Integer. Number of \code{AUUUA} pentamers
#'     planted in the 3' UTR as AU-rich elements, counting overlapping
#'     pentamers separately: 0, 3, 6 or 9. Each planted run contributes
#'     three overlapping pentamers.}
#'   \item{has_canonical_polya}{Logical. \code{TRUE} if the 3' UTR
#'     carries the canonical \code{AAUAAA} polyadenylation signal,
#'     \code{FALSE} if that signal was disrupted to \code{AAGAAA}.}
#'   \item{has_premature_stop}{Logical. \code{TRUE} if a premature
#'     in-frame stop codon was planted about 40\% of the way into the
#'     CDS. Recorded for completeness and plotted as context; CDS
#'     features never contribute to a ranking score.}
#' }
#'
#' @details
#' Features are spread across the twelve constructs so that each one has
#' both positive and negative examples, and so that most constructs
#' carry more than one planted feature. Background sequence is screened
#' so that it cannot imitate a planted feature: the 5' UTR backgrounds
#' contain no \code{ATG}, no four-G run and no G-quadruplex-prone
#' arrangement of four G-triplets, and the 3' UTR backgrounds contain
#' no \code{AUUUA} pentamer and no \code{AAUAAA} hexamer. Every
#' planted feature is then counted back out of the finished sequences
#' and checked against the recorded value, so the truth columns cannot
#' drift from the sequences they describe.
#'
#' Note that \code{n_upstream_aug} and \code{n_are_pentamer} count
#' overlapping occurrences separately, which is the convention the
#' package's own detection functions follow.
#'
#' @source
#' Simulated by \code{data-raw/makeSimulatedConstructs.R} in the
#' UTRunway package sources, under a fixed random seed so that the
#' dataset is exactly reproducible. No external or experimental data
#' were used.
#'
#' @references
#' Kozak, M. (1987). An analysis of 5'-noncoding sequences from 699
#' vertebrate messenger RNAs. \emph{Nucleic Acids Research} 15(20),
#' 8125-8148. \doi{10.1093/nar/15.20.8125}
#'
#' Shaw, G. and Kamen, R. (1986). A conserved AU sequence from the 3'
#' untranslated region of GM-CSF mRNA mediates selective mRNA
#' degradation. \emph{Cell} 46(5), 659-667.
#' \doi{10.1016/0092-8674(86)90341-7}
#'
#' Proudfoot, N. J. (2011). Ending the message: poly(A) signals then and
#' now. \emph{Genes & Development} 25(17), 1770-1782.
#' \doi{10.1101/gad.17268411}
#'
#' @examples
#' data(simulatedConstructs)
#'
#' # The construct sequences and their planted ground truth.
#' str(simulatedConstructs)
#'
#' # Which constructs carry upstream open reading frames?
#' simulatedConstructs[
#'   simulatedConstructs$n_upstream_orf > 0,
#'   c("construct_id", "n_upstream_aug", "n_upstream_orf")
#' ]
#'
#' # Region lengths, to show the shape of the candidate designs.
#' data.frame(
#'   construct_id = simulatedConstructs$construct_id,
#'   utr5_nt = nchar(x = simulatedConstructs$utr5),
#'   cds_nt = nchar(x = simulatedConstructs$cds),
#'   utr3_nt = nchar(x = simulatedConstructs$utr3)
#' )
"simulatedConstructs"
