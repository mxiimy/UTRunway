#' UTRunway: Profile and Rank Candidate mRNA Construct Designs
#'
#' @description
#' UTRunway profiles candidate messenger RNA (mRNA) construct designs
#' supplied as three parts: a 5' untranslated region (UTR), a coding
#' sequence (CDS) and a 3' UTR. It computes sequence features known to
#' affect translation initiation and transcript stability, ranks candidate
#' designs against one another, and draws the whole construct as a single
#' annotated multi-track map.
#'
#' @section Scope:
#' UTRunway works on sequence alone. It requires no sequencing data, no
#' genome annotation and no genomic coordinates.
#'
#' \itemize{
#'   \item Ranking scores use 5' UTR and 3' UTR features only. CDS features
#'     are reported and plotted as context but never contribute to a score.
#'   \item RNA secondary structure is estimated with a lightweight in-R
#'     proxy (sliding-window GC content plus a simple complementarity
#'     scan). It is not a thermodynamic folding prediction.
#'   \item All computation is done in R using CRAN and Bioconductor
#'     packages; no external programs are required.
#' }
#'
#' @section Relation to other packages:
#' ORFik analyses translation using ribosome profiling data; cubar
#' analyses and optimises codon usage in coding sequences; RegVar
#' interprets 3' UTR variants at genomic positions. UTRunway instead
#' compares whole candidate constructs from their sequence, with scoring
#' focused on the untranslated regions.
#'
#' @author Melody Mo \email{melody.mo@@mail.utoronto.ca}
#'
#' @references
#' Liu, M., Zi, B., Zhang, H. and Zhang, H. (2025). cubar: a versatile
#' package for codon usage bias analysis in R. \emph{Genetics},
#' iyaf191. \doi{10.1093/genetics/iyaf191}
#'
#' @keywords internal
"_PACKAGE"
