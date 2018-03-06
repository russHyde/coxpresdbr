###############################################################################
# Statistics functions for `coxpresdbr` package
#
###############################################################################

# Statistical tests should take a DGELRT (or similar), a MArrayLM or a
# data-frame containing gene_id, two-tailed p-values, direction-of-change and
# a coex_partners data-frame (as returned by get_coex_partners)

# Where a `fit` object from edgeR or limma is used, the stats functions should
# convert this into a gene-per-row data-frame containing the two-tailed pvals
# and directions and pass this into the data_frame statistics function

# As of 2018-03-06, we only combine genes by using metap::sumz approach

###############################################################################

#' Summarise whether the coexpression partners of each source-gene in the input
#' behave consistently
#'
#' TODO: description
#'
#' @param       x              A data-frame containing columns \code{gene_id},
#' \code{p_value} and \code{direction} (at least). The \code{p_value} column
#' should contain two-tailed p-values.
#'
#' @param       coex_partners   A subset of the coexpresDB.jp database
#' containing the coexpression partners of a set of source-genes. As returned
#' by get_coex_partners. Must contain columns \code{source_id} and
#' \code{target_id}.
#'
#' @param       ...            Further arguments for passing to
#' \code{metap::sumz}.
#'
#' @return      A data-frame containing a row for each source gene in the input
#' and a summary of the number of partner-genes, the average z-score across
#' all partner genes and the p-value equivalent to this z-score.
#'
#' @importFrom   methods       is
#'
#' @export
#'
evaluate_coex_partners <- function(
                                   x,
                                   coex_partners,
                                   ...) {
  stopifnot(methods::is(x, "data.frame"))
  stopifnot(all(c("gene_id", "p_value", "direction") %in% colnames(x)))
  stopifnot(methods::is(coex_partners, "data.frame"))
  stopifnot(all(c("source_id", "target_id") %in% colnames(coex_partners)))
  NULL
}
