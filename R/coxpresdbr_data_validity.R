###############################################################################

#' Checks that a User-provided gene_statistics data-frame conforms to the
#' required standard
#'
#' The data-frame should have columns `gene_id:char`, `p_value:numeric` and
#' `direction:numeric` and should have at most one row for each gene_id.
#'
#' @param        x             A putative gene_statistics data-frame
#'
.is_gene_statistics_df <- function(x) {
  is.data.frame(x) &&
    ncol(x) > 0 &&
    nrow(x) > 0 &&
    all(c("gene_id", "p_value", "direction") %in% colnames(x)) &&
    all(!duplicated(x$gene_id))
}

###############################################################################
