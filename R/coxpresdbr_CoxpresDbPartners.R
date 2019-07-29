###############################################################################

#' Validity checker for \code{CoxpresDbPartners} objects
#'
#' @param        object        A putative \code{CoxpresDbPartners} object
#'
#' @include      coxpresdbr_data_validity.R
#'
.validity_coxpresdb_partners <- function(object) {
  # S3 tests:
  if (!is.null(object@cluster_graph) &&
    !is(object@cluster_graph, "igraph")
  ) {
    return(
      paste(
        "`cluster_graph` should be `NULL` or inherit from `igraph` in",
        "`CoxpresDbPartners`"
      )
    )
  }

  if (
    any(dim(object@gene_statistics) > 0) &&
      !.is_gene_statistics_df(object@gene_statistics)
  ) {
    return(
      paste(
        "`gene_statistics` should be empty or have `gene_id`, `p_value` and",
        "`direction` columns"
      )
    )
  }

  if (
    any(dim(object@partners) > 0) &&
      !all(c("source_id", "target_id") %in% colnames(object@partners))
  ) {
    return(
      paste(
        "`partners` should be empty or (minimally) have `source_id` and",
        "`target_id` as column names"
      )
    )
  }
}

###############################################################################

#' @title        Template for \code{CoxpresDbPartners} class
#'
#' @description   This class stores the gene-partners of a set of genes, and
#'   some summary statistics over those partners that are obtained from
#'   analysis of a user-provided \code{gene_statistics} object.
#'
#' @param        gene_statistics   A data-frame. Must be non-empty and have
#'   columns `gene_id`, `p_value` and `direction`.
#' @param        partners      A data-frame with a `source_id` and a
#'   `target_id` column.
#' @param        partner_summaries   A data-frame.
#' @param        cluster_graph   A igraph object (or NULL if not defined).
#'
#' @name         CoxpresDbPartners-class
#' @rdname       CoxpresDbPartners-class
#'
#' @exportClass   CoxpresDbPartners
#'
methods::setClass(
  "CoxpresDbPartners",
  slots = list(
    gene_statistics = "data.frame",
    partners = "data.frame",
    partner_summaries = "data.frame",
    cluster_graph = "ANY"
  ),
  validity = function(object) .validity_coxpresdb_partners(object)
)

###############################################################################

#' Generic function for returning the set of all genes considered in the
#' analysis of gene-partners
#'
#' @param        x             A CoxpresDbPartners object, as returned by
#' \code{run_coex_partner_workflow}
#'
setGeneric("get_gene_universe", valueClass = "character", function(x) {
  standardGeneric("get_gene_universe")
})

#' Obtains the set of all genes considered in the analysis of gene-partners
#'
#' @param        x             A CoxpresDbPartners object, as returned by
#' \code{run_coex_partner_workflow}
#'
#' @export
#'
setMethod("get_gene_universe", signature("CoxpresDbPartners"), function(x) {
  x@gene_statistics$gene_id
})

###############################################################################

#' Generic function for obtaining the set of source genes for which partner
#' genes were assessed
#'
#' @param        x             A CoxpresDbPartners object, as returned by
#' \code{run_coex_partner_workflow}

setGeneric("get_source_genes", valueClass = "character", function(x) {
  standardGeneric("get_source_genes")
})

#' Obtains the set of source genes for which partner genes were assessed
#'
#' @param        x             A CoxpresDbPartners object, as returned by
#' \code{run_coex_partner_workflow}
#
#' @export
#'
setMethod("get_source_genes", signature("CoxpresDbPartners"), function(x) {
  if (nrow(x@partners) == 0) {
    character(0)
  } else {
    unique(x@partners$source_id)
  }
})

###############################################################################

.compute_z_scores <- function(p_values, directions) {
  if (
    missing(p_values) || missing(directions) ||
      any(p_values < 0) || any(p_values > 1)
  ) {
    stop()
  }
  qnorm(p_values / 2) * -1 * sign(directions)
}

setGeneric(".add_z_scores", valueClass = "CoxpresDbPartners", function(x) {
  standardGeneric(".add_z_scores")
})

setMethod(".add_z_scores", signature("CoxpresDbPartners"), function(x) {
  statistics <- x@gene_statistics
  statistics$z_score <- .compute_z_scores(
    statistics$p_value,
    statistics$direction
  )
  x@gene_statistics <- statistics
  x
})

###############################################################################
