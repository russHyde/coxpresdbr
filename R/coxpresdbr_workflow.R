###############################################################################

# Workflows for use in `coxpresdbr` package

###############################################################################

#' Extract partners of a set of gene-ids, generate summary statistics over
#' those neighbours and then cluster the initial gene-ids
#'
#' @param        gene_ids     A set of gene identifiers. Should be a subset of
#' the genes in \code{gene_statistics} and (if defined) \code{gene_universe}.
#'
#' @param        gene_statistics   A table containing the two-tailed p-values,
#' direction-of-change and gene-ids for all genes that are to be considered in
#' the analysis (any genes that are not in the `gene_universe` are
#' disregarded).
#'
#' @param        importer      A CoxpresDbImporter object.
#'
#' @param        gene_universe   A set of gene identifiers or NULL. If defined,
#' each of the gene-identifiers should have data present in a row of
#' \code{gene_statistics}. If undefined, the universe is taken to be the
#' intersection of all genes in \code{gene_statistics} and all genes accessible
#' through the CoxpresDB dataset in \code{importer}. Once defined, all genes in
#' \code{gene_ids} should be present in \code{gene_universe}.
#'
#' @param        n_partners    The maximum number of partners to pull out for
#' a given source gene from the CoxpresDb database.
#'
#' @param        ...           Other arguments passed to
#' \code{get_coex_partners}.
#'
#' @importFrom   dplyr         filter_
#'
#' @export
#'
run_coex_partner_workflow <- function(
                                      gene_ids,
                                      gene_statistics,
                                      importer,
                                      gene_universe = NULL,
                                      n_partners = 100, ...) {
  if (missing(gene_ids) || is.null(gene_ids)) {
    stop("`gene_ids` should be defined in `run_coex_partner_workflow`")
  }

  if (missing(gene_statistics) ||
    is.null(gene_statistics) ||
    !.is_gene_statistics_df(gene_statistics)
  ) {
    stop(paste(
      "`gene_statistics` should be defined and pass `.is_gene_statistics_df`",
      "in `run_coex_partner_workflow`"
    ))
  }

  if (missing(importer) || !(is(importer, "CoxpresDbImporter"))) {
    stop("`importer` should be defined in `run_coex_partner_workflow`")
  }

  if (is.null(gene_universe)) {
    gene_universe <- intersect(
      gene_statistics$gene_id, get_gene_ids(importer)
    )
  }

  if (!all(gene_universe %in% gene_statistics[["gene_id"]])) {
    stop(paste(
      "If `gene_universe` is provided, all it's elements should be",
      "present in `gene_statistics$gene_id`"
    ))
  }

  if (!all(gene_ids %in% gene_universe)) {
    stop(paste(
      "If `gene_universe` is provided, ensure every one of `gene_ids` is",
      "present inside `gene_universe`"
    ))
  }

  gene_statistics <- dplyr::filter_(
      gene_statistics, ~ gene_id %in% gene_universe
    )

  partners <- get_coex_partners(
      gene_ids, importer, gene_universe, n_partners, ...
    )

  partner_summaries <- evaluate_coex_partners(
    gene_statistics, partners
    )

  new(
    "CoxpresDbPartners",
    gene_statistics = gene_statistics,
    partners = partners,
    partner_summaries = partner_summaries
  )
}

###############################################################################
