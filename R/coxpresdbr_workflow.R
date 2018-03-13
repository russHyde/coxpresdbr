###############################################################################

# Workflows for use in `coxpresdbr` package

###############################################################################

# TODO: `run_coex_partner_workflow`
run_coex_partner_workflow <- function(
                                      gene_ids,
                                      gene_statistics,
                                      importer,
                                      gene_universe = NULL) {
  if (missing(gene_ids) || is.null(gene_ids)) {
    stop("gene_ids should be deined in run_coex_partner_workflow")
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
    stop("importer should be defined in run_coex_partner_workflow")
  }
}

###############################################################################
