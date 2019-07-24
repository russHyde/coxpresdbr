###############################################################################
# Parsing functions for `coxpresdbr` package
#
###############################################################################

.define_default_gene_universe <- function(coex_df) {
  stopifnot(
    all(c("source_id", "target_id") %in% colnames(coex_df))
  )
  unique(
    c(
      coex_df[["source_id"]],
      coex_df[["target_id"]]
    )
  )
}

###############################################################################

#' Function to filter the contents of a 'coxpresdb.jp' gene-coexpression
#' dataset
#'
#' @param        coex_df       A dataframe containing coexpression data. As
#'   returned by \code{import_all_coex_partners}.
#'
#' @param        gene_universe   The genes in the dataframe should be filtered
#'   to ensure they are all present in this set. Note that both the entries in
#'   \code{source_id} and \code{target_id} are filtered to be in this set.
#'
#' @param        n_partners    The maximum number of partners to return for the
#'   source gene(s) in \code{coex_df}. Gene partners will be sorted in order of
#'   mutual-rank before selecting the top partners.
#'
#' @param        mr_threshold   All gene-pairs in the returned dataset will
#'   have a mutual-rank of at most this value.
#'
#' @importFrom   dplyr         arrange_   filter_
#' @importFrom   utils         head
#'
.filter_coex_partners <- function(
                                  coex_df,
                                  gene_universe,
                                  n_partners,
                                  mr_threshold) {
  if (length(unique(coex_df[["source_id"]])) > 1) {
    stop(
      "`.filter_coex_partners` has not yet been implemented for >1 source gene"
    )
  }

  gene_universe <- if (missing(gene_universe) || is.null(gene_universe)) {
    .define_default_gene_universe(coex_df)
  } else {
    gene_universe
  }

  n_partners <- if (missing(n_partners) || is.null(n_partners)) {
    nrow(coex_df)
  } else {
    n_partners
  }

  mr_threshold <- if (missing(mr_threshold) || is.null(mr_threshold)) {
    max(coex_df[["mutual_rank"]])
  } else {
    mr_threshold
  }

  coex_df %>%
    dplyr::filter_(
      ~ source_id %in% gene_universe &
        target_id %in% gene_universe &
        mutual_rank <= mr_threshold
    ) %>%
    dplyr::arrange_(~mutual_rank) %>%
    head(n = n_partners)
}

###############################################################################

#' Returns the coexpression partners of a given set of genes
#'
#' Several filters are applied while identifying the coexpression partners of a
#' given set of genes. Only those gene-pairs with a mutual rank below- the
#' user-selected threshold, and for a given input gene _at most_
#' \code{n_partners} are returned.
#'
#' If the user provides a \code{gene_universe}, the dataframe will only contain
#' genes that are present in this set (the default universe is the entirety of
#' the genes present in the coxpresDB archive).
#'
#' If the thresholds, number-of-partners or gene_universe are set to
#' \code{NULL}, the corresponding filters are not applied to the dataset.
#'
#' @param        gene_ids      A vector of gene identifiers. All of these
#'   should be annotated within the 'coxpresDb' dataset.
#'
#' @param        importer      A \code{CoxpresDbImporter}. This allows access
#'   to a 'coxpresDb' archive.
#'
#' @include      coxpresdbr_io.R
#'
#' @inheritParams   .filter_coex_partners
#'
#' @importFrom   dplyr         bind_rows
#'
#' @export
#'
get_coex_partners <- function(
                              gene_ids,
                              importer,
                              gene_universe = NULL,
                              n_partners = 100,
                              mr_threshold = NULL) {
  .import_fn <- function(x) {
    import_all_coex_partners(gene_id = x, importer = importer)
  }

  .filter_fn <- function(x) {
    .filter_coex_partners(
      x,
      gene_universe = gene_universe, n_partners = n_partners,
      mr_threshold = mr_threshold
    )
  }

  imported <- Map(.import_fn, gene_ids)
  filtered <- Map(.filter_fn, imported)

  dplyr::bind_rows(filtered)
}

###############################################################################
