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
#' @param        coex_df       A dataframe containing coexpression data. Should
#'   contain columns \code{source_id}, \code{target_id} and \code{mutual_rank}.
#'   May contain multiple distinct \code{source_id}s.
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
#' @importFrom   dplyr         arrange
#' @importFrom   rlang         .data
#' @importFrom   utils         head
#'
.filter_coex_partners <- function(
                                  coex_df,
                                  gene_universe,
                                  n_partners,
                                  mr_threshold) {
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

  keep_rows <- with(
    coex_df,
    which(
      source_id %in% gene_universe &
        target_id %in% gene_universe &
        source_id != target_id &
        mutual_rank <= mr_threshold
    )
  )

  coex_df[keep_rows, ] %>%
    dplyr::group_by(
      .data[["source_id"]]
    ) %>%
    dplyr::top_n(
      n = -n_partners, wt = .data[["mutual_rank"]]
    ) %>%
    dplyr::arrange(
      .data[["mutual_rank"]]
    ) %>%
    dplyr::ungroup()
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
#' @param        importer      A \code{CoxpresDbAccessor}. This allows access
#'   to a 'coxpresDb' archive.
#'
#' @include      coxpresdbr_io.R
#'
#' @inheritParams   .filter_coex_partners
#'
#' @importFrom   dplyr         bind_rows
#' @importFrom   methods       is
#' @importFrom   purrr         map_df
#'
#' @export
#'
get_coex_partners <- function(
                              gene_ids,
                              importer,
                              gene_universe = NULL,
                              n_partners = 100,
                              mr_threshold = NULL) {
  .are_genes_valid <- function() {
    !any(duplicated(gene_ids)) &
      all(gene_ids %in% get_source_ids(importer))
  }

  stopifnot(
    is(importer, "CoxpresDbAccessor") &&
      .are_genes_valid()
  )

  imported <- if (is(importer, "CoxpresDbArchiveAccessor")) {
    get_all_coex_partners(
      gene_ids,
      importer = importer
    )
  } else {
    # must be a CoxpresDbDataframeAccessor
    # - extract the source-gene rows directly from the dataframe
    rows <- which(importer@df[["source_id"]] %in% gene_ids)

    importer@df[rows, ]
  }

  .filter_coex_partners(
    imported,
    gene_universe = gene_universe, n_partners = n_partners,
    mr_threshold = mr_threshold
  )
}

###############################################################################
