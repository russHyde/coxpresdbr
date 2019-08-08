###############################################################################
# Statistics functions for `coxpresdbr` package
#
###############################################################################

# Statistical tests should take a DGELRT (or similar), a MArrayLM or a
# data-frame containing gene_id, two-tailed p-values, direction-of-change and
# a coex_partners data-frame (as returned by get_coex_partners)

# Where a `fit` object from edgeR or limma is used, the stats functions should
# convert this into a gene-per-row data-frame containing the two-tailed pvals
# and directions and pass this into the tibble statistics function

# As of 2018-03-06, we only combine genes by using metap::sumz approach

###############################################################################

#' Summarise whether the coexpression partners of each source-gene in the input
#' behave consistently
#'
#' For each source gene in the dataset, there is a set of target genes (the
#' neighbourhood of that source gene). A two-tailed p-value for both the
#' source and target genes should be defined (with respect to some experimental
#' comparison).
#'
#' For each source gene, this function combines the p-values observed in its
#' neighbourhood into a single summary score. We use the sum-of-z-scores method
#' as used in the `metap` package (although we have reimplemented it for
#' numerical stability); for a source gene with `n` neighbouring genes, a
#' z-score is computed for each neighbouring gene, and then we take sumz =
#' sum(z_i : i = 1 .. n) / sqrt(n) as the combined score. A two-tailed p-value
#' corresponding to `sumz` is also returned. (the z-score for a given
#' neighbour is obtained by comparing its p-value against the standard normal
#' distribution).
#'
#' @param       x              A data-frame containing columns \code{gene_id},
#'   \code{p_value} and \code{direction} (at least). The \code{p_value} column
#'   should contain two-tailed p-values. There should be no duplicate gene
#'   identifiers in the data-frame.
#'
#' @param       coex_partners   A subset of the coexpresDB.jp database
#'   containing the coexpression partners of a set of source-genes. As returned
#'   by \code{get_coex_partners}. Must contain columns \code{source_id} and
#'   \code{target_id}.
#'
#' @return      A data-frame containing a row for each source gene in the input
#'   and a summary of the number of partner-genes, the average z-score across
#'   all partner genes and the p-value equivalent to this z-score.
#'
#' @include      coxpresdbr_data_validity.R
#'
#' @importFrom   dplyr         group_by   mutate   n   summarise   ungroup   rename
#' @importFrom   magrittr      extract   %>%
#' @importFrom   methods       is
#' @importFrom   rlang         .data
#' @importFrom   stats         pchisq   qnorm
#' @importFrom   tibble        tibble
#'
#' @export
#'
evaluate_coex_partners <- function(
                                   x,
                                   coex_partners) {
  # An earlier implementation used `metap::sumz` to summarise over p-values.
  # We have decided against using that package: it required one-tailed p-values
  # and this led to numerical inaccuracies, eg, when an input one-tailed
  # p-value was > 1 - 1e-17, this was numerically rounded to 1 before
  # conversion to a z-score

  # We now do the same calculation as metap::sumz, but we construct z-scores
  # from the two-tailed p-values and direction-of-change

  # The sum-of-z-scores method for summarising p-values takes the form
  # sumz = sum(z_i : i = 1..n) / sqrt(n)
  # The resulting sumz value is converted back to a two-tailed p-value by
  # comparing to the standard normal distribution (or equivalently, it's square
  # is compared to ChiSquare with 1-df)

  if (!.is_gene_statistics_df(x)) {
    stop(
      "`x` should contain columns `gene_id`, `p_value`, `direction`",
      "and should have no duplicate gene-IDs in `evaluate_coex_partners`"
    )
  }

  stopifnot(methods::is(coex_partners, "data.frame"))
  stopifnot(all(c("source_id", "target_id") %in% colnames(coex_partners)))

  .add_z_scores <- function(.df) {
    .df %>%
      dplyr::mutate(
        z_score = qnorm(.data[["p_value"]] / 2) * (-1 * .data[["direction"]])
      )
  }

  .summarise_neighbours <- function(.df) {
    .df %>%
      dplyr::group_by(.data[["source_id"]]) %>%
      dplyr::summarise(
        n_partners = as.integer(n()),
        z_score = sum(.data[["z_score"]]) / sqrt(n()),
        p_value = pchisq(
          .data[["z_score"]] * .data[["z_score"]],
          df = 1, lower.tail = FALSE
        )
      ) %>%
      dplyr::ungroup()
  }

  res <- x %>%
    .add_z_scores() %>%
    merge(
      coex_partners,
      by.x = "gene_id", by.y = "target_id"
    ) %>%
    .summarise_neighbours() %>%
    dplyr::rename(gene_id = .data[["source_id"]]) %>%
    magrittr::extract(
      c("gene_id", "n_partners", "z_score", "p_value")
    )

  res
}

###############################################################################

#' .format_unsorted_nodes_for_tidygraph
#'
#' @importFrom   dplyr         rename
#' @importFrom   magrittr      %>%   extract
#' @importFrom   rlang         .data
#'
#' @noRd
#'
.format_unsorted_nodes_for_tidygraph <- function(
                                                 coex_partners) {
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@gene_statistics) > 0))

  if (!("z_score" %in% colnames(coex_partners@gene_statistics))) {
    coex_partners <- .add_z_scores(coex_partners)
  }

  coex_partners@gene_statistics %>%
    dplyr::rename(name = .data[["gene_id"]]) %>%
    magrittr::extract(c("name", "z_score", "p_value", "direction"))
}

###############################################################################

#' .format_coex_edges_for_tidygraph
#'
#' @importFrom   dplyr         rename
#' @importFrom   methods       is
#' @importFrom   rlang         .data
#'
#' @noRd
#'
.format_coex_edges_for_tidygraph <- function(
                                             coex_partners,
                                             cluster_source_nodes_only = TRUE) {
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@partners) > 0))

  relabelled <- dplyr::rename(
    coex_partners@partners,
    from = .data[["source_id"]], to = .data[["target_id"]]
  )

  if (cluster_source_nodes_only) {
    keep_rows <- relabelled[["to"]] %in% relabelled[["from"]]
    relabelled[keep_rows, ]
  } else {
    relabelled
  }
}

###############################################################################

#' .add_direction_parities_to_coex_edges
#'
#' @importFrom   dplyr      inner_join   mutate
#' @importFrom   magrittr   %>%   extract
#' @importFrom   methods    is
#' @importFrom   rlang      .data
#'
#' @noRd
#'
.add_direction_parities_to_coex_edges <- function(
                                                  coex_partners) {
  # both gene-statistics and partners should be defined in coex_partners
  # this modifies the contents of coex_partners@partners, adding the column
  # `direction_parity`
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@gene_statistics) > 0))
  stopifnot(all(dim(coex_partners@partners) > 0))

  gene_to_gene <- coex_partners@partners[
    c("source_id", "target_id")
  ]

  gene_sub_statistics <- coex_partners@gene_statistics[
    c("gene_id", "direction")
  ]

  direction_matches <- gene_to_gene %>%
    merge(y = gene_sub_statistics, by.x = "source_id", by.y = "gene_id") %>%
    merge(y = gene_sub_statistics, by.x = "target_id", by.y = "gene_id") %>%
    dplyr::mutate(
      direction_parity = .data[["direction.x"]] == .data[["direction.y"]]
    ) %>%
    magrittr::extract(c("source_id", "target_id", "direction_parity"))

  coex_partners@partners <- dplyr::inner_join(
    x = coex_partners@partners,
    y = direction_matches,
    by = c("source_id", "target_id")
  )

  coex_partners
}

###############################################################################

#' Obtain a graph containing clusters of genes, where known co-expression
#' partners are linked if their expression is consistent in the user's data
#'
#' @param        coex_partners   A \code{CoxpresDbPartners} object containing
#'   gene-statistics from the user's experiment and gene-gene associations
#'   between highly correlated partners as derived from `coxpresdb.jp`.
#'
#' @param        drop_disparities   Boolean. Should the clustering of known
#'   coexpression partners disregard associations between genes that are
#'   differentially expressed in the opposite direction from each other?
#'   Default: TRUE.
#'
#' @importFrom   dplyr         filter   left_join
#' @importFrom   igraph        get.vertex.attribute   vertex_attr
#' @importFrom   rlang         .data
#' @importFrom   tibble        tibble   as_tibble
#' @importFrom   tidygraph     as_tbl_graph
#'
#' @export

cluster_by_coex_partnership <- function(
                                        coex_partners,
                                        drop_disparities = TRUE) {
  stopifnot(
    is.logical(drop_disparities) && length(drop_disparities) == 1
  )

  .edge_filter <- if (drop_disparities) {
    function(x) dplyr::filter(x, .data[["direction_parity"]])
  } else {
    identity
  }

  node_attributes <- .format_unsorted_nodes_for_tidygraph(
    coex_partners
  )

  edges <- coex_partners %>%
    .add_direction_parities_to_coex_edges() %>%
    .format_coex_edges_for_tidygraph(
      cluster_source_nodes_only = TRUE
    ) %>%
    dplyr::filter(.data[["to"]] %in% .data[["from"]]) %>%
    .edge_filter()

  graph <- tidygraph::as_tbl_graph(
    edges
  )

  nodes <- if (nrow(edges) == 0) {
    tibble::tibble(name = character(0))
  } else {
    igraph::get.vertex.attribute(graph) %>%
      tibble::as_tibble() %>%
      dplyr::left_join(node_attributes, by = "name")
  }

  igraph::vertex_attr(graph) <- as.list(nodes)

  coex_partners@cluster_graph <- graph
  coex_partners
}

###############################################################################
