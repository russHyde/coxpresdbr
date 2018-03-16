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
#' @return      A data-frame containing a row for each source gene in the input
#' and a summary of the number of partner-genes, the average z-score across
#' all partner genes and the p-value equivalent to this z-score.
#'
#' @include      coxpresdbr_data_validity.R
#'
#' @importFrom   dplyr         group_by_   mutate_   n   summarise_   ungroup
#' @importFrom   magrittr      extract
#' @importFrom   metap         sumz   two2one
#' @importFrom   methods       is
#' @importFrom   tibble        data_frame
#'
#' @export
#'
evaluate_coex_partners <- function(
                                   x,
                                   coex_partners) {
  stopifnot(.is_gene_statistics_df(x))

  stopifnot(methods::is(coex_partners, "data.frame"))
  stopifnot(all(c("source_id", "target_id") %in% colnames(coex_partners)))

  res <- x %>%
    dplyr::mutate_(
      invert = ~ ifelse(direction < 0, TRUE, FALSE),
      p_value_onetail_forward = ~ metap::two2one(p_value, invert = invert),
      p_value_onetail_reversed = ~ metap::two2one(p_value, invert = !invert)
    ) %>%
    merge(
      coex_partners, by.x = "gene_id", by.y = "target_id"
    ) %>%
    dplyr::group_by_(~ source_id) %>%
    dplyr::summarise_(
      n_partners = ~ n(),
      z_score_forward = ~ ifelse(
        n_partners > 1,
        metap::sumz(p_value_onetail_forward)$z,
        qnorm(p_value_onetail_forward, lower.tail = FALSE)
      ),
      z_score_reversed = ~ ifelse(
        n_partners > 1,
        metap::sumz(p_value_onetail_reversed)$z,
        qnorm(p_value_onetail_reversed, lower.tail = FALSE)
      ),
      z_score = ~ (z_score_forward - z_score_reversed) / 2,
      p_value_onesided = ~ pnorm(z_score, lower.tail = FALSE),
      p_value = ~ ifelse(
        z_score > 0,
        2 * p_value_onesided,
        2 * (1 - p_value_onesided)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_(
      gene_id = ~ source_id,
      n_partners = ~ as.integer(n_partners),
      z_score = ~ as.numeric(z_score),
      p_value = ~ as.numeric(p_value)
    ) %>%
    magrittr::extract(
      c("gene_id", "n_partners", "z_score", "p_value")
    )

  res
}

###############################################################################

#' @importFrom   dplyr         rename_   select_
#' @importFrom   magrittr      %>%
#'
.format_unsorted_nodes_for_tidygraph <- function(
                                                 coex_partners) {
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@gene_statistics) > 0))
  stopifnot("z" %in% colnames(coex_partners@gene_statistics))

  coex_partners@gene_statistics %>%
    dplyr::rename_(.dots = list(name = "gene_id")) %>%
    dplyr::select_(.dots = c("name", "z", "p_value", "direction"))
}

###############################################################################

#' @importFrom   dplyr         filter_   transmute_
#' @importFrom   magrittr      %>%
#'
.format_coex_edges_for_tidygraph <- function(
                                             coex_partners,
                                             cluster_source_nodes_only = TRUE) {
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@partners) > 0))

  relabelled <- coex_partners@partners %>%
    dplyr::rename_(from = ~ source_id, to = ~ target_id)

  if (cluster_source_nodes_only) {
    dplyr::filter_(relabelled, ~ to %in% from)
  } else {
    relabelled
  }
}

###############################################################################

.add_direction_parities_to_coex_edges <- function(
                                                  coex_partners) {
  # both gene-statistics and partners should be defined in coex_partners
  # this modifies the contents of coex_partners@partners, adding the column
  # `direction_parity`
  stopifnot(methods::is(coex_partners, "CoxpresDbPartners"))
  stopifnot(all(dim(coex_partners@gene_statistics) > 0))
  stopifnot(all(dim(coex_partners@partners) > 0))

  gene_to_gene <- coex_partners@partners[c("source_id", "target_id")]

  gene_sub_statistics <- extract(
    coex_partners@gene_statistics, c("gene_id", "direction")
  )

  direction_matches <- gene_to_gene %>%
    merge(y = gene_sub_statistics, by.x = "source_id", by.y = "gene_id") %>%
    merge(y = gene_sub_statistics, by.x = "target_id", by.y = "gene_id") %>%
    dplyr::mutate_(direction_parity = ~ direction.x == direction.y) %>%
    dplyr::select_(.dots = c("source_id", "target_id", "direction_parity"))

  coex_partners@partners <- dplyr::inner_join(
    x = coex_partners@partners,
    y = direction_matches,
    by = c("source_id", "target_id")
  )

  coex_partners
}

###############################################################################

#' @importFrom   dplyr         filter_   left_join
#' @importFrom   igraph        get.vertex.attribute   vertex_attr
#' @importFrom   tibble        as_data_frame
#' @importFrom   tidygraph     as_tbl_graph
cluster_by_coex_partnership <- function(
                                        coex_partners,
                                        drop_disparities = TRUE) {
  stopifnot(is.logical(drop_disparities))

  node_attributes <- .format_unsorted_nodes_for_tidygraph(
    coex_partners
  )

  edges <- coex_partners %>%
    .add_direction_parities_to_coex_edges() %>%
    .format_coex_edges_for_tidygraph(cluster_source_nodes_only = TRUE) %>%
    dplyr::filter_(~ to %in% from) %>%
    dplyr::filter_(~ direction_parity)

  graph <- tidygraph::as_tbl_graph(
    edges
  )

  nodes <- igraph::get.vertex.attribute(graph) %>%
    tibble::as_data_frame() %>%
    dplyr::left_join(node_attributes, by = "name")

  igraph::vertex_attr(graph) <- as.list(nodes)

  coex_partners@cluster_graph <- graph
  coex_partners
}

###############################################################################
