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
