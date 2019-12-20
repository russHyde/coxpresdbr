#' Compute the (partial) correlation between columns of the matrix `x`.
#'
#' If the matrix `control` is defined, the partial correlation is computed
#' with-respect-to this matrix.
#'
#' From the outside, this is vectorised over source_id and target_id; each
#' source_id / target_id pair will be used in the correlation analysis.
#'
#' For each source_id, s:
#' - For each target_id, t, of s:
#'   - Compute the (partial) correlation between s and t
#'
#' [For partial correlation, compute the correlation between the residuals of
#' s and the residuals of t, after regressing-out the `control` matrix]
#'
#' Then compute the z-score for each correlation coefficient based on Fisher's
#' transformation (using Fisher's z-score transformation for correlation
#' coefficients)
#'
#' @param   x   A matrix. We want to perform correlation between pairs of
#'   columns of \code{x}.
#' @param   source_id,target_id   Colnames within the matrix `x`; each pair of
#'   (source_id, target_id) values corresponds to a pair of column vectors in
#'   `x` that we want to perform correlation analysis on.
#' @param   control A matrix. this provides a numeric control for partial
#'   correlation (eg, allowing a different baseline mean for different
#'   datasets). See ppcor::pcor for details of how to specify control.
#'
#' @return       A tibble (source_id, target_id, partial_cor, z_score)
#' @export

compute_cor <- function(x, source_id, target_id, control = NULL) {
  stopifnot(is.matrix(x) && is.numeric(x))
  stopifnot(is.character(source_id) && all(source_id %in% colnames(x)))
  stopifnot(is.character(target_id) && all(target_id %in% colnames(x)))
  stopifnot(is.null(control) || is.numeric(control))

  if (any(source_id == target_id)) {
    stop("source and target id should be distinct in each pair")
  }

  # If no 'control' variables are present, we extract the standard correlation
  # coefficient; if a control matrix is present, we perform partial-correlation
  # including the explanatory control
  #
  cor_fn <- if (is.null(control)) {
    function(x, y) stats::cor(x, y)
  } else {
    function(x, y) ppcor::pcor(cbind(x, y, control))[["estimate"]][1, 2]
  }

  pcors <- purrr::map2_dbl(
    source_id, target_id, function(s0, t0) {
      cor_fn(x[, s0], x[, t0])
    }
  )

  tibble::tibble(
    source_id = source_id,
    target_id = target_id,
    partial_cor = pcors,
    # Fisher's z-transformation of the correlation coefficients
    z_score = 0.5 * log((1 + pcors) / (1 - pcors))
  )
}
