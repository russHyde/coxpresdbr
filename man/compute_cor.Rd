% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_cor.R
\name{compute_cor}
\alias{compute_cor}
\title{Compute the (partial) correlation between columns of the matrix `x`.}
\usage{
compute_cor(x, source_id, target_id, control = NULL)
}
\arguments{
\item{x}{A matrix. We want to perform correlation between pairs of
columns of \code{x}.}

\item{source_id, target_id}{Colnames within the matrix `x`; each pair of
(source_id, target_id) values corresponds to a pair of column vectors in
`x` that we want to perform correlation analysis on.}

\item{control}{A matrix. this provides a numeric control for partial
correlation (eg, allowing a different baseline mean for different
datasets). See ppcor::pcor for details of how to specify control.}
}
\value{
A tibble (source_id, target_id, partial_cor, z_score)
}
\description{
If the matrix `control` is defined, the partial correlation is computed
with-respect-to this matrix.
}
\details{
From the outside, this is vectorised over source_id and target_id; each
source_id / target_id pair will be used in the correlation analysis.

For each source_id, s:
- For each target_id, t, of s:
  - Compute the (partial) correlation between s and t

[For partial correlation, compute the correlation between the residuals of
s and the residuals of t, after regressing-out the `control` matrix]

Then compute the z-score for each correlation coefficient based on Fisher's
transformation (using Fisher's z-score transformation for correlation
coefficients)
}
