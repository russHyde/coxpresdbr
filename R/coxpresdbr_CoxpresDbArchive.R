###############################################################################

#' Template for \code{CoxpresDbAccessor} class
#'
#' `CoxpresDbAccessor` is a virtual class, to construct an object of this
#' class, use `CoxpresDbArchiveAccessor` or `CoxpresDbDataframeAccessor`.
#'
#' @name         CoxpresDbAccessor-class
#' @rdname       CoxpresDbAccessor-class
#'
#' @exportClass       CoxpresDbAccessor
#'
methods::setClass(
  "CoxpresDbAccessor"
)

#' Constructor for the concrete class `CoxpresDbArchiveAccessor`
#'
#' This is a datastructure that stores the location of a CoxpresDb archive and
#' the relative paths of all subfiles within that archive. If the archive is
#' compressed, on construction of \code{CoxpresDbAccessor} a temporary copy of
#' the uncompressed archive is constructed. So construction may take a couple
#' of minutes depending on the size of the archive. You may want to wrap a call
#' to this function in a \code{future()} block.
#'
#' @param        archive       The path to a CoxpresDb archive (either .tar.bz2
#'   or .tar).
#' @param        archive_uncompressed   The path to an uncompressed copy of the
#'   CoxpresDb archive.
#' @param        file_paths    The relative paths of all file present in
#'   subdirectories of the archive. As a data-frame of `gene_id` -> `file_path`
#'   pairs.
#'
#' @name         CoxpresDbArchiveAccessor-class
#' @rdname       CoxpresDbArchiveAccessor-class
#'
methods::setClass(
  "CoxpresDbArchiveAccessor",
  slots = list(
    archive = "character",
    archive_uncompressed = "character",
    file_paths = "data.frame"
  ),
  contains = "CoxpresDbAccessor"
)

###############################################################################
