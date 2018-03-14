###############################################################################

#' Template for CoxpresDbImporter class
#'
#' This is a datastructure that stores the location of a CoxpresDb archive and
#' the relative paths of all subfiles within that archive. If the archive is
#' compressed, on construction of CoxpresDbImporter a temporary copy of the
#' uncompressed archive is constructed. So construction may take a couple of
#' minutes depending on the size of the archive. You may want to wrap a call to
#' this function in a \code{future()} block.
#'
#' @param        archive       The path to a CoxpresDb archive (either .tar.bz2
#'  or .tar).
#' @param        archive_uncompressed   The path to an uncompressed copy of the
#' CoxpresDb archive.
#' @param        file_paths    The relative paths of all file present in
#' subdirectories of the archive. As a data-frame of gene_id -> file_path
#' pairs.
#'
#' @name         CoxpresDbImporter-class
#' @rdname       CoxpresDbImporter-class
#'
#' @export       CoxpresDbImporter
#'
methods::setClass(
  "CoxpresDbImporter",
  slots = list(
    archive = "character",
    archive_uncompressed = "character",
    file_paths = "data.frame"
  )
)

###############################################################################

#'
.validity_coxpresdb_partners <- function(object) {
  # S3 tests:
  if (!is.null(object@cluster_graph) &&
    !is(object@cluster_graph, "igraph")
  ) {
    return(
      paste(
        "`cluster_graph` should be NULL or inherit from `igraph` in",
        "`CoxpresDbPartners`"
      )
    )
  }

  # TODO: data.frame tests
}

###############################################################################

#' Template for CoxpresDbPartners class
#'
methods::setClass(
  "CoxpresDbPartners",
  slots = list(
    gene_statistics = "data.frame",
    partners = "data.frame",
    partner_summaries = "data.frame",
    cluster_graph = "ANY"
  ),
  validity = .validity_coxpresdb_partners
)

###############################################################################

setGeneric("get_gene_universe", valueClass = "character", function(x) {
  standardGeneric("get_gene_universe")
})

setMethod("get_gene_universe", signature("CoxpresDbPartners"), function(x) {
  x@gene_statistics$gene_id
})

###############################################################################
