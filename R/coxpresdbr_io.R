###############################################################################
# import functions for `coxpresdbr` package
#
###############################################################################

###############################################################################

#' Checks whether a provided file-path points to a valid CoxpresDB.jp archive
#'
#' @param        db_archive    A file path. Should be a single \code{*.tar.bz2}
#' or \code{*.tar} file as downloaded from the coxpresdb.jp website.
#'
#' @importFrom   stringr       str_detect
#'
.is_coxpresdb_archive <- function(
                                  db_archive) {
  length(db_archive) == 1 &&
    file.exists(db_archive) &&
    stringr::str_detect(db_archive, "\\.tar(.bz2)*$")
}


###############################################################################

#' Constructor for the CoxpresDbImporter class - Use this to obtain data from
#' a CoxpresDB.jp archive
#'
#' The db_archive should be a .tar or a .tar.bz2, if it's a .tar.bz2 then
#' \code{temp_dir} should be defined and a random copy of the uncompressed
#' archive will be made. All access to the stored data will be made via the
#' uncompressed copy of the archive, so make a CoxpresDbImporter _once_ during
#' any given script.
#'
#' @inheritParams   .is_coxpresdb_archive
#'
#' @param        temp_dir      A directory into which a compressed CoxpresDB
#' archive will be decompressed. By default this is the temp_dir for the
#' current R session. Only relevant if the coxpresdb archive is compressed.
#'
#' @param        overwrite_in_bunzip2   Boolean. If the CoxpresDB archive is
#' compressed, and a decompressed copy of the archive is found in the target
#' directory, should the function throw an exception? See \code{overwrite} in
#' \code{R.utils::bunzip2}
#'
#' @param        remove_in_bunzip2   Boolean. If a compressed CoxpresDB archive
#' is provided, should the compressed version be deleted after decompression?
#' See \code{remove} in \code{R.utils::bunzip2}
#'
#' @importFrom   methods       new
#' @importFrom   R.utils       isBzipped   bunzip2
#' @importFrom   tibble        data_frame
#' @importFrom   utils         untar
#'
#' @include      coxpresdbr_classes.R
#'
#' @return       A CoxpresDbImporter object.
#' @export
#'
CoxpresDbImporter <- function(
                              db_archive,
                              temp_dir = tempdir(),
                              overwrite_in_bunzip2 = FALSE,
                              remove_in_bunzip2 = FALSE) {
  stopifnot(.is_coxpresdb_archive(db_archive))

  db_uncompressed <- if (R.utils::isBzipped(db_archive)) {
    stopifnot(dir.exists(temp_dir))

    out_file <- file.path(
      temp_dir,
      gsub("[.]bz2$", "", basename(db_archive))
    )

    R.utils::bunzip2(
      filename = db_archive,
      destname = out_file,
      remove = remove_in_bunzip2,
      overwrite = overwrite_in_bunzip2
    )

    out_file
  } else {
    db_archive
  }

  gene_files <- sort(grep(
    pattern = ".*/$",
    utils::untar(db_uncompressed, list = TRUE),
    value = TRUE,
    invert = TRUE
  ))
  gene_file_df <- tibble::data_frame(
    gene_id = basename(gene_files),
    file_path = gene_files
  )

  methods::new(
    "CoxpresDbImporter",
    archive = db_archive,
    archive_uncompressed = db_uncompressed,
    file_paths = gene_file_df
  )
}

###############################################################################

#' Obtains all the file-paths that are present in a CoxpresDb archive
#'
#' @param        x             A \code{CoxpresDbImporter} object corresponding
#' to a CoxpresDb.jp archive
#'
#' @return       A data_frame of gene_id:file_path pairs, each file is present
#' in the archive (referenced within CoxpresDbImporter) and represents the
#' coexpression data for a particular gene.
#'
#' @importFrom   methods       .valueClassTest
#'
setGeneric("get_file_paths", valueClass = "data.frame", function(x) {
  standardGeneric("get_file_paths")
})

setMethod("get_file_paths", signature("CoxpresDbImporter"), function(x) {
  x@file_paths
})

###############################################################################

setGeneric(
  "get_file_path_for_gene",
  valueClass = "character",
  function(gene_id, importer) {
    standardGeneric("get_file_path_for_gene")
  }
)

setMethod(
  "get_file_path_for_gene",
  signature(gene_id = "character", importer = "CoxpresDbImporter"),
  function(gene_id, importer) {
    paths <- get_file_paths(importer)
    rows <- which(paths$gene_id %in% gene_id)
    paths[["file_path"]][rows]
  }
)

###############################################################################

#' Generic method for getting gene identifiers.
#'
#' @param        x             A \code{CoexpresDbImporter} corresponding to a
#' CoxpresDB archive.
#'
setGeneric("get_gene_ids", valueClass = "character", function(x) {
  standardGeneric("get_gene_ids")
})

#' Imports the gene identifiers that are represented within a given
#' CoxpresDB.jp archive
#'
#' @param        x             A \code{CoexpresDbImporter} corresponding to a
#' CoxpresDB archive.
#'
#' @return       A vector of gene-ids; data for each such gene is present in
#' the user-supplied archive
#'
#' @export
#'

setMethod("get_gene_ids", signature("CoxpresDbImporter"), function(x) {
  get_file_paths(x)[["gene_id"]]
})

###############################################################################

setGeneric("get_raw_archive", valueClass = "character", function(x) {
  standardGeneric("get_raw_archive")
})

setMethod("get_raw_archive", signature("CoxpresDbImporter"), function(x) {
  x@archive
})

###############################################################################

setGeneric("get_uncompressed_archive", valueClass = "character", function(x) {
  standardGeneric("get_uncompressed_archive")
})

setMethod(
  "get_uncompressed_archive",
  signature("CoxpresDbImporter"),
  function(x) {
    x@archive_uncompressed
  }
)

###############################################################################

setGeneric(
  "import_all_coex_partners",
  valueClass = "data.frame",
  function(gene_id, importer) {
    standardGeneric("import_all_coex_partners")
  }
)

#' Import all the coexpression partner data for a single gene from a given
#' CoxpresDb archive.
#'
#' Users should use get_coex_partners(gene_ids, importer) rather than this
#' method.
#'
#' @param        gene_id       A gene-identifier in the same format as present
#' throughout the CoxpresDb archive.
#' @param        importer      CoxpresDbImporter object.
#'
#' @return       A single dataframe containing the source -> target mappings,
#' and both the mutual ranks and correlations between the gene pairs.
#'
#' @importFrom   data.table    fread
#' @importFrom   dplyr         filter_   mutate_
#' @importFrom   magrittr      %>%   set_colnames
#' @importFrom   tibble        as_data_frame   data_frame
#'
setMethod(
  "import_all_coex_partners",
  signature("character", "CoxpresDbImporter"),
  function(gene_id, importer) {
    if (length(gene_id) != 1) {
      stop(
        "`gene_id` should be a single gene-id in `import_all_coex_partners`"
      )
    }

    if (!(gene_id %in% get_gene_ids(importer))) {
      stop(
        sprintf(
          "Gene-ID %s is not present in the CoxpresDB archive %s",
          gene_id, get_raw_archive(importer)
        )
      )
    }

    gene_file <- get_file_path_for_gene(gene_id, importer)

    stopifnot(length(gene_file) == 1)

    expected_colnames <- c(
      "source_id", "target_id", "mutual_rank", "correlation"
    )

    coex_db <- data.table::fread(
      paste(
        "tar --to-stdout -xf", get_uncompressed_archive(importer), gene_file
      )
    ) %>%
      magrittr::set_colnames(value = expected_colnames[-1]) %>%
      tibble::as_data_frame() %>%
      dplyr::mutate_(
        source_id = ~gene_id,
        target_id = ~ as.character(target_id)
      ) %>%
      magrittr::extract(expected_colnames) %>%
      dplyr::filter_(~ target_id != gene_id)

    coex_db
  }
)

###############################################################################
