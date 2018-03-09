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

methods::setClass(
  "CoxpresDbImporter",
  slots = list(
    archive = "character",
    archive_uncompressed = "character",
    file_paths = "data.frame"
  )
)

###############################################################################
# TODO: write these notes up in the rox for CoxpresDbImporter

# - check validity of db_archive
#     - must be a file, may be .tar or .tar.bz2
#     - if it's a .tar.bz2
#         - temp_dir should be defined
#         - uncompressed archive is set to temp_dir ++ basename(db_archive)
#     - otherwise pass archive to archive_uncompressed

# - set up uncompressed archive
#     - call to system's "bunzip2 --keep --stdout {archive.tar.bz2} > \
#         {temp_dir/archive.tar}"
#     - and set archive_uncompressed to {temp_dir/archive.tar}

# - obtain all file-paths and associated genes from the uncompressed archive
#   and store it in a data.frame as file_paths

#' Constructor for the CoxpresDbImporter class - Use this to obtain data from
#' a CoxpresDB.jp archive
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
#' @importFrom   R.utils       isBzipped   bunzip2
#'
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

# Validate a CoxpresDbImporter
# - archive and archive_uncompressed should be files
# - archive_uncompressed should be a .tar
# - file_paths should have colnames "gene_id" and "file_path"

###############################################################################

# Accessors
# - get_gene_ids - done
# - get_file_paths - done
# - get_raw_archive - done
# - get_uncompressed_archive - done

# Methods
# - import_coex_partners(CoxpresDbImporter, gene_id)

#' @importFrom   methods       .valueClassTest
#'
setGeneric("get_file_paths", valueClass = "data.frame", function(x) {
  standardGeneric("get_file_paths")
})

setMethod("get_file_paths", signature("CoxpresDbImporter"), function(x) {
  x@file_paths
})

setGeneric(
  "get_file_path_for_gene", valueClass = "character",
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

setGeneric("get_gene_ids", valueClass = "character", function(x) {
  standardGeneric("get_gene_ids")
})

setMethod("get_gene_ids", signature("CoxpresDbImporter"), function(x) {
  get_file_paths(x)[["gene_id"]]
})

setGeneric("get_raw_archive", valueClass = "character", function(x) {
  standardGeneric("get_raw_archive")
})

setMethod("get_raw_archive", signature("CoxpresDbImporter"), function(x) {
  x@archive
})

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

#' Extracts the relative file-paths for each file within a given CoxpresDB.jp
#' archive
#'
#' @inheritParams   .is_coxpresdb_archive
#'
#' @importFrom   gtools        mixedsort
#' @importFrom   stringr       str_interp
#'
#' @return       A vector of file-paths, sorted alphanumerically, each file is
#' present in the given archive and represents the coexpression data for a
#' particular gene.
#'
#' @export
#'
.get_coxpresdb_file_paths <- function(
                                      db_archive) {
  # TODO: remove this function
  get_file_paths(
    CoxpresDbImporter(
      db_archive = db_archive, overwrite_in_bunzip2 = TRUE
    )
  )$file_path
}

###############################################################################

#' Import the coxpression dataset for a single gene
#'
#' @param        gene_id       A single identifier for a gene. This should be
#' present in the database (an error is thrown if not).
#'
#' @inheritParams   .is_coxpresdb_archive
#'
#' @return       A single dataframe containing the source -> target mappings,
#' and both the mutual ranks and correlations between the gene pairs.
#'
#' @importFrom   data.table    fread
#' @importFrom   dplyr         filter_   mutate_
#' @importFrom   magrittr      %>%   set_colnames
#' @importFrom   tibble        as_data_frame   data_frame
#'
#' @export
#'
import_coex_db <- function(
                           gene_id,
                           db_archive) {
  importer <- CoxpresDbImporter(db_archive, overwrite_in_bunzip2 = TRUE)

  # gene_files <- .get_coxpresdb_file_paths(db_archive)
  if (length(gene_id) != 1) {
    stop("`gene_id` should be a single gene-id in `import_coex_db`")
  }

  if (!(gene_id %in% get_gene_ids(importer))) {
    stop(
      sprintf(
        "Gene-ID %s is not present in the CoxpresDB archive %s",
        gene_id, db_archive
      )
    )
  }

  expected_colnames <- c(
    "source_id", "target_id", "mutual_rank", "correlation"
  )

  gene_file <- get_file_path_for_gene(gene_id, importer)

  stopifnot(length(gene_file) == 1)

  coex_db <- data.table::fread(
    paste("tar --to-stdout -xf", get_uncompressed_archive(importer), gene_file)
  ) %>%
    magrittr::set_colnames(expected_colnames[-1]) %>%
    tibble::as_data_frame() %>%
    dplyr::mutate_(
      source_id = ~ gene_id,
      target_id = ~ as.character(target_id)
    ) %>%
    magrittr::extract(expected_colnames) %>%
    dplyr::filter_(~ target_id != gene_id)

  coex_db
}

###############################################################################

#' Imports the gene identifiers that are represented within a given
#' CoxpresDB.jp archive
#'
#' @inheritParams   .is_coxpresdb_archive
#'
#' @return       A vector of gene-ids; data for each such gene is present in
#' the user-supplied archive
#'
#' @importFrom   gtools        mixedsort
#' @export
#'
import_coex_db_universe <- function(
                                    db_archive) {
  get_gene_ids(CoxpresDbImporter(db_archive, overwrite_in_bunzip2 = TRUE))
}

###############################################################################
