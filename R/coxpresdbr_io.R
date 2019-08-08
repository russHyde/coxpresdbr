###############################################################################
# import functions for `coxpresdbr` package
#
###############################################################################

###############################################################################

#' Checks if a file has a .zip extension
#'
#' @noRd
.is_zip_file <- function(x) {
  grepl(pattern = ".*\\.zip$", x = x)
}

#' Checks if a file has a .tar or .tar.bz2 extension
#'
#' @noRd
.is_tar_file <- function(x) {
  grepl(pattern = ".*\\.tar(\\.bz2){0,1}", x)
}

#' Checks whether a provided file-path points to a valid CoxpresDB.jp archive
#'
#' @param        db_archive    A file path. Should be a single
#'   \code{`*`.tar.bz2}, \code{`*`.tar}, or \code{`*`.zip} file as downloaded
#'   from the _coxpresdb.jp_ website.
#'
.is_coxpresdb_archive <- function(
                                  db_archive) {
  length(db_archive) == 1 &&
    file.exists(db_archive) &&
    (
      .is_zip_file(db_archive) ||
        .is_tar_file(db_archive)
    )
}

###############################################################################

#' Constructor for the CoxpresDbAccessor class - Use this to obtain data from
#' a CoxpresDB.jp archive
#'
#' The \code{db_archive} should be a \code{.zip}, \code{.tar} or a
#' \code{.tar.bz2}, if it's a \code{.tar.bz2} then \code{temp_dir} should be
#' defined and a random copy of the uncompressed archive will be made. All
#' access to the stored data will be made via the uncompressed copy of the
#' archive, so make a \code{CoxpresDbAccessor} _once_ during any given script.
#'
#' @inheritParams   .is_coxpresdb_archive
#'
#' @param        temp_dir      A directory into which a compressed CoxpresDB
#'   archive will be decompressed. By default this is the temp_dir for the
#'   current R session. Only relevant if the coxpresdb archive is compressed.
#'
#' @param        overwrite_in_bunzip2   Boolean. If the CoxpresDB archive is
#'   compressed, and a decompressed copy of the archive is found in the target
#'   directory, should the function throw an exception? See \code{overwrite} in
#'   \code{R.utils::bunzip2}
#'
#' @param        remove_in_bunzip2   Boolean. If a compressed CoxpresDB archive
#'   is provided, should the compressed version be deleted after decompression?
#'   See \code{remove} in \code{R.utils::bunzip2}
#'
#' @importFrom   methods       new
#' @importFrom   R.utils       isBzipped   bunzip2
#' @importFrom   tibble        tibble
#' @importFrom   utils         untar
#'
#' @include      coxpresdbr_CoxpresDbAccessor.R
#'
#' @return       A \code{CoxpresDbAccessor} object.
#' @export
#'
CoxpresDbAccessor <- function(
                              db_archive,
                              temp_dir = tempdir(),
                              overwrite_in_bunzip2 = FALSE,
                              remove_in_bunzip2 = FALSE) {
  # Make an object for in-memory access to the coxpresDb data
  if (is.data.frame(db_archive)) {
    return(
      new("CoxpresDbDataframeAccessor", df = db_archive)
    )
  }

  # Make an object for from-file access to the coxpresDb adata
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

  file_listing_function <- if (.is_zip_file(db_archive)) {
    function(x) utils::unzip(x, list = TRUE)[["Name"]]
  } else {
    function(x) utils::untar(x, list = TRUE)
  }

  gene_files <- sort(
    grep(
      pattern = ".*/$",
      file_listing_function(db_uncompressed),
      value = TRUE,
      invert = TRUE
    )
  )

  gene_file_df <- tibble::tibble(
    gene_id = basename(gene_files),
    file_path = gene_files
  )

  methods::new(
    "CoxpresDbArchiveAccessor",
    archive = db_archive,
    archive_uncompressed = db_uncompressed,
    file_paths = gene_file_df
  )
}

###############################################################################

#' Obtains all the file-paths that are present in a CoxpresDb archive
#'
#' @param        x             A \code{CoxpresDbArchiveAccessor} object
#'   corresponding to a CoxpresDb.jp archive
#'
#' @return       A tibble of gene_id:file_path pairs, each file is present in
#'   the archive (referenced within \code{CoxpresDbArchiveAccessor}) and
#'   represents the coexpression data for a particular gene.
#'
#' @importFrom   methods       .valueClassTest
#'
setGeneric("get_file_paths", valueClass = "data.frame", function(x) {
  standardGeneric("get_file_paths")
})

setMethod("get_file_paths", signature("CoxpresDbArchiveAccessor"), function(x) {
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
  signature(gene_id = "character", importer = "CoxpresDbArchiveAccessor"),
  function(gene_id, importer) {
    paths <- get_file_paths(importer)
    rows <- which(paths$gene_id %in% gene_id)
    paths[["file_path"]][rows]
  }
)

###############################################################################

#' Generic method for getting gene identifiers.
#'
#' The returned vector contains all gene IDs present in the dataset: this
#' should include all source- and target- genes.
#'
#' @param        x             A \code{CoexpresDbAccessor} corresponding to a
#'   CoxpresDB archive.
#'
setGeneric("get_gene_ids", valueClass = "character", function(x) {
  standardGeneric("get_gene_ids")
})

#' Imports the gene identifiers that are represented within a given
#' CoxpresDB.jp archive
#'
#' @param        x             A \code{CoexpresDbArchiveAccessor} corresponding
#'   to a CoxpresDB archive.
#'
#' @return       A vector of gene-ids; data for each such gene is present in
#'   the user-supplied archive.
#'
#' @export
#'

setMethod("get_gene_ids", signature("CoxpresDbArchiveAccessor"), function(x) {
  # Note that without opening every file in the archive, there is no way of
  # knowing whether there are some genes mentioned inside one of those files
  # which is not a source gene (for an archive, the source genes are the
  # filenames present in the archive)
  get_source_ids(x)
})

#' Imports the gene identifiers that are represented within a data-frame-based
#' CoxpresDB.jp archive
#'
#' @param        x             A \code{CoxpresDbDataframeAccessor}
#'   corresponding to a CoxpresDB archive.
#'
#' @return       A vector of gene-ids; data for each such gene is present in
#'   the user-supplied archive.
#'
#' @export
#'

setMethod("get_gene_ids", signature("CoxpresDbDataframeAccessor"), function(x) {
  sort(union(get_source_ids(x), x@df[["target_id"]]))
})

###############################################################################

#' Generic method for getting source-gene identifiers.
#'
#' The returned vector contains all gene IDs for which coexpression partners
#' can be obtained from the dataset. There may be target-genes that are not
#' source-genes.
#'
#' @param        x             A \code{CoexpresDbAccessor} corresponding to a
#'   CoxpresDB archive.
#'
setGeneric("get_source_ids", valueClass = "character", function(x) {
  standardGeneric("get_source_ids")
})

#' Imports the source-gene identifiers that are represented within a given
#' CoxpresDB.jp archive
#'
#' @param        x             A \code{CoexpresDbArchiveAccessor} corresponding
#'   to a CoxpresDB archive.
#'
#' @return       A vector of gene-ids; coexpression partners for each such gene
#'   are available in the user-supplied archive.
#'
#' @export
#'

setMethod("get_source_ids", signature("CoxpresDbArchiveAccessor"), function(x) {
  get_file_paths(x)[["gene_id"]]
})

#' Imports the source gene identifiers that are represented within a
#' data-frame-based CoxpresDB.jp archive
#'
#' @param        x             A \code{CoxpresDbDataframeAccessor}
#'   corresponding to a CoxpresDB archive.
#'
#' @return       A vector of gene-ids; coexpression partners for each such gene
#'   are present in the user-supplied archive.
#'
#' @export
#'

setMethod(
  "get_source_ids",
  signature("CoxpresDbDataframeAccessor"),
  function(x) {
    sort(unique(x@df[["source_id"]]))
  }
)

# TODO: `get_target_ids`

###############################################################################

setGeneric("get_raw_archive", valueClass = "character", function(x) {
  standardGeneric("get_raw_archive")
})

setMethod(
  "get_raw_archive", signature("CoxpresDbArchiveAccessor"),
  function(x) {
    x@archive
  }
)

###############################################################################

setGeneric("get_uncompressed_archive", valueClass = "character", function(x) {
  standardGeneric("get_uncompressed_archive")
})

setMethod(
  "get_uncompressed_archive", signature("CoxpresDbArchiveAccessor"),
  function(x) {
    x@archive_uncompressed
  }
)

###############################################################################

setGeneric(
  "get_all_coex_partners",
  valueClass = "data.frame",
  function(gene_id, importer) {
    standardGeneric("get_all_coex_partners")
  }
)

#' Import all the coexpression partner data for a single gene from a given
#' CoxpresDb archive.
#'
#' Users should use \code{get_coex_partners(gene_ids, importer)} rather than
#' this (unexported) method.
#'
#' @param        gene_id       A gene-identifier in the same format as present
#'   throughout the CoxpresDb archive. If this isn't a single identifier, or
#'   it isn't present in the database, the function will throw an error.
#'
#' @param        importer      \code{CoxpresDbArchiveAccessor} object.
#'
#' @return       A single dataframe containing the source -> target mappings,
#'   and the mutual ranks between the gene pairs.
#'
#' @importFrom   data.table    fread
#' @importFrom   dplyr         arrange
#' @importFrom   tibble        tibble
#'
setMethod(
  "get_all_coex_partners",
  signature("character", "CoxpresDbArchiveAccessor"),
  function(gene_id, importer) {
    # This function is not exported.

    # We assume that any function which calls this has passed in a single gene
    # and has already checked that the gene_id passed in is a valid identifier
    # for this coexpression dataset

    # TODO: implement for multiple 'gene_id's

    .drop_self_edges <- function(df) {
      # returned data-frames should contain a single gene in the source column,
      # and that gene should be absent from the target column, and the rows
      # should be ordered by increasing mutual-rank value.
      rows <- which(df[["source_id"]] != df[["target_id"]])

      df[rows, ]
    }

    gene_file <- get_file_path_for_gene(gene_id, importer)

    if (length(gene_file) != 1) {
      stop(
        "Either > 1 gene_id passed to `get_all_coex_partners`,",
        "or no file was found in the archive for this gene,",
        "or > 1 file was found in the archive for this gene"
      )
    }

    archive <- get_uncompressed_archive(importer)

    import_command <- if (.is_tar_file(archive)) {
      paste("tar --to-stdout -xf", archive, gene_file)
    } else if (.is_zip_file(archive)) {
      paste("unzip -p", archive, gene_file)
    } else {
      stop("archive should be either .zip or .tar in `get_all_coex_partners`")
    }

    initial_db <- data.table::fread(cmd = import_command)

    # Original versions of coxpresdb files were of the format (target_id,
    # mutual_rank, correlation), with each seed gene having a separate file
    #
    # In 2019, versions of the coxpresdb files do not contain a correlation
    # coefficient column (ie, target_id, mutual_rank)
    #
    # We disregard the correlation-coefficient column since it isn't
    # consistently presented across different releases of coxpresdb

    tibble::tibble(
      source_id = gene_id,
      target_id = as.character(initial_db[[1]]),
      mutual_rank = initial_db[[2]]
    ) %>%
      .drop_self_edges()
  }
)

###############################################################################
