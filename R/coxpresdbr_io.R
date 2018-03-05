###############################################################################
# import functions for `coxpresdbr` package
#
###############################################################################

#' Checks whether a provided file-path points to a valid CoxpresDB.jp archive
#'
#' @param        db_archive    A file path. Should be a single \code{*.tar.bz2}
#' file as downloaded from the coxpresdb.jp website.
#'
#' @importFrom   stringr       str_detect
#'
.is_coxpresdb_archive <- function(
                                  db_archive) {
  length(db_archive) == 1 &&
    file.exists(db_archive) &&
    stringr::str_detect(db_archive, "tar.bz2$")
}

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
  stopifnot(.is_coxpresdb_archive(db_archive))
  files <- system(
    # nolint start
    stringr::str_interp("tar --list -jf ${db_archive} | grep -v '/$'"),
    # nolint end
    intern = TRUE
  )
  gtools::mixedsort(files)
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
  gene_files <- .get_coxpresdb_file_paths(db_archive)
  stopifnot(
    length(gene_id) == 1 &&
      gene_id %in% basename(gene_files)
  )

  expected_colnames <- c(
    "source_id", "target_id", "mutual_rank", "correlation"
  )

  gene_file <- gene_files[basename(gene_files) == gene_id]
  stopifnot(length(gene_file) == 1)

  coex_db <- data.table::fread(
    paste("tar --to-stdout -xjf", db_archive, gene_file)
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
  gene_files <- .get_coxpresdb_file_paths(db_archive)
  basename(gene_files)
}

###############################################################################
