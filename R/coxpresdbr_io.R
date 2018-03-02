###############################################################################
# import functions for `coxpresdbr` package
#
###############################################################################

# data.table::fread("tar --to-stdout -xjf {archive_name} {file_name}")
#
#' @export
#'
import_coex_db <- function(gene_id, db_archive) {
  stopifnot(file.exists(db_archive))
}

# system("tar --list -jf {archive_name} | grep -v '/$'", intern = TRUE)
#
get_coex_db_universe <- function(db_archive) {
}
