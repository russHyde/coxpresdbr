###############################################################################

# Helper functions for use while developing

#' Nonexported function
#'
#' @importFrom   styler        style_pkg
#' @importFrom   devtools      build   check   document
#'
#' @noRd

.packit <- function(pkg = ".") {
  styler::style_pkg(pkg)

  devtools::document(pkg)

  devtools::build(pkg)

  devtools::check(pkg)
}

###############################################################################
