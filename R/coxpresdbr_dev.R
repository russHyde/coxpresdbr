# Helper functions for use while developing

#' Nonexported function
#'
#' @importFrom   styler        style_pkg
#' @importFrom   lintr         lint_package
#' @importFrom   devtools      build   check
#'
.packit <- function(pkg = ".") {
  styler::style_pkg(pkg)

  devtools::document(pkg)

  devtools::build(pkg)

  devtools::check(pkg)
}
