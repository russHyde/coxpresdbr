# Helper functions for use while developing

#' Nonexported function
#'
#' @importFrom   styler        style_pkg
#' @importFrom   lintr         lint_package
#' @importFrom   devtools      build   check
#'
.packit <- function(pkg = ".") {
  styler::style_pkg(pkg)

  # lintr doesn't push it's results into rstudio when called inside .packit()
  # lintr::lint_package(pkg)

  devtools::build(pkg)

  devtools::check(pkg)
}
