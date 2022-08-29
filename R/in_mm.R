#' Convert inches to mm
#'
#' @param inches
#'
#' @return
#' @export
#'
#' @examples
in_mm <- function(inches) {
  grid::convertUnit(unit(inches, "in"), "mm", valueOnly = TRUE)
}
