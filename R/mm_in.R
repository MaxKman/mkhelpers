#' Convert mm to inches
#'
#' @param mm
#'
#' @return
#' @export
#'
#' @examples
mm_in <- function(mm) {
  grid::convertUnit(unit(mm, "mm"), "in", valueOnly = TRUE)
}
