#' Extract variable name as string
#'
#' @param var_x A variable
#'
#' @return
#' @export
var_to_str <- function(var_x) {
  deparse(substitute(var_x))
}
