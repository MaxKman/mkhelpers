#' Remove NAs from vector or list
#'
#' @param A vector or list
#'
#' @return A vector or list with NA entries removed
#' @export
#'
#' @examples
#' x <- c("something", "something else", NA, "and another thing", NA)
#' rm_na(x)
rm_na <- function(x) {
  x[!is.na(x)]
}
