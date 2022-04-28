#' Select fields in a string
#'
#' @param x A string to select fields from
#' @param sep The separator between fields
#' @param fields Which fields to choose, e.g. c(1,2,5) or 1:6
#' @param leading_sep Whether to add/preserve a leading separator
#' @param trailing_sep Whether to add/preserve a trailing separator

#'
#' @return A string
#' @export
#'
#' @examples
#' x <- "apple-banana-pear-mango-peach"
#' str_select_fields(x, "-", c(1,5))
str_select_fields <- function(x, sep, fields, leading_sep = FALSE, trailing_sep = FALSE) {
  tmp <- x %>% stringr::str_split(sep) %>% unlist %>% .[fields] %>% stringr::str_c(collapse = sep)
  if(leading_sep) {
    tmp <- tmp %>% stringr::str_c(sep, .)
  }
  if(trailing_sep) {
    tmp <- tmp %>% stringr::str_c(sep)
  }
  return(tmp)
}
