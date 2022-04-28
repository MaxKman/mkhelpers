str_select_fields <- function(x, sep, fields, trailing_sep = TRUE) {
  tmp <- x %>% stringr::str_split(sep) %>% unlist %>% .[fields] %>% stringr::str_c(collapse = sep)
  if(trailing_sep) {
    return(tmp %>% stringr::str_c(sep))
  } else {
    return(tmp)
  }
}
