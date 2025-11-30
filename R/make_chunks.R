#' Dividing a vector into chunks
#'
#' @param x A vector (e.g. a character vector or numeric vector)
#' @param n_chunks Set this to an integer to divide x into n chunks. Setting n_chunks to zero ignores the parameter
#' @param size_chunks Set this to the desired size of chunks and n_chunks will be automatically determined
#'
#' @return A list of chunked vectors
#' @export
#'
#' @examples
#' library(mkhelpers)
#' char_vector <- purrr::map_chr(49:122, intToUtf8)
#' make_chunks(char_vector, n_chunks = 5)
#' make_chunks(char_vector, size_chunks = 3)
make_chunks <- function(x, n_chunks = 0, size_chunks = 1000) {
  if(n_chunks > 0) {
    chunks <- sort(rep_len(1:n_chunks, length(x))) %>% split(x, .)
    size_chunks <- length(chunks[[1]])
  } else {
    n_chunks <- ceiling(length(x) / size_chunks)
    chunks <- sort(rep_len(1:n_chunks, length(x))) %>% split(x, .)
  }
    return(chunks)
}
