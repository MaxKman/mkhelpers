#' Add additional ggplot layers to a list of existing ggplot objects
#'
#' @param gg_list A list of ggplot objects
#' @param arguments_add ggplot layers to add provided as a comma-separated list. Do not use + signs between layers.
#' @param subset Define a subset of the list to add the layers to, eg. 1:3, or c(1,5,7)
#'
#' @return
#' @export
#'
#' @examples
add_to_gg_list <- function(gg_list, arguments_add, subset = NULL) {
  if(is.null(subset)) {
    subset <- 1:length(gg_list)
  }
  for(k in subset) {
    gg_list[[k]] <- gg_list[[k]] + arguments_add
  }
  return(gg_list)
}
