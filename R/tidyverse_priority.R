#' Giving all tidyverse functions priority over conflicting functions from other packages
#'
#' @return
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(GenomicRanges)
#' tidyverse_priority()
tidyverse_priority <- function() {
  library(conflicted)
  conflicts <- conflict_scout()
  conflicts_tidyverse <- purrr::map(conflicts, ~any(. %in% tidyverse_packages())) %>% unlist
  conflicts <- conflicts[conflicts_tidyverse]
  iwalk(conflicts, ~conflict_prefer(.y, dplyr::intersect(.x, tidyverse_packages())[[1]]))
}
