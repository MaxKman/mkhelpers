#' Add percentile column to tibble
#'
#' @description Adds a new column to the tibble, which shows the percentiles for a column of choice
#' @param tbl_x A tibble
#' @param var_x The tibble column for which to calculate percentiles
#'
#' @return A tibble
#' @export
#'
#' @examples
#' library(mkhelpers)
#' tbl_add_percentiles(mtcars, disp)
tbl_add_percentiles <- function(tbl_x, var_x) {
  tbl_x %>%
    mutate(id_tmp = row_number()) %>%
    arrange({{var_x}}) %>%
    mutate(rank_tmp = row_number(),
           "{{var_x}}_percentiles" := rank_tmp / length(rank_tmp) * 100) %>%
    arrange(id_tmp) %>%
    select(-id_tmp, -rank_tmp)
}
