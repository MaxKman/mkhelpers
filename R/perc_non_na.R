#' Report percentage of non-na values in tbl
#'
#' @param tbl_x A tibble
#'
#' @return A summarised tibble reporting the percentage of non-na values for all columns in the input tibble
#' @export
#'
#' @examples
perc_non_na <- function(tbl_x) {
  tbl_x %>%
    summarise_all(~!is.na(.)) %>%
    summarise_all(sum) %>%
    mutate_all(~round(./nrow(tbl_x) * 100, 2)) %>%
    pivot_longer(everything(), names_to = "vars", values_to = "perc_non_na") %>%
    arrange(-perc_non_na)
}
