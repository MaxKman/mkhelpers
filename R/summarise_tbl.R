#' Summarise the variables of a tibble
#'
#' @param tbl_x A tibble
#'
#' @return A tibble containing the data types of all columns of the input tibble, the number of non-na values, the number of unique values, the top ten values, three randomly drawn examples (cases 1-3, corresponding to entire rows of the original tibble) and one randomly drawn non-na value for each column.
#' @export
#'
#' @examples
#' summarise_tbl(mtcars)
summarise_tbl <- function(tbl_x) {
  n_unique_values <- tbl_x %>% map(~unique(.) %>% length) %>% unlist
  top_3_values <- colnames(tbl_x) %>%
    map(~count_(tbl_x,.) %>%
          slice_max(order_by = n, n = 10, with_ties = FALSE) %>%
          pull(1) %>%
          as.character %>%
          str_replace_na %>%
          str_c(collapse = "; ")) %>%
    unlist
  # The as.character conversion needs to be done for each column separately
  # otherwise errors are easily introduced, e.g. for dates
  draw_case <- function() {
    tbl_x %>% slice_sample %>% map(~.[[1]] %>% as.character) %>% unlist
  }
  case_1 <- draw_case()
  case_2 <- draw_case()
  case_3 <- draw_case()
  non_na_examples <- map_chr(tbl_x, function(var_x) {
    var_x_rm_na <- rm_na(var_x)
    if(length(var_x_rm_na) > 0) {
      return(sample(var_x_rm_na, 1) %>% as.character())
    }
    return(NA)
  })
  tibble(vars = colnames(tbl_x),
         types = tbl_x %>% map(class) %>% unlist,
         n_unique_values = n_unique_values,
         top_3_values = top_3_values,
         case_1 = case_1,
         case_2 = case_2,
         case_3 = case_3,
         non_na_examples = non_na_examples) %>%
    left_join(perc_non_na(tbl_x)) %>%
    relocate(perc_non_na, .after = types)
}