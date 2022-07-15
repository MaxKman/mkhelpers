#' Summarise the column data types of a tibble
#'
#' @param tbl_x A tibble
#'
#' @return A tibble containing the data types of all columns of the input tibble as well as three randomly drawn examples (cases 1-3, corresponding to entire rows of the original tibble) and one randomly drawn non-na value for each column.
#' @export
#'
#' @examples
#' summarise_data_types(mtcars)
summarise_data_types <- function(tbl_x) {
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
  tibble(var_names = colnames(tbl_x),
         var_types = tbl_x %>% map(class) %>% unlist,
         case_1 = case_1,
         case_2 = case_2,
         case_3 = case_3,
         non_na_examples = non_na_examples)
}
