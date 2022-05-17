#' Split a dataframe into a list, where each entry is one row of the dataframe
#'
#' @param df The input dataframe
#'
#' @return A list of rows
#' @export
#'
#' @examples
row_list <- function(df) {
  split(df, 1:nrow(df))
}
