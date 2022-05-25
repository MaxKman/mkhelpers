#' Assign the rownames of a data frame to the values in a column without removing the column
#'
#' @param df_x A data frame
#' @param col_name The column name that contains rownames
#'
#' @return A data frame with row names
#' @export
#'
#' @examples
#' car_df <- mtcars %>% rownames_to_column("car_name")
#' car_df <- car_df %>% copy_column_to_rownames("car_name")
copy_column_to_rownames <- function(df_x, col_name) {
  rownames(df_x) <- df_x %>% pull(col_name)
  return(df_x)
}
