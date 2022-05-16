#' Title Filter on genes expressed in more than n samples (e.g. cells)
#'
#' @param m A count matrix (can be a sparse matrix)
#' @param min_perc_exp The minimum percentage of samples expressing a given gene for it to be included in the output
#'
#' @return A filtered count matrix
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#' filter_genes(m = pbmc3k@assays$RNA@counts, min_perc_exp = 5)
#' }
filter_genes <- function(m, min_perc_exp) {
  m[Matrix::rowSums(m > 0) > ncol(m)*min_perc_exp/100,]
}
