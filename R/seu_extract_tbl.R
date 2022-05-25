#' Extracting information stored in seurat objects into a tibble
#'
#' @param seu_x Seurat object
#' @param reduction Name of reduction to be extracted
#' @param metadata_cols Names of metadata columns to be extracted
#' @param extract_expr Whether to extract expression information (default: FALSE)
#' @param genes_extract Genes to extract, when extracting expression information (default: "all")
#' @param assay_extract Assay to extract expression information from (default: "RNA")
#' @param slot_extract Slot to extract expression information from (default: "data")
#' @param expr_format Return gene expression information in wide or long format (default: "wide")
#'
#' @return A tibble
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#'
#' seu_extract_tbl(seu_x = seu_pbmc,
#'   reduction = "umap",
#'   metadata_cols = c("nCount_RNA", "nFeature_RNA", "seurat_clusters"),
#'   extract_expr = TRUE,
#'   genes_extract = c("NCAM1", "SELL", "CD3D"),
#'   assay_extract = "RNA",
#'   slot_extract = "data",
#'   expr_format = "wide")
#' }
seu_extract_tbl <- function(seu_x, reduction = NULL, metadata_cols = "all", extract_expr = FALSE, genes_extract = "all", assay_extract = "RNA", slot_extract = "data", expr_format = c("wide", "long")) {

  if("cell_name" %in% colnames(seu_x@meta.data)) {
    seu_x@meta.data <- seu_x@meta.data %>% select(-cell_name)
  }

  # make a list for appending dataframes
  df_list <- list(tibble(cell_name = colnames(seu_x)))

  # extract reduction
  if(!is.null(reduction)) {
    if(reduction %in% names(seu_x@reductions)) {
      df_list <- seu_x@reductions[[reduction]]@cell.embeddings %>%
        as.data.frame() %>%
        rownames_to_column("cell_name") %>%
        as_tibble() %>%
        list %>%
        append(df_list, .)
    } else {
      stop("Reduction not found in seurat object!")
    }
  }

  # extract metadata
  if(metadata_cols[[1]] == "all") {metadata_cols <- colnames(seu_x@meta.data)}
  if(!is.null(metadata_cols)) {
    if(all(metadata_cols %in% colnames(seu_x@meta.data))) {
      df_list <- seu_x@meta.data %>%
        as.data.frame() %>%
        rownames_to_column("cell_name") %>%
        .[,c("cell_name", metadata_cols)] %>%
        as_tibble() %>%
        list %>%
        append(df_list, .)
    } else {
      stop(glue::glue("Error: The following columns are not found in metadata: {metadata_cols[!(metadata_cols %in% colnames(seu_x@meta.data))] %>% str_c(collapse = ', ')}\n"))
    }

  }

  # extract expression data
  if(extract_expr) {
    if(genes_extract[[1]] == "all") {
      genes_extract <- rownames(seu_x)
    }
    if(all(genes_extract %in% (seu_x@assays[[assay_extract]] %>% slot(slot_extract) %>% rownames))) {
      m <- seu_x@assays[[assay_extract]] %>% slot(slot_extract)
      m <- m[rownames(m) %in% genes_extract,]
      if(length(genes_extract) == 1) {
        tbl_temp <- tibble(cell_name = names(m), !!genes_extract := m)
      } else {
        tbl_temp <- m %>% as.matrix %>% t %>% as.data.frame() %>% rownames_to_column("cell_name") %>% as_tibble()
      }

      if(expr_format[[1]] == "wide") {
        df_list <- tbl_temp %>%
          list %>%
          append(df_list, .)
      } else if(expr_format[[1]] == "long") {
        df_list <- tbl_temp %>%
          pivot_longer(-cell_name, names_to = "gene", values_to = "expr") %>%
          list %>%
          append(df_list, .)
      }
    } else {
        stop(glue::glue("Error: The following genes are not found in the expression data: {genes_extract[!(genes_extract %in% (seu_x@assays[[assay_extract]] %>% slot(slot_extract) %>% rownames))] %>% str_c(collapse = ', ')}\n"))
    }
  }

  # turn list into one joined tibble
  suppressMessages(reduce(df_list, left_join))
}
