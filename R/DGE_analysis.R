#' DEseq2 wrapper for pseudobulk differential gene expression analysis in single cell data
#' @description This function can be run in parallel by setting up a future environment. See examples.
#' @param m A single cell matrix with genes as rows and cells as columns
#' @param md A metadata dataframe (must contain cluster_col, sample_col and group_col)
#' @param cluster_col The metadata column identifying the cluster of each cell
#' @param sample_col The metadata column identifying the sample to which each cell belongs (e.g. the donor)
#' @param group_col The metadata column identifying the group / condition to which each cell belongs (e.g. treatment vs. control)
#' @param add_var An additional column that should be preserved, for example because it is used in the design formula (e.g. batch)
#' @param title The title of comparison (will be used for file naming)
#' @param group1 The name of group 1 in group_col
#' @param group2 The name of group 2 in group_col
#' @param design The DEseq2 design formula to be used
#' @param n_cells_normalize The number of cells to which the pseudobulking will be standardized. E.g. if set to 10000, the counts aggregated from a cluster/sample combination will be normalized as follows: count_matrix / n_cells * n_cells_normalize.
#' @param n_cells_min The minimum number of cells for a cluster/sample combination to be included in the analysis. E.g. if set to 20, samples will be excluded for clusters, for which they contain less than 20 cells
#' @param min_n_samples_group Minimum number of samples in each group to perform DGE analysis. E.g. if set to 3 only 3 versus 3 comparisons will be carried out.
#' @param cell_name_col The name of the column that contains the cell identifier
#' @param exp_percentage Minimum percentage of cells in a cluster from one group, which need to express a gene for it to be included in the DGE analysis
#' @param save_results Whether to save results to an RDS file (default: FALSE)
#' @param savepath Path where the results of the analysis should be saved
#'
#' @return A dataframe with the DEseq2 results with an additional column indicating the cluster
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratData)
#' InstallData("ifnb")
#' pbmc_ifnb <- LoadData("ifnb")
#'
#' # Simulate some donors
#' add_donors <- pbmc_ifnb@meta.data %>%
#'   rownames_to_column("cell_name") %>%
#'   group_by(stim) %>%
#'   mutate(simulated_donors = sample(3, n(), replace = TRUE)) %>%
#'   ungroup %>%
#'   mutate(simulated_donors = str_c(simulated_donors, "_", stim)) %>%
#'   column_to_rownames("cell_name") %>%
#'   select(simulated_donors)
#'add_donors$cell_name <- rownames(add_donors)
#'pbmc_ifnb <- AddMetaData(pbmc_ifnb, add_donors)
#'# run sequential:
#'options(future.globals.maxSize = 20*1024^3)
#'future::plan(future::sequential())
#'tictoc::tic()
#'test_dge <- DGE_analysis(m = pbmc_ifnb@assays$RNA@counts,
#'  md = pbmc_ifnb@meta.data,
#'  cluster_col = seurat_annotations,
#'  sample_col = simulated_donors,
#'  group_col = stim,
#'  group1 = "CTRL",
#'  group2 = "STIM",
#'  design = ~stim,
#'  n_cells_normalize = 10000,
#'  n_cells_min = 10,
#'  min_n_samples_group = 3,
#'  cell_name_col = cell_name,
#'  exp_percentage = 1)
#'tictoc::toc()
#'# run in parallel:
#'options(future.globals.maxSize = 20*1024^3)
#'future::plan(future::multisession(workers = 10, gc = TRUE))
#'tictoc::tic()
#'test_dge <- DGE_analysis(m = pbmc_ifnb@assays$RNA@counts,
#'  md = pbmc_ifnb@meta.data,
#'  cluster_col = seurat_annotations,
#'  sample_col = simulated_donors,
#'  group_col = stim,
#'  group1 = "CTRL",
#'  group2 = "STIM",
#'  design = ~stim,
#'  n_cells_normalize = 10000,
#'  n_cells_min = 10,
#'  min_n_samples_group = 3,
#'  cell_name_col = cell_name,
#'  exp_percentage = 1)
#'tictoc::toc()
#'future:::ClusterRegistry("stop")
#'  }
DGE_analysis <- function(m, md, cluster_col, sample_col, group_col, add_var = NULL, title, group1, group2, design = ~batch_pair + group, n_cells_normalize, n_cells_min, min_n_samples_group, cell_name_col = cell_name, exp_percentage = 5, save_results = FALSE, savepath = "../../RDS/DGE") {
  sample_col_str <- deparse(substitute(sample_col))
  md <- md %>% mutate({{group_col}} := factor({{group_col}}, levels = c(group1, group2)))

  md_persample <- md %>% select({{sample_col}}, {{group_col}}, {{ add_var }}) %>% distinct() %>% ungroup %>% as.data.frame()
  rownames(md_persample) <- md_persample %>% pull({{sample_col}})

  samples_g1 <- md %>% filter({{group_col}} == group1) %>% pull({{sample_col}}) %>% unique
  samples_g2 <- md %>% filter({{group_col}} == group2) %>% pull({{sample_col}}) %>% unique

  md <- md %>% group_by({{cluster_col}})
  cells_cluster_sample <- md %>% group_split() %>% map(function(x) x %>% split(.[,sample_col_str])) %>% map_depth(~pull(., {{cell_name_col}}), .depth = 2)
  names(cells_cluster_sample) <- group_keys(md) %>% pull %>% as.character()

  #Remove samples seperately for each cluster with less than n_cells_min
  cells_cluster_sample <- map(cells_cluster_sample, function(x) {
    x[x %>% map(length) > n_cells_min]
  })

  # Remove clusters, for which less than min_n_samples_group samples remain in a group
  for (k in seq_along(cells_cluster_sample)) {
    n_g1 <- names(cells_cluster_sample[[k]]) %in% samples_g1 %>% sum
    n_g2<- names(cells_cluster_sample[[k]]) %in% samples_g2 %>% sum
    if(n_g1 < min_n_samples_group || n_g2 < min_n_samples_group) {cells_cluster_sample[[k]] <- NA}
  }
  cells_cluster_sample <- cells_cluster_sample[!is.na(cells_cluster_sample)]

  cells_cluster_sample_expr <- map(cells_cluster_sample, function(x) {
    map(x, ~m[,.])
  })

  tictoc::tic()
  results <- suppressMessages(furrr::future_map(cells_cluster_sample_expr, function(x) {
    x <- map(x, as.matrix) # somehow the parallel execution needs this to work
    # Subset on genes, which are at least expressed in x% of cells of each sample within a group
    sample_genes <- map(x, ~rownames(.[Matrix::rowSums(. > 0) >= ncol(.)*exp_percentage/100,]))

    rows_1 <- sample_genes[names(sample_genes) %in% samples_g1] %>% purrr::reduce(.x = ., .f = intersect)
    rows_2 <- sample_genes[names(sample_genes) %in% samples_g2] %>% purrr::reduce(., intersect)
    rows <- c(rows_1, rows_2) %>% unique()
    print(str_c("Performing DGE analysis on ", length(rows), " genes!"))
    if(length(rows) == 0) {
      return(data.frame())
    }
    x <- map(x, ~.[rows,])
    cnames <- names(x)
    n_cells <- map(x, ncol)
    #aggregate counts
    x <- map(x, ~Matrix::rowSums(.))
    #normalise
    x <- map2(x, n_cells, function(x_x, n_cells_x) {
      x_x <- x_x / n_cells_x * n_cells_normalize
      round(x_x)
    }) %>% purrr::reduce(cbind)
    colnames(x) <- cnames

    fData <- data.frame(names = rows)
    rownames(fData) <- rows
    cData <- md_persample %>% filter(get(sample_col_str) %in% colnames(x))
    x <- x[,rownames(cData)]
    DEG <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                  colData = cData,
                                  design = design)

    DEG <- DESeq2::DESeq(DEG)
    res <- DESeq2::results(DEG, pAdjustMethod = "fdr") %>% as.data.frame() %>% mutate(gene = rownames(.))
    return(res)
  }))

  results <- map2(.x = names(results), .y = results, ~mutate(.y, {{ cluster_col }} := .x))
  tictoc::toc()

  # the following is way faster than using purrr:reduce calls
  results <- results %>% do.call(rbind, .) %>% as_tibble
  if(save_results) {
    dir.create(savepath, recursive = TRUE)
    saveRDS(results, str_c(savepath, "/DGE_", n_cells_normalize, "_cells_", deparse(substitute(cluster_col)), "_", title, ".RDS"))
  }
  return(results)
}
