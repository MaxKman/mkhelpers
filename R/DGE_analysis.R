#' DEseq2 wrapper for pseudobulk differential gene expression analysis in single cell data
#' @param m A single cell matrix with genes as rows and cells as columns
#' @param md A metadata dataframe (must contain cluster_col, sample_col and group_col)
#' @param cluster_col The metadata column identifying the cluster of each cell
#' @param sample_col The metadata column identifying the sample to which each cell belongs (e.g. the donor)
#' @param group_col The metadata column identifying the group / condition to which each cell belongs (e.g. treatment vs. control)
#' @param batch_col The metadata column identifying the batch to which each cell belongs. Note, that to correct for a batch effect in the analysis the design formula needs to be adapted accordingly
#' @param balance_batches If set to TRUE, will ensure that only samples are included in the comparison for which a corresponding sample from the same batch is available in the other group. Requires batch_col to be defined.
#' @param add_var An additional column that should be preserved, for example because it is used in the design formula
#' @param title The title of comparison (will be used for file naming)
#' @param group1 The name of group 1 in group_col
#' @param group2 The name of group 2 in group_col
#' @param design The DEseq2 design formula to be used
#' @param n_cells_normalize The number of cells to which the pseudobulking will be standardized. E.g. if set to 10000, the counts aggregated from a cluster/sample combination will be normalized as follows: count_matrix / n_cells * n_cells_normalize.
#' @param n_cells_min The minimum number of cells for a cluster/sample combination to be included in the analysis. E.g. if set to 20, samples will be excluded for clusters, for which they contain less than 20 cells
#' @param min_n_samples_group Minimum number of samples in each group to perform DGE analysis. E.g. if set to 3 only 3 versus 3 comparisons will be carried out.
#' @param cell_name_col The name of the column that contains the cell identifier
#' @param exp_percentage Minimum percentage of cells in a cluster from one group, which need to express a gene for it to be included in the DGE analysis
#' @param exp_percentage_strict If set to TRUE, only those genes will be included in the DGE analysis, which are expressed in exp_percentage cells in each sample in a cluster from one group. If set to FALSE the inclusion criterion is an average expression in exp_percentage cells across samples (default: FALSE)
#' @param exp_percentage_type Whether to select genes, which meet the exp_percentage criterion in both groups ('intersect') or in one or both groups ('union'). Default: 'intersect'.
#' @param aggr_counts_non_zero_percentage This parameter is used to exclude genes from the DGE analysis, for which less than a given percentage of samples have zero expression after count aggregation. This is helpful in stimulation experiments, where some genes are hardly detectable in the unstimulated condition, which can lead to spurious DGE results due to inflated fold changes and low variation. Default: 90%.
#' @param lfc_threshold Log2 fold change treshold passed to DESeq2. See ?DESeq2::results.
#' @param save_results Whether to save results to an RDS file and a log file containing the settings used in this run (default: FALSE)
#' @param save_raw_deseq_objs Save raw DEseq 2 objects, to allow for costum contrast extraction
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
#' # Simulate some batches and donors
#' add_donors <- pbmc_ifnb@meta.data %>%
#'   rownames_to_column("cell_name") %>%
#'   group_by(stim) %>%
#'   mutate(simulated_batches = sample(4, n(), replace = TRUE)) %>%
#'   ungroup %>%
#'   mutate(simulated_donors = simulated_batches) %>%
#'   mutate(simulated_donors = str_c(simulated_donors, "_", stim)) %>%
#'   column_to_rownames("cell_name") %>%
#'   select(simulated_donors, simulated_batches)
#' add_donors$cell_name <- rownames(add_donors)
#' pbmc_ifnb <- AddMetaData(pbmc_ifnb, add_donors)
#'
#' # Remove a donor from one batch to test whether the entire batch is
#' # correctly removed from the analysis (if balance_batches is set to true)
#' donor_remove <- pbmc_ifnb@meta.data$simulated_donors %>% unique %>% .[[1]]
#' pbmc_ifnb <- subset(pbmc_ifnb, simulated_donors != donor_remove)
#'
#' test_dge <- DGE_analysis(m = pbmc_ifnb@assays$RNA@counts,
#'                          md = pbmc_ifnb@meta.data,
#'                          title = "testrun",
#'                          cluster_col = seurat_annotations,
#'                          sample_col = simulated_donors,
#'                          group_col = stim,
#'                          group1 = "CTRL",
#'                          group2 = "STIM",
#'                          design = ~simulated_batches + stim,
#'                          balance_batches = TRUE,
#'                          batch_col = simulated_batches,
#'                          n_cells_normalize = 10000,
#'                          n_cells_min = 10,
#'                          min_n_samples_group = 3,
#'                          cell_name_col = cell_name,
#'                          exp_percentage = 1)
#'  }
DGE_analysis <- function(m, md, m_norm, cluster_col, sample_col, group_col, batch_col = NULL, balance_batches = FALSE, add_var = NULL, title, group1, group2, design = ~group, n_cells_normalize, n_cells_min, min_n_samples_group, cell_name_col = cell_name, exp_percentage = 5, exp_percentage_strict = FALSE, exp_percentage_type = c('intersect', 'union'), aggr_counts_non_zero_percentage = 90, lfc_threshold = 0, save_results = FALSE, save_raw_deseq_objs = FALSE, savepath = "../../RDS/DGE") {
  sample_col_str <- deparse(substitute(sample_col))
  md <- md %>% mutate({{group_col}} := factor({{group_col}}, levels = c(group1, group2)))

  md_persample <- md %>% select({{sample_col}}, {{group_col}}, {{add_var}}, {{batch_col}}) %>% distinct() %>% ungroup %>% as.data.frame()
  rownames(md_persample) <- md_persample %>% pull({{sample_col}})

  samples_g1 <- md %>% filter({{group_col}} == group1) %>% pull({{sample_col}}) %>% unique
  samples_g2 <- md %>% filter({{group_col}} == group2) %>% pull({{sample_col}}) %>% unique

  md <- md %>% group_by({{cluster_col}})
  cells_cluster_sample <- md %>% group_split() %>% map(function(x) x %>% split(.[,sample_col_str])) %>% map_depth(~pull(., {{cell_name_col}}), .depth = 2)
  names(cells_cluster_sample) <- group_keys(md) %>% pull %>% as.character()
  md <- md %>% ungroup

  # Remove samples separately for each cluster with less than n_cells_min
  cells_cluster_sample <- map(cells_cluster_sample, function(x) {
    x[x %>% map(length) > n_cells_min]
  })

  # ensure that only samples are included in the comparison for which a corresponding sample from the same batch is available in the other group. Requires batch_col to be defined.
  if(balance_batches) {
    for (k in seq_along(cells_cluster_sample)) {
      batch_info <- md %>%
        filter({{cluster_col}} == names(cells_cluster_sample)[[k]]) %>%
        select({{sample_col}}, {{group_col}}, {{batch_col}}) %>%
        distinct %>%
        filter({{sample_col}} %in% names(cells_cluster_sample[[k]]))
      suppressMessages(
        batches_both_groups <- batch_info %>%
          group_by({{group_col}}, {{batch_col}}) %>%
          summarise(n_batch = n()) %>%
          pivot_wider(names_from = {{group_col}}, values_from = n_batch) %>%
          drop_na %>%
          left_join(batch_info)
      )
      donors_keep <- batches_both_groups %>% pull({{sample_col}})
      if(nrow(batches_both_groups) >= 1) {
        cells_cluster_sample[[k]] <- cells_cluster_sample[[k]][donors_keep]
      } else {
        cells_cluster_sample[[k]] <- NA
      }
    }
  }

  # Remove clusters, for which less than min_n_samples_group samples remain in a group
  for (k in seq_along(cells_cluster_sample)) {
    n_g1 <- names(cells_cluster_sample[[k]]) %in% samples_g1 %>% sum
    n_g2<- names(cells_cluster_sample[[k]]) %in% samples_g2 %>% sum
    if(n_g1 < min_n_samples_group || n_g2 < min_n_samples_group) {cells_cluster_sample[[k]] <- NA}
  }
  # Print message to note which samples are kept
  for (k in seq_along(cells_cluster_sample)) {
    if(!is.na(cells_cluster_sample[[k]][[1]][[1]])) {
      print(glue::glue("Cluster {names(cells_cluster_sample)[[k]]}:\n{names(cells_cluster_sample[[k]]) %>% str_c(collapse = ', ')}\nare kept!\n\n"))
    } else {
      print(glue::glue("Cluster {names(cells_cluster_sample)[[k]]}: too few samples, removed from analysis!\n\n"))
    }
  }
  cells_cluster_sample <- cells_cluster_sample[!is.na(cells_cluster_sample)]

  print("Extracting expression values from count matrix. This can take a while...")
  cells_cluster_sample_expr <- map(cells_cluster_sample, .progress = TRUE, function(x) {
    map(x, ~m[,.])
  })

  tictoc::tic()
  results <- suppressMessages(imap(cells_cluster_sample_expr, function(x, name_x) {
    x <- map(x, as.matrix)

    if(exp_percentage_strict) {
      # Subset on genes, which are at least expressed in x% of cells of each sample within a group
      sample_genes <- map(x, ~rownames(.[Matrix::rowSums(. > 0) >= ncol(.)*exp_percentage/100,]))
      genes_1 <- sample_genes[names(sample_genes) %in% samples_g1] %>% purrr::reduce(.x = ., .f = intersect)
      genes_2 <- sample_genes[names(sample_genes) %in% samples_g2] %>% purrr::reduce(., intersect)
      if(exp_percentage_type[[1]] == 'union') {
        genes_select <- union(genes_1, genes_2)
      } else {
        genes_select <- intersect(genes_1, genes_2)
      }

    } else {
      # Subset on genes, which are on average expressed in exp_percentage cells in each group
      rowsumsdf <- map_dfc(x, ~Matrix::rowSums(. > 0) / ncol(.))
      genes_select <- map(list(samples_g1, samples_g2), function(samples_x) {
        genes_select <- rowsumsdf[,colnames(rowsumsdf) %in% samples_x] %>% rowMeans()
        names(genes_select) <- rownames(x[[1]])
        genes_select <- genes_select[genes_select >= exp_percentage/100]
        return(names(genes_select))
      })
      if(exp_percentage_type[[1]] == 'union') {
        genes_select <- genes_select %>% purrr::reduce(union)
      } else {
        genes_select <- genes_select %>% purrr::reduce(intersect)
      }
    }

    if(length(genes_select) == 0) {
      return(data.frame())
    }

    x <- map(x, ~.[genes_select,])

    cnames <- names(x)
    n_cells <- map(x, ncol)
    #aggregate counts
    x <- map(x, ~Matrix::rowSums(.))
    #normalize
    x <- map2(x, n_cells, function(x_x, n_cells_x) {
      x_x <- x_x / n_cells_x * n_cells_normalize
      round(x_x)
    }) %>% purrr::reduce(cbind)
    colnames(x) <- cnames

    # Select genes, whose aggregated counts are non-zero in > x% of samples
    genes_select_greater_zero <- map(list(samples_g1, samples_g2), function(samples_x) {
      greater_zero <- (x[,colnames(x) %in% samples_x] > 0) %>% as.matrix()
      greater_zero <- rowSums(greater_zero) / ncol(greater_zero)
      greater_zero[greater_zero > aggr_counts_non_zero_percentage / 100] %>% names
    }) %>% purrr::reduce(intersect)
    x <- x[genes_select_greater_zero,]

    # Calculating log2fc and applying cutoff
    if(lfc_threshold > 0) {
      aggr_exp <- map_dfr(colnames(x), function(sample_x) {
        cells_aggr <- md %>% filter({{sample_col}} == sample_x, {{cluster_col}} == name_x) %>% pull({{cell_name_col}})
        m_norm_sub <- m_norm[rownames(m_norm) %in% rownames(x), colnames(m_norm) %in% cells_aggr]
        tibble(gene = rownames(m_norm_sub), exp =  Matrix::rowMeans (m_norm_sub), sample = sample_x) %>%
          left_join(md_persample %>% select({{sample_col}}, {{group_col}}))
      }) %>%
        summarise(mean_exp = mean(exp), .by = c(gene, {{group_col}})) %>%
        pivot_wider(names_from = {{group_col}}, values_from = mean_exp) #%>%
        #mutate(log2fc = log2(group1 / group2))
      View(aggr_exp)
    }
    return(NULL)

    print(str_c("Cluster ", name_x, ": Performing DGE analysis on ", nrow(x), " genes!"))

    fData <- data.frame(names = genes_select)
    rownames(fData) <- genes_select
    cData <- md_persample %>% filter(get(sample_col_str) %in% colnames(x))
    x <- x[,rownames(cData)]
    DEG <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                  colData = cData,
                                  design = design)

    DEG <- DESeq2::DESeq(DEG)
    if(save_raw_deseq_objs) {
      dir.create(str_c(savepath, "/deseq2_objects"), showWarnings = FALSE)
      saveRDS(DEG, str_c(savepath, "/deseq2_objects/", "/DGE_", n_cells_normalize, "_cells_", title, "_", name_x, ".RDS"))
    }
    res <- DESeq2::results(DEG, pAdjustMethod = "fdr", lfcThreshold = lfc_threshold) %>% as.data.frame() %>% mutate(gene = rownames(.))
    return(res)
  }))

  results <- map2(.x = names(results), .y = results, ~mutate(.y, {{ cluster_col }} := .x))
  tictoc::toc()

  # the following is way faster than using purrr:reduce calls
  results <- results %>% do.call(rbind, .) %>% as_tibble
  if(save_results) {
    dir.create(savepath, recursive = TRUE)
    saveRDS(results, str_c(savepath, "/DGE_", n_cells_normalize, "_cells_", title, ".RDS"))
    cluster_col_str <- deparse(substitute(cluster_col))
    group_col_str <- deparse(substitute(group_col))
    write_lines(file = str_c(savepath, "/DGE_", n_cells_normalize, "_cells_", title, ".log"),
                glue::glue("Title: {title}",
                     "Cluster column: {cluster_col_str}",
                     "Sample column: {sample_col_str}",
                     "Group column: {group_col_str}",
                     "Group 1: {group1}",
                     "Group 2: {group2}",
                     "Design: {deparse1(design)}",
                     "Balance batches: {balance_batches}",
                     "Number of cells to normalize to: {n_cells_normalize}",
                     "Minimum cells per sample: {n_cells_min}",
                     "Minimum number of samples per group: {min_n_samples_group}",
                     "Minimum expression percentage: {exp_percentage}",
                     .sep = "\n"))
  }
  return(results)
}
