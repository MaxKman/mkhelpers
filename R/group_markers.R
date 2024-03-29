#' Group marker gene identification
#' @description This function tries to identify marker genes for each group defined in the provided data. Importantly, it takes the identity of biological replicates (e.g. donors) into consideration.
#' @param m A count matrix, where rows are genes and columns are sample/group combinations. Single cell data should not be put in directly but aggregated first using aggregate_exp_data().
#' @param auto_extract_groups_samples If the count matrix was generated using aggregate_exp_data() groups and sample identities can be automatically extracted from the column names. Set to TRUE by default.
#' @param groups If auto_extract_groups_samples is set to FALSE, supply a character vector with group identities (e.g. clusters) for each column in count matrix.
#' @param samples If auto_extract_groups_samples is set to FALSE, supply a character vector with sample identity (e.g. donors) for each column in count matrix.
#' @param max_padj Upper cutoff for FDR-adjusted p value during marker selection.
#' @param min_log2f Lower cutoff for log2 fold change during marker selection.
#' @param log2f_dist This parameter is designed to prioritise 'exclusive' markers, which only mark one group in the data set. It determines how much stronger the log2 fold change for a given group must be in comparison to the group with second highest log2 fold change for a given gene. E.g. setting it to 2 will filter on group markers, whose log2 fold change is at least twice as strong in comparison to other groups.
#' @param max_n_markers Maximum number of marker genes returned for each group.
#' @param marker_direction Whether to return 'positive' or 'negative' markers for groups.
#' @param n_workers Number of cores to use when executing DEseq2
#' @param output_path All output will be stored in this folder
#' @param previous_output Path to a folder containing the output of a previous analysis. If this is provided only the different filtering criteria and cutoffs will be applied without re-running the differential gene expression analysis.
#'
#' @return The follogwing files are created in the chosen output path:
#' - deseq_dds.rds (contains the results of running DESeq2::DESeq())
#' - raw_results.csv (contains the formatted results for all input genes)
#' - selected_markers.csv (selected markers in csv format)
#' - selected_markers_list.rds (selected markers in a list format)
#' @export
#'
#' @examples
#' \dontrun{
#' # In this example we identify markers for different clusters in single cell PBMC data
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#'
#' # Simulate some donors
#' add_donors <- seu_pbmc@meta.data %>%
#'   rownames_to_column("cell_name") %>%
#'   group_by(seurat_clusters) %>%
#'   mutate(simulated_donors = sample(5, n(), replace = TRUE)) %>%
#'   ungroup %>%
#'   column_to_rownames("cell_name") %>%
#'   select(simulated_donors)
#' add_donors$cell_name <- rownames(add_donors)
#' seu_pbmc <- AddMetaData(seu_pbmc, add_donors)
#'
#' # Make pseudobulks
#' m_aggr <- aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = seu_pbmc@meta.data, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, n_cells_normalize = 10000, min_n_samples_aggr = 3, mode = "count")
#'
#' # Identify markers
#' group_markers(m = m_aggr, output_path = "~")
#'
#' # Plot results
#' genes_plot <- read_csv("~/selected_markers.csv") %>% group_by(group) %>% slice_max(log2fc_diff_max_2nd) %>% pull(gene)
#' tbl_x <- seu_extract_tbl(seu_x = seu_pbmc,
#'   reduction = "umap",
#'   metadata_cols = "seurat_clusters",
#'   extract_expr = TRUE,
#'   genes_extract = genes_plot,
#'   assay_extract = "RNA",
#'   slot_extract = "data",
#'   expr_format = "wide")
#' plot_umap_grid(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~/test_plot_cluster_markers_umap.png", feature_list = c("seurat_clusters", genes_plot), show_labels = TRUE, point_size = 0.7)
#' m_aggr_df <- aggregate_exp_data(m = seu_pbmc@assays$RNA@data, md = seu_pbmc@meta.data, aggr_col = seurat_clusters, sample_col = simulated_donors, n_cells_min = 20, min_n_samples_aggr = 3, mode = "mean", return_matrix = FALSE, expr_format = "long")
#' plot_feature_comparison_grid(tbl_x = m_aggr_df, feature_list = genes_plot, output_path = "~/test_plot_cluster_markers_boxplot.png", expr_col = "exp", group_col = "seurat_clusters", n_cols = 5)
#' }
group_markers <- function(m, auto_extract_groups_samples = TRUE, groups, samples, max_padj = 0.1, min_log2f = 1, log2f_dist = 2, max_n_markers = 0, marker_direction = c("positive", "negative"), n_workers = 4, output_path, previous_output = NULL) {

  if(is.null(previous_output)){
    if(auto_extract_groups_samples) {
      groups_samples <- colnames(m)
      groups <- str_remove(groups_samples, "_\\._.*")
      samples <- str_remove(groups_samples, ".*_\\._")
    } else {
      groups_samples <- str_c(groups, "_._", samples)
    }
    coldata <- data.frame(group = groups, sample = samples)
    rownames(coldata) <- groups_samples
    print("Running DESeq2")
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = m,
                                          colData = coldata,
                                          design = ~ group) ## see https://support.bioconductor.org/p/105087/
    dds <- DESeq2::DESeq(dds, betaPrior = TRUE, BPPARAM = MulticoreParam(workers = n_workers))

    print("Building results")

    result_groups <- groups %>% unique
    group_names <- result_groups
    result_groups <- result_groups %>% str_c("group", .)
    names(result_groups) <- result_groups

    buildresults_distinct <- function(deseq_obj, group) {
      results <- imap(group, function(group_x, group_x_name) {
        print(str_c("Group ", group_x_name))
        compare_against <- group[!(group %in% group_x)]
        DESeq2::results(deseq_obj, contrast = list(group_x, compare_against), listValues = c(1/length(group_x), -1/length(compare_against))) %>% as.data.frame() %>% mutate(gene = rownames(.), group = group_x_name)
      })
      bind_rows(results)
    }

    results <- buildresults_distinct(deseq_obj = dds, group = result_groups)
  } else {
    results <- read_csv(glue::glue("{previous_output}/raw_results.csv"))
  }
  print("Selecting markers.")

  top_markers <- results %>%
    group_by(gene) %>%
    mutate(base_mean_all = sum(baseMean)) %>%
    filter(base_mean_all > 0) %>%
    mutate(log2fc_max = max(log2FoldChange, na.rm = TRUE), log2fc_2nd = sort(log2FoldChange, decreasing = TRUE) %>% .[2]) %>%
    ungroup %>%
    mutate(log2fc_diff_max_2nd = log2fc_max - log2fc_2nd) %>%
    filter(log2FoldChange == log2fc_max) %>%
    filter(log2FoldChange >= min_log2f) %>%
    filter(padj <= max_padj)
  if(marker_direction[[1]] == "positive") {
    top_markers <- top_markers %>%
      filter(stat > 0)
  } else if(marker_direction[[1]] == "negative") {
    top_markers <- top_markers %>%
      filter(stat < 0)
  } else {
    stop("marker_direction needs to be either set to positive or negative")
  }

  top_markers <- top_markers %>% group_by(group) %>%
    filter(is.na(log2fc_diff_max_2nd) | log2fc_diff_max_2nd > log2(log2f_dist))
  if(max_n_markers > 0) {
    top_markers <- top_markers %>%
      slice_max(log2fc_diff_max_2nd, n = max_n_markers)
  }

  top_markers_list <- top_markers %>%
    split(.$group) %>%
    map("gene")
  if(is.null(previous_output)){
    saveRDS(dds, glue::glue("{output_path}/deseq_dds.rds"))
    write_csv(results, glue::glue("{output_path}/raw_results.csv"))
  }
  write_csv(top_markers, glue::glue("{output_path}/selected_markers.csv"))
  saveRDS(top_markers_list, glue::glue("{output_path}/selected_markers_list.rds"))

  iwalk(top_markers_list, function(markers_x, name_x) {
    print(glue::glue("Identified {length(markers_x)} marker genes for group {name_x}"))
  })
  print(glue::glue("Results can be found in {output_path}"))
}
