#' Test for gene set enrichment using AUCell
#' @description This is a wrapper for AUCell (see https://bioconductor.org/packages/release/bioc/html/AUCell.html). It first ranks all genes for each individual sample / cell by their expression. In a second step gene set enrichment is determined as described here https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html. When rankings have been calculated once for a given data set they don't have to be re-run each time a different gene set with different parameters is tested. Therefore it makes sense to store these rankings and reload them when needed
#' @param m A gene expression matrix with genes as rows and samples / cells as columns
#' @param geneset_list A named list of gene sets to test
#' @param chunk_data Whether to process the data in multiple chunks, which can be useful to for very large matrices and allows for parallel processing (set future options accordingly in your environment, e.g. plan(future::multisession(workers = 4)))
#' @param chunk_size The number of samples / cells assigned to each chunk. Default: 1000
#' @param rankings_save_path Where to save the results of the expression ranking.
#' @param saved_rankings If a file path is provided here saved rankings are loaded and expression ranking is skipped.
#' @param rankings_object The rankings can also be directly loaded from an object in the current environment provided here
#' @param aucMaxRank_perc Cutoff percentage of top genes in the ranking considered for gene set enrichment analysis
#' @param expr_format Return gene expression information in wide or long format (default: "wide")
#' @param min_n_genes_detected A minimum number of genes found in the expression data for inclusion of a geneset. This only considers whether a gene is reported in the expression matrix and filtering such as removing genes with all zeroes needs to be performed beforehand.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' seu_pbmc <- LoadData("pbmc3k", "pbmc3k.final")
#'
#' celltype_markers <- list(
#'   b_cells = c("CD79B", "SNHG7", "MS4A1", "BANK1", "CD79A", "RCSD1", "CYB561A3", "HLA-DQB1", "HVCN1", "RPL22L1"),
#'   nk_cells = c("GNLY", "PRF1", "APMAP", "FGFBP2", "GZMB", "C5orf56", "AKR1C3", "SPON2", "TTC38", "IL2RB"),
#'   monocytes = c("LYZ", "S100A9", "GSTP1", "S100A8", "GPX1", "GRN", "CEBPD", "CAPG", "NFKBIA", "MIR24-2")
#'   )
#' results <- geneset_enrichment(m = seu_pbmc@assays$RNA@data, geneset_list = celltype_markers, chunk_data = FALSE, rankings_save_path = "~/saved_rankings.rds")
#'
#' # Plot results
#' tbl_x <- seu_extract_tbl(seu_x = seu_pbmc,
#'                          reduction = "umap",
#'                          metadata_cols = "seurat_clusters",
#'                          extract_expr = FALSE,
#'                          expr_format = "wide")
#' tbl_x <- tbl_x %>% left_join(results %>% rename(cell_name = sample))
#' plot_umap_grid(tbl_x = tbl_x, umap_dim_col_1 = UMAP_1, umap_dim_col_2 = UMAP_2, output_path = "~/test_plot_geneset_enrichment.png", feature_list = c("seurat_clusters", names(celltype_markers)), show_labels = TRUE, point_size = 0.7, n_cols = 2)
#'
#' tbl_x_long <- tbl_x %>% pivot_longer(all_of(names(celltype_markers)), names_to = "geneset", values_to = "AUC")
#' plot_feature_comparison_grid(tbl_x = tbl_x_long, feature_list = names(celltype_markers), feature_col = "geneset", output_path = "~/test_plot_geneset_enrichment_boxplot.png", expr_col = "AUC", group_col = "seurat_clusters", n_cols = 3, plot_width = 40, plot_height = 30, point_size = 0.5)
#' }
geneset_enrichment <- function(m, geneset_list, chunk_data = FALSE, chunk_size = 10000, rankings_save_path, saved_rankings = "none", rankings_object = NULL, aucMaxRank_perc = 20, expr_format = c("wide", "long"), min_n_genes_detected = 0) {
  set.seed(123)

  # Apply cutoff for number of detected genes
  if(min_n_genes_detected > 0) {
    genesets_keep <- imap_chr(geneset_list, function(genes_x, name_x) {
      n_genes_detected <- (rownames(m) %in% genes_x) %>% sum()
      if(n_genes_detected >= min_n_genes_detected) {
        return(name_x)
      }
      return(NA)
    }) %>% rm_na()

    geneset_list <- geneset_list[genesets_keep]

    if(length(geneset_list) == 0) {
      print("No genesets left for testing after applying min_n_genes_detected")
      return(NA)
    }
  }

  if(!chunk_data) {
    if(saved_rankings == "none" & is.null(rankings_object)) {
      print("...building rankings...")
      AU_rankings <- AUCell::AUCell_buildRankings(m, verbose = TRUE, nCores = 2)
      print("done!")
      saveRDS(AU_rankings, rankings_save_path)
    } else {
      if(is.null(rankings_object)) {
        AU_rankings <- readRDS(saved_rankings)
      } else {
        AU_rankings <- rankings_object
      }
    }
    print("...calculating AUC...")
    AUC <- AUCell::AUCell_calcAUC(geneSets = geneset_list, rankings = AU_rankings, aucMaxRank=nrow(m)*aucMaxRank_perc/100)
    print("done!")
    AUC_results <- AUC@assays@data$AUC %>% t %>% as.data.frame %>% rownames_to_column("sample") %>% as_tibble()
  } else {
    if(saved_rankings == "none" & is.null(rankings_object)) {
      #Make chunks
      print("...chunking the data...")
      chunk_list <- make_chunks(colnames(m), size_chunks = chunk_size)
      m_chunks <- map(chunk_list, ~m[,.])
      print("done!")
      print("...building rankings...")
      AU_rankings_list <- furrr::future_map(m_chunks, function(x) {
        print("...another chunk ranked...")
        AUCell::AUCell_buildRankings(x, verbose = TRUE, nCores = 2) #nCores = 1 uses all cores for some reason!
      })
      print("done!")

      #I would have liked to bind the chunks again, but cbind doesn't really work as described here https://github.com/aertslab/AUCell/issues/5 Therefore AUC calculations need to be mapped as well
      saveRDS(AU_rankings_list, rankings_save_path)
    } else {
      if(is.null(rankings_object)) {
        AU_rankings_list <- readRDS(saved_rankings)
      } else {
        AU_rankings_list <- rankings_object
      }
    }
    print("...calculating AUC...")
    AUC_results <- map(AU_rankings_list, function(ranking_x) {
      AUC <- AUCell::AUCell_calcAUC(geneSets = geneset_list, rankings = ranking_x, aucMaxRank=nrow(m)*aucMaxRank_perc/100)
      AUC@assays@data$AUC %>% t %>% as.data.frame %>% rownames_to_column("sample") %>% as_tibble()
    })
    AUC_results <- bind_rows(AUC_results)
    print("done!")
  }
  if(expr_format[[1]] == "long") {
    AUC_results <- AUC_results %>% pivot_longer(-sample, names_to = "geneset", values_to = "AUC")
  }
  return(AUC_results)
}
