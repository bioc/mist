#' Differential methylation evaluation over time for scDNA-seq data under the 2-group scenario.
#'
#' This function performs differential methylation (DM) analysis to identify genomic features
#' showing significant changes between two groups along pseudotime. The function models
#' methylation changes and compares the fitted curves for each group, calculating the integral
#' of the differences between the curves.
#'
#' @param Dat_sce_g1 A `SingleCellExperiment` object for group 1, containing:
#' - `mist_pars`: A numeric matrix in `rowData` with estimated parameters for all genomic features (generated by \code{estiParamTwo}).
#' @param Dat_sce_g2 A `SingleCellExperiment` object for group 2, containing:
#' - `mist_pars`: A numeric matrix in `rowData` with estimated parameters for all genomic features (generated by \code{estiParamTwo}).
#' 
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend for computations, as used in `bplapply()`. Defaults to `MulticoreParam()` for parallel processing.
#' 
#' @return A named numeric vector where each value corresponds to a genomic feature (e.g., a gene).
#' The values represent the integral of the differences between the fitted curves of scDNA methylation
#' levels for the two groups. The vector is sorted in descending order, with larger values indicating
#' more drastic differences between the groups.
#'
#' @import MCMCpack BiocParallel car mvtnorm
#' @importFrom S4Vectors subjectHits queryHits Rle
#' @importFrom methods is
#' @importFrom rtracklayer start offset end
#' @importFrom Matrix cov2cor toeplitz update
#' @importFrom stats pgamma poly qgamma rnorm runif
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data_g1 <- readRDS(system.file("extdata", "group1_sampleData_sce.rds", package = "mist"))
#' data_g2 <- readRDS(system.file("extdata", "group2_sampleData_sce.rds", package = "mist"))
#' Dat_sce_g1 <- estiParam(
#'     Dat_sce = data_g1,
#'     Dat_name = "Methy_level_group1",
#'     ptime_name = "pseudotime"
#' )
#'
#' Dat_sce_g2 <- estiParam(
#'     Dat_sce = data_g2,
#'     Dat_name = "Methy_level_group2",
#'     ptime_name = "pseudotime"
#' ) 
#' # Run differential methylation analysis
#' dm_results <- dmTwoGroups(
#'     Dat_sce_g1 = Dat_sce_g1,
#'     Dat_sce_g2 = Dat_sce_g2
#' )
dmTwoGroups <- function(Dat_sce_g1,
                        Dat_sce_g2,
                        BPPARAM = MulticoreParam()) {
  ######## 1. Input Validation
  # Check if inputs are SingleCellExperiment objects
  if (!methods::is(Dat_sce_g1, "SingleCellExperiment") || !methods::is(Dat_sce_g2, "SingleCellExperiment")) {
    stop("Both Dat_sce_g1 and Dat_sce_g2 must be SingleCellExperiment objects.",
         call. = TRUE, domain = NULL)
  }
  
  ######## 2. Extract Parameter Matrices
  mist_pars_matrix1 <- rowData(Dat_sce_g1)$mist_pars
  mist_pars_matrix2 <- rowData(Dat_sce_g2)$mist_pars
  
  if (is.null(mist_pars_matrix1) || is.null(mist_pars_matrix2)) {
    stop("Missing mist_pars matrices in one or both SCE objects.")
  }
  
  ######## 3. Identify Common Features
  common_features <- intersect(rownames(mist_pars_matrix1), rownames(mist_pars_matrix2))
  
  if (length(common_features) == 0) {
    stop("No common features found between the two groups.")
  }
  
  mist_pars_matrix1 <- mist_pars_matrix1[common_features, , drop = FALSE]
  mist_pars_matrix2 <- mist_pars_matrix2[common_features, , drop = FALSE]
  
  ######## 4. Split into Lists for Parallel Processing
  beta_sigma_list_group <- list(
    Group1 = split(as.data.frame(mist_pars_matrix1), rownames(mist_pars_matrix1)),
    Group2 = split(as.data.frame(mist_pars_matrix2), rownames(mist_pars_matrix2))
  )
  
  ######## 5. Validate Features
  valid_features <- vapply(common_features, function(feature) {
    all(is.finite(as.numeric(beta_sigma_list_group$Group1[[feature]]))) &&
      all(is.finite(as.numeric(beta_sigma_list_group$Group2[[feature]])))
  }, logical(1))
  
  beta_mu_mean_group1 <- beta_sigma_list_group$Group1[valid_features]
  beta_mu_mean_group2 <- beta_sigma_list_group$Group2[valid_features]
  
  ######## 6. Compute Integrals for Valid Features
  int_list <- bplapply(seq_along(beta_mu_mean_group1), function(i) {
    calculate_integral(
      as.numeric(beta_mu_mean_group1[[i]]),
      as.numeric(beta_mu_mean_group2[[i]])
    )
  }, BPPARAM = BPPARAM)
  
  ######## 7. Aggregate Results
  names(int_list) <- names(beta_mu_mean_group1)
  int_res <- sort(unlist(int_list), decreasing = TRUE)
  
  ######## 8. Return Results
  return(int_res)
}
