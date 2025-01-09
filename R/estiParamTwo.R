#' Parameter Estimation for Two-Group Scenario
#'
#' This function performs the Gibbs sampling procedure based on hierarchical Bayesian modeling,
#' separately for two groups, to produce the parameters required for differential methylation analysis.
#'
#' @param Dat_sce A `SingleCellExperiment` object containing the single-cell DNA methylation level.
#'   Methylation levels should be stored as assays, with genomic feature (gene) names in rownames
#'   and cells in colnames.
#' @param Dat_name_g1 A character string specifying the name of the assay to extract the group 1 methylation level data.
#' @param Dat_name_g2 A character string specifying the name of the assay to extract the group 2 methylation level data.
#' @param ptime_name_g1 A character string specifying the name of the column in `colData` for the group 1 pseudotime vector.
#' @param ptime_name_g2 A character string specifying the name of the column in `colData` for the group 2 pseudotime vector.
#' @param verbose A logical value indicating whether to print progress messages to the console.
#'   Defaults to \code{TRUE}. Set to \code{FALSE} to suppress messages.
#'   
#' @return The updated sce object with two matrices in the rowData, one for each group, where each matrix contains estimated
#'   parameters for all genomic features, including:
#'   - \eqn{\beta_0} to \eqn{\beta_4}: Estimated coefficients for the polynomial of degree 4.
#'   - \eqn{\sigma^2_1} to \eqn{\sigma^2_4}: Estimated variances for each stage along the pseudotime.
#'
#' @import MCMCpack BiocParallel car mvtnorm SummarizedExperiment SingleCellExperiment BiocGenerics
#' @importFrom stats pgamma poly qgamma rnorm runif
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data <- readRDS(system.file("extdata", "small_sampleData_sce.rds", package = "mist"))
#' Dat_sce_new <- estiParamTwo(
#'     Dat_sce = data,
#'     Dat_name_g1 = "Methy_level_group1",
#'     Dat_name_g2 = "Methy_level_group2",
#'     ptime_name_g1 = "pseudotime",
#'     ptime_name_g2 = "pseudotime_g2"
#' )
estiParamTwo <- function(Dat_sce,
                         Dat_name_g1,
                         Dat_name_g2,
                         ptime_name_g1,
                         ptime_name_g2,
                         verbose = TRUE) {
  ######## 1. Input Validation
  # Ensure input validity
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }
  
  if (is.null(Dat_name_g1)) stop("Missing Dat_name_g1: Specify the assay name for group 1.")
  if (is.null(Dat_name_g2)) stop("Missing Dat_name_g2: Specify the assay name for group 2.")
  if (is.null(ptime_name_g1)) stop("Missing ptime_name_g1: Specify pseudotime column for group 1.")
  if (is.null(ptime_name_g2)) stop("Missing ptime_name_g2: Specify pseudotime column for group 2.")
  
  # Extract assays and pseudotime data
  scDNAm_mat_group1 <- assay(Dat_sce, Dat_name_g1)
  scDNAm_mat_group2 <- assay(Dat_sce, Dat_name_g2)
  ptime_group1 <- colData(Dat_sce)[[ptime_name_g1]]
  ptime_group2 <- colData(Dat_sce)[[ptime_name_g2]]
  
  if (ncol(scDNAm_mat_group1) != length(ptime_group1)) stop("Mismatch in group 1 data dimensions.")
  if (ncol(scDNAm_mat_group2) != length(ptime_group2)) stop("Mismatch in group 2 data dimensions.")
  
  ######## 2. Normalize Pseudotime
  ptime_group1_all <- ptime_group1[is.finite(ptime_group1)] / max(ptime_group1, na.rm = TRUE)
  ptime_group2_all <- ptime_group2[is.finite(ptime_group2)] / max(ptime_group2, na.rm = TRUE)
  
  if (verbose) message("Pseudotime normalization completed.")
  
  ######## 3. Combined Parallel Processing for Filtering and Parameter Estimation
  combined_results <- bplapply(seq_len(nrow(scDNAm_mat_group1)), function(i) {
    # Process group 1
    group1_valid <- rmBad(scDNAm_mat_group1[i, ], ptime_all = ptime_group1_all)
    group1_estimation <- if (!group1_valid) {
      run_bayesian_estimation(i, dat_ready = scDNAm_mat_group1, ptime_all = ptime_group1_all)
    } else {
      NULL
    }
    
    # Process group 2
    group2_valid <- rmBad(scDNAm_mat_group2[i, ], ptime_all = ptime_group2_all)
    group2_estimation <- if (!group2_valid) {
      run_bayesian_estimation(i, dat_ready = scDNAm_mat_group2, ptime_all = ptime_group2_all)
    } else {
      NULL
    }
    
    list(
      group1 = if (!is.null(group1_estimation)) setNames(group1_estimation, c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")) else NULL,
      group2 = if (!is.null(group2_estimation)) setNames(group2_estimation, c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")) else NULL
    )
  }, BPPARAM = SnowParam())
  
  ######## 4. Post-Processing
  # Extract results for group 1 and group 2
  beta_sigma_list_group1 <- lapply(combined_results, `[[`, "group1")
  beta_sigma_list_group2 <- lapply(combined_results, `[[`, "group2")
  
  # Remove NULL entries (invalid genomic features)
  beta_sigma_list_group1 <- beta_sigma_list_group1[!sapply(beta_sigma_list_group1, is.null)]
  beta_sigma_list_group2 <- beta_sigma_list_group2[!sapply(beta_sigma_list_group2, is.null)]
  
  # Assign genomic feature names
  names(beta_sigma_list_group1) <- rownames(scDNAm_mat_group1)[!sapply(beta_sigma_list_group1, is.null)]
  names(beta_sigma_list_group2) <- rownames(scDNAm_mat_group2)[!sapply(beta_sigma_list_group2, is.null)]
  
  ######## 5. Assign Results to Dat_sce
  rowData(Dat_sce)$mist_pars_group1 <- do.call(rbind, beta_sigma_list_group1)
  rowData(Dat_sce)$mist_pars_group2 <- do.call(rbind, beta_sigma_list_group2)
  
  if (verbose) message("Parameter estimation completed.")
  return(Dat_sce)
}
