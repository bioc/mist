#' Parameter Estimation for Single-Group
#'
#' This function performs the Gibbs sampling procedure based on hierarchical Bayesian modeling
#' to produce the parameters required for differential methylation analysis.
#'
#' @param Dat_sce A `SingleCellExperiment` object containing the single-cell DNA methylation level.
#'   Methylation levels should be stored as an assay, with genomic feature (gene) names in rownames
#'   and cells in colnames.
#' @param Dat_name A character string specifying the name of the assay to extract the methylation level data.
#' @param ptime_name A character string specifying the name of the column in `colData` containing the pseudotime vector.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend for computations, as used in `bplapply()`. Defaults to `SnowParam()` for cluster-based parallel processing.
#' @param verbose A logical value indicating whether to print progress messages to the console.
#'   Defaults to \code{TRUE}. Set to \code{FALSE} to suppress messages.
#'   
#' @return The updated sce object with A numeric matrix of estimated parameters for all genomic features in the rowData, including:
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
#' Dat_sce_new <- estiParamSingle(
#'     Dat_sce = data,
#'     Dat_name = "Methy_level_group1",
#'     ptime_name = "pseudotime"
#' )
estiParamSingle <- function(Dat_sce,
                            Dat_name,
                            ptime_name,
                            BPPARAM = SnowParam(),
                            verbose = TRUE) {
  ######## 1. Input Validation
  # Check if Dat_sce is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }

  # Check if Dat_name is provided
  if (is.null(Dat_name)) {
    stop("Missing Dat_name: Specify the assay name to extract data.",
         call. = TRUE, domain = NULL)
  }

  # Check if ptime_name is provided
  if (is.null(ptime_name)) {
    stop("Missing ptime_name: Specify the column name in colData for pseudotime.",
         call. = TRUE, domain = NULL)
  }

  # Extract the assay and pseudotime data
  scDNAm_mat <- assay(Dat_sce, Dat_name)
  ptime <- colData(Dat_sce)[[ptime_name]]

  # Check if the number of cells (columns) matches the length of pseudotime
  if (ncol(scDNAm_mat) != length(ptime)) {
    stop("The number of cells in the data matrix and the pseudotime vector must match.",
         call. = TRUE, domain = NULL)
  }

  ###### 2. Normalize pseudotime to 0 - 1
  ptime <- ptime[is.finite(ptime) & !is.na(ptime)]
  ptime_all <- c(ptime / max(ptime))
  if (verbose) message("Pseudotime cleaning and normalization to [0, 1] completed.")

  ###### 3. Remove Genomic Features Containing Only 0/1 Values in All Timepoints
  rmRes <- BiocParallel::bplapply(scDNAm_mat, rmBad, ptime_all = ptime_all,
                                  BPPARAM = BPPARAM)
  rmIndex <- which(unlist(rmRes) == 1)
  scDNAm_mat_clean <- scDNAm_mat[!seq_len(nrow(scDNAm_mat)) %in% rmIndex, ]
  if (verbose) message("Removal of genomic features with too many 0/1 values completed.")

  ###### 4. Parameter Estimation using Gibbs Sampling
  beta_sigma_list <- BiocParallel::bplapply(seq_len(nrow(scDNAm_mat_clean)), run_bayesian_estimation,
                                            dat_ready = scDNAm_mat_clean, ptime_all = ptime_all,
                                            BPPARAM = BPPARAM)

  # Assign names to the parameters
  name_vector <- c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")
  beta_sigma_list <- lapply(beta_sigma_list, function(x) {
    names(x) <- name_vector
    return(x)
  })

  # Assign genomic feature names to the list
  names(beta_sigma_list) <- rownames(scDNAm_mat_clean)

  rowData(Dat_sce)$mist_pars <- do.call(rbind, beta_sigma_list)
  
  # Return the final sce object
  return(Dat_sce)
}
