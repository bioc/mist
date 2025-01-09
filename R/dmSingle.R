#' Differential methylation evaluation over time for scDNA-seq data under the 1-group scenario.
#'
#' This function performs DM analysis to identify genomic features showing drastic
#' changes along pseudotime, by modeling the methylation changes along pseudotime and estimating
#' the developmental-stage-specific biological variations.
#'
#' @param Dat_sce The updated sce object with A numeric matrix of estimated parameters for all genomic features in the rowData, including:
#'   - \eqn{\beta_0} to \eqn{\beta_4}: Estimated coefficients for the polynomial of degree 4.
#'   - \eqn{\sigma^2_1} to \eqn{\sigma^2_4}: Estimated variances for each stage along the pseudotime.
#' 
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend for computations, as used in `bplapply()`. Defaults to `SnowParam()` for cluster-based parallel processing.
#'   
#' @return The updated sce object with A named numeric vector where each value corresponds to a genomic feature
#' (e.g., a gene) in the rowData. The values represent the minimum area between the fitted curve
#' of scDNA methylation levels along pseudotime and a constant horizontal line for
#' each feature. This metric measures the magnitude of methylation level changes
#' along pseudotime. The values are sorted in descending order, with larger values
#' indicating more drastic changes.
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
#' library(mist)
#' data <- readRDS(system.file("extdata", "small_sampleData_sce.rds", package = "mist"))
#' Dat_sce_new <- estiParamSingle(
#'     Dat_sce = data,
#'     Dat_name = "Methy_level_group1",
#'     ptime_name = "pseudotime"
#' )
#' dm_sce <- dmSingle(Dat_sce_new)
dmSingle <- function(Dat_sce,
                     BPPARAM = SnowParam()) {
  ######## 1. Input Validation
  # Ensure input is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }
  
  ######## 2. Remove Genomic Features with Invalid Values
  # Extract the parameter matrix from rowData
  mist_pars_matrix <- rowData(Dat_sce)$mist_pars
  
  # Split matrix rows into a list of numeric vectors
  beta_sigma_list <- lapply(split(as.data.frame(mist_pars_matrix), rownames(mist_pars_matrix)), as.numeric)
  
  # Identify valid features (no NA/Inf values) in a single step
  valid_features <- vapply(beta_sigma_list, function(x) all(is.finite(x)), logical(1))
  
  # Subset valid features only
  beta_mu_mean <- beta_sigma_list[valid_features]
  
  ######## 3. Integral Calculation for Each Genomic Feature
  # Compute integrals in parallel with efficient chunking
  int_list <- bplapply(beta_mu_mean, calculate_integral, BPPARAM = BPPARAM)
  
  # Assign names of valid features to the result
  names(int_list) <- names(beta_mu_mean)
  
  ######## 4. Convert to Named Numeric Vector and Sort
  # Flatten the list of integrals and sort in descending order
  int_res <- sort(unlist(int_list), decreasing = TRUE)
  
  ######## 5. Save Results and Return
  rowData(Dat_sce)$mist_int <- int_res
  return(Dat_sce)
}
