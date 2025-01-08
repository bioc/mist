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
dmSingle <- function(Dat_sce) {
    ######## 1. Input Validation
  # Check if Dat_sce is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }

    ######## 2. Remove Genomic Features with Invalid Values
    # Extract the matrix from rowData
    mist_pars_matrix <- rowData(Dat_sce)$mist_pars
  
    # Convert the matrix back to a list
    beta_sigma_list <- split(as.data.frame(mist_pars_matrix), rownames(mist_pars_matrix))
  
    # Convert each row back to a numeric vector
    beta_sigma_list <- lapply(beta_sigma_list, as.numeric)
  
    # Identify and remove genomic features with NA or infinite values
    bad_elements <- lapply(beta_sigma_list, contains_inf_or_na) # contains_inf_or_na checks for NA/Inf values
    valid_features <- !unlist(bad_elements) # Features without NA/Inf values
    beta_mu_mean <- beta_sigma_list[valid_features] # Keep only valid features

    ######## 3. Integral Calculation for Each Genomic Feature
    # Parallelized computation of integrals for each genomic feature
    int_list <- bplapply(beta_mu_mean, calculate_integral,
                         BPPARAM = SnowParam())
    # Assign names to the integrals based on the genomic feature names
    names(int_list) <- names(beta_mu_mean)

    ######## 4. Convert to Named Numeric Vector and Sort
    # Convert the list of integrals to a named numeric vector
    int_res <- unlist(int_list)
    # Sort the integrals in descending order (larger values indicate more drastic changes)
    int_res <- sort(int_res, decreasing = TRUE)
    rowData(Dat_sce)$mist_int <- int_res
    ######## 5. Return the Result
    return(Dat_sce)
}
