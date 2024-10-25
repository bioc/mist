#### Parameter Estimation function for single-group
#' Parameter estimation for scDNAm data under the 1-group scenario.
#'
#' This function performs the Gibbs sampling procedure based on the hierarchical Bayesian modeling,
#' to produce the parameters required in the differential methylation analysis.
#'
#' @param scDNAm_mat A numeric matrix containing the single-cell DNA methylation level. Methylation levels
#' can be calculated as the ratio of methylation coverage over the total coverage. Genomic feature (gene) names
#' are labeled in rownames, and cells are labeled in colnames.
#' @param ptime A numeric vector containing the estimated pseudotime of all cells in "scDNAm_mat".
#'
#' @return A numeric list of estimated parameters for all genomic features that will be used as the input
#' of the function \code{\link{dmSingle}}. Parameters include:
#' - \eqn{\beta_0}: The estimated intercept of the polynomial of degree 4.
#' - \eqn{\beta_1}: The estimated coefficient for the linear term of the polynomial.
#' - \eqn{\beta_2}: The estimated coefficient for the quadratic term of the polynomial.
#' - \eqn{\beta_3}: The estimated coefficient for the cubic term of the polynomial.
#' - \eqn{\beta_4}: The estimated coefficient for the quartic term of the polynomial.
#' - \eqn{\sigma^2_1}: The estimated variance for the first stage along the pseudotime.
#' - \eqn{\sigma^2_2}: The estimated variance for the second stage along the pseudotime.
#' - \eqn{\sigma^2_3}: The estimated variance for the third stage along the pseudotime.
#' - \eqn{\sigma^2_4}: The estimated variance for the fourth stage along the pseudotime.
#'
#' @import stats MCMCpack BiocParallel car mvtnorm
#' @importFrom S4Vectors subjectHits queryHits Rle
#' @importFrom methods is
#' @importFrom rtracklayer start offset end
#' @importFrom Matrix cov2cor toeplitz update
#' @export
estiParamSingle <- function(scDNAm_mat,
                            ptime) {

  ######## 1. Input Validation
  # Check if scDNAm_mat is provided
  if (is.null(scDNAm_mat)) {
    stop("Missing Data Matrix!",
         call. = TRUE, domain = NULL
    )
  }

  # Check if ptime is provided
  if (is.null(ptime)) {
    stop("Missing Pseudotime Vector!",
         call. = TRUE, domain = NULL
    )
  }

  # Check if the number of cells (columns) matches the length of pseudotime
  if (ncol(scDNAm_mat) != length(ptime)) {
    stop("Cell Numbers of Data Matrix and Pseudotime Vector are not the same!",
         call. = TRUE, domain = NULL
    )
  }

  ###### 2. Normalize pseudotime to 0 - 1
  # Normalize pseudotime to a scale between 0 and 1, keeping only finite and non-NA values
  ptime <- ptime[is.finite(ptime) & !is.na(ptime)]
  ptime_all <- c(ptime / max(ptime))  # Normalize to [0, 1]
  rm(ptime)  # Remove original pseudotime
  message("Pseudotime cleaning and normalization to [0, 1] completed.")
  ###### 3. Remove Genomic Features Containing Only 0/1 Values in All Timepoints
  # Remove rows (genomic features) with no variation (only 0s or 1s) across all time points
  rmRes <- BiocParallel::bplapply(scDNAm_mat, rmBad, ptime_all = ptime_all)
  rmIndex <- which(unlist(rmRes) == 1)
  scDNAm_mat_clean <- scDNAm_mat[!1:nrow(scDNAm_mat) %in% rmIndex, ]  # Keep only valid rows
  rm(scDNAm_mat)  # Remove the original matrix to save memory
  rm(rmRes)  # Remove temporary result
  message("Removal of genomic features with too many 0/1 values completed.")
  ###### 4. Parameter Estimation using Gibbs Sampling
  # Perform parameter estimation for each genomic feature (row) using Gibbs sampling

  beta_sigma_list <- bplapply(1:nrow(scDNAm_mat_clean), run_bayesian_estimation,
                              dat_ready = scDNAm_mat_clean, ptime_all = ptime_all)

  # Define names for the parameters
  name_vector <- c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")

  # Assign parameter names to each genomic feature's result
  beta_sigma_list <- lapply(beta_sigma_list, function(x) {
    names(x) <- name_vector
    return(x)
  })

  # Assign genomic feature names to the list
  names(beta_sigma_list) <- rownames(scDNAm_mat_clean)

  # Return the final list of estimated parameters for each genomic feature
  return(beta_sigma_list)
}
