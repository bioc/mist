#### Parameter Estimation function for two-group scenario
#' Parameter estimation for scDNAm data under the 2-group scenario.
#'
#' This function performs the Gibbs sampling procedure based on hierarchical Bayesian modeling,
#' separately for two groups, to produce the parameters required in the differential methylation analysis.
#'
#' @param scDNAm_mat_group1 A numeric matrix containing the single-cell DNA methylation levels for Group 1.
#' Genomic feature (gene) names are labeled in rownames, and cells are labeled in colnames.
#' @param scDNAm_mat_group2 A numeric matrix containing the single-cell DNA methylation levels for Group 2.
#' Genomic feature (gene) names are labeled in rownames, and cells are labeled in colnames.
#' @param ptime_group1 A numeric vector containing the estimated pseudotime of all cells in Group 1.
#' @param ptime_group2 A numeric vector containing the estimated pseudotime of all cells in Group 2.
#'
#' @return A list containing two elements, one for each group, where each element is a numeric list of estimated
#' parameters for all genomic features. Each element includes:
#' - \eqn{\beta_0}, \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\beta_4}: The estimated coefficients for
#'   the polynomial of degree 4.
#' - \eqn{\sigma^2_1}, \eqn{\sigma^2_2}, \eqn{\sigma^2_3}, \eqn{\sigma^2_4}: The estimated variances for each
#'   stage along the pseudotime.
#'
#' @import stats MCMCpack BiocParallel car mvtnorm
#' @importFrom S4Vectors subjectHits queryHits Rle
#' @importFrom methods is
#' @importFrom rtracklayer start offset end
#' @importFrom Matrix cov2cor toeplitz update
#' @export
estiParamTwo <- function(scDNAm_mat_group1,
                               scDNAm_mat_group2,
                               ptime_group1,
                               ptime_group2) {
  # 1. Input Validation
  if (is.null(scDNAm_mat_group1) || is.null(scDNAm_mat_group2)) {
    stop("Missing Data Matrix for one of the groups!", call. = TRUE)
  }

  if (is.null(ptime_group1) || is.null(ptime_group2)) {
    stop("Missing Pseudotime Vector for one of the groups!", call. = TRUE)
  }

  if (ncol(scDNAm_mat_group1) != length(ptime_group1)) {
    stop("Cell Numbers of Data Matrix and Pseudotime Vector for Group 1 are not the same!", call. = TRUE)
  }

  if (ncol(scDNAm_mat_group2) != length(ptime_group2)) {
    stop("Cell Numbers of Data Matrix and Pseudotime Vector for Group 2 are not the same!", call. = TRUE)
  }

  # 2. Normalize pseudotime to 0-1 for both groups
  ptime_group1 <- ptime_group1[is.finite(ptime_group1)&!is.na(ptime_group1)]
  ptime_group2 <- ptime_group2[is.finite(ptime_group2)&!is.na(ptime_group2)]

  ptime_group1_all <- ptime_group1 / max(ptime_group1)
  ptime_group2_all <- ptime_group2 / max(ptime_group2)
  message("Pseudotime cleaning and normalization to [0, 1] completed.")
  # 3. Remove Genomic Features with only 0/1 values in both groups
  rmRes_group1 <- BiocParallel::bplapply(scDNAm_mat_group1, rmBad, ptime_all = ptime_group1_all)
  rmRes_group2 <- BiocParallel::bplapply(scDNAm_mat_group2, rmBad, ptime_all = ptime_group2_all)

  rmIndex_group1 <- which(unlist(rmRes_group1) == 1)
  rmIndex_group2 <- which(unlist(rmRes_group2) == 1)

  scDNAm_mat_clean_group1 <- scDNAm_mat_group1[!1:nrow(scDNAm_mat_group1) %in% rmIndex_group1, ]
  scDNAm_mat_clean_group2 <- scDNAm_mat_group2[!1:nrow(scDNAm_mat_group2) %in% rmIndex_group2, ]
  message("Removal of genomic features with too many 0/1 values completed.")
  # 4. Parameter Estimation using Gibbs Sampling for both groups
  beta_sigma_list_group1 <- bplapply(1:nrow(scDNAm_mat_clean_group1), run_bayesian_estimation,
                                     dat_ready = scDNAm_mat_clean_group1,
                                     ptime_all = ptime_group1_all)

  beta_sigma_list_group2 <- bplapply(1:nrow(scDNAm_mat_clean_group2), run_bayesian_estimation,
                                     dat_ready = scDNAm_mat_clean_group2,
                                     ptime_all = ptime_group2_all)

  # 5. Name the results for both groups
  name_vector <- c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")

  beta_sigma_list_group1 <- lapply(beta_sigma_list_group1, function(x) {
    names(x) <- name_vector
    return(x)
  })

  beta_sigma_list_group2 <- lapply(beta_sigma_list_group2, function(x) {
    names(x) <- name_vector
    return(x)
  })

  # 6. Assign genomic feature names
  names(beta_sigma_list_group1) <- rownames(scDNAm_mat_clean_group1)
  names(beta_sigma_list_group2) <- rownames(scDNAm_mat_clean_group2)

  # 7. Return a list containing both groups' parameter estimates
  return(list(Group1 = beta_sigma_list_group1, Group2 = beta_sigma_list_group2))
}
