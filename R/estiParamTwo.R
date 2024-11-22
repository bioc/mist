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
#'
#' @return A list containing two elements, one for each group, where each element is a numeric list of estimated
#'   parameters for all genomic features. Each element includes:
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
#' beta_sigma_list_group <- estiParamTwo(
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
                         ptime_name_g2) {
  ######## 1. Input Validation
  # Check if Dat_sce is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }

  # Check if Dat_name_g1 and Dat_name_g2 are provided
  if (is.null(Dat_name_g1)) {
    stop("Missing Dat_name_g1: Specify the assay name for group 1.",
         call. = TRUE, domain = NULL)
  }
  if (is.null(Dat_name_g2)) {
    stop("Missing Dat_name_g2: Specify the assay name for group 2.",
         call. = TRUE, domain = NULL)
  }

  # Check if ptime_name_g1 and ptime_name_g2 are provided
  if (is.null(ptime_name_g1)) {
    stop("Missing ptime_name_g1: Specify the column name in colData for group 1 pseudotime.",
         call. = TRUE, domain = NULL)
  }
  if (is.null(ptime_name_g2)) {
    stop("Missing ptime_name_g2: Specify the column name in colData for group 2 pseudotime.",
         call. = TRUE, domain = NULL)
  }

  # Extract assays and pseudotime data
  scDNAm_mat_group1 <- assay(Dat_sce, Dat_name_g1)
  scDNAm_mat_group2 <- assay(Dat_sce, Dat_name_g2)
  ptime_group1 <- colData(Dat_sce)[[ptime_name_g1]]
  ptime_group2 <- colData(Dat_sce)[[ptime_name_g2]]

  # Check for matching dimensions
  if (ncol(scDNAm_mat_group1) != length(ptime_group1)) {
    stop("The number of cells in the data matrix and pseudotime vector for group 1 must match.",
         call. = TRUE, domain = NULL)
  }
  if (ncol(scDNAm_mat_group2) != length(ptime_group2)) {
    stop("The number of cells in the data matrix and pseudotime vector for group 2 must match.",
         call. = TRUE, domain = NULL)
  }

  ######## 2. Normalize Pseudotime
  ptime_group1 <- ptime_group1[is.finite(ptime_group1) & !is.na(ptime_group1)]
  ptime_group2 <- ptime_group2[is.finite(ptime_group2) & !is.na(ptime_group2)]

  ptime_group1_all <- ptime_group1 / max(ptime_group1)
  ptime_group2_all <- ptime_group2 / max(ptime_group2)
  message("Pseudotime cleaning and normalization to [0, 1] completed.")

  ######## 3. Remove Genomic Features Containing Only 0/1 Values
  rmRes_group1 <- BiocParallel::bplapply(scDNAm_mat_group1, rmBad, ptime_all = ptime_group1_all)
  rmRes_group2 <- BiocParallel::bplapply(scDNAm_mat_group2, rmBad, ptime_all = ptime_group2_all)

  rmIndex_group1 <- which(unlist(rmRes_group1) == 1)
  rmIndex_group2 <- which(unlist(rmRes_group2) == 1)

  scDNAm_mat_clean_group1 <- scDNAm_mat_group1[!seq_len(nrow(scDNAm_mat_group1)) %in% rmIndex_group1, ]
  scDNAm_mat_clean_group2 <- scDNAm_mat_group2[!seq_len(nrow(scDNAm_mat_group2)) %in% rmIndex_group2, ]
  message("Removal of genomic features with too many 0/1 values completed.")

  ######## 4. Parameter Estimation using Gibbs Sampling
  beta_sigma_list_group1 <- BiocParallel::bplapply(seq_len(nrow(scDNAm_mat_clean_group1)),
                                                   run_bayesian_estimation,
                                                   dat_ready = scDNAm_mat_clean_group1,
                                                   ptime_all = ptime_group1_all)

  beta_sigma_list_group2 <- BiocParallel::bplapply(seq_len(nrow(scDNAm_mat_clean_group2)),
                                                   run_bayesian_estimation,
                                                   dat_ready = scDNAm_mat_clean_group2,
                                                   ptime_all = ptime_group2_all)

  ######## 5. Name the Results
  name_vector <- c("Beta_0", "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4")

  beta_sigma_list_group1 <- lapply(beta_sigma_list_group1, function(x) {
    names(x) <- name_vector
    return(x)
  })
  beta_sigma_list_group2 <- lapply(beta_sigma_list_group2, function(x) {
    names(x) <- name_vector
    return(x)
  })

  ######## 6. Assign Genomic Feature Names
  names(beta_sigma_list_group1) <- rownames(scDNAm_mat_clean_group1)
  names(beta_sigma_list_group2) <- rownames(scDNAm_mat_clean_group2)

  ######## 7. Return the Results
  return(list(Group1 = beta_sigma_list_group1, Group2 = beta_sigma_list_group2))
}
