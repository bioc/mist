#### Parameter Estimation function for two-group scenario
#' Parameter estimation for scDNAm data under the 2-group scenario.
#'
#' This function performs the Gibbs sampling procedure based on hierarchical Bayesian modeling,
#' separately for two groups, to produce the parameters required in the differential methylation analysis.
#'
#' @param Dat_sce A character string specifying the file path to the SingleCellExperiment object containing
#' the single-cell DNA methylation level. Methylation levels can be calculated as the ratio of methylation
#' coverage over the total coverage. Genomic feature (gene) names are labeled in rownames, and cells are
#' labeled in colnames.
#' @param Dat_name_g1 A character string specifying the name of the assay to extract the group 1 methylation
#' level data from the SCE obejct. Data matrix should be stored in the 'Assay'.
#' @param Dat_name_g2 A character string specifying the name of the assay to extract the group 2 methylation
#' level data from the SCE obejct. Data matrix should be stored in the 'Assay'.
#' @param ptime_name_g1 A character string specifying the name of the colData to extract the group 1 pseudotime
#' vector from the SCE obejct. Pseudotime should be stored in the 'colData'.
#' @param ptime_name_g2 A character string specifying the name of the colData to extract the group 2 pseudotime
#' vector from the SCE obejct. Pseudotime should be stored in the 'colData'.
#'
#' @return A list containing two elements, one for each group, where each element is a numeric list of estimated
#' parameters for all genomic features. Each element includes:
#' - \eqn{\beta_0}, \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\beta_4}: The estimated coefficients for
#'   the polynomial of degree 4.
#' - \eqn{\sigma^2_1}, \eqn{\sigma^2_2}, \eqn{\sigma^2_3}, \eqn{\sigma^2_4}: The estimated variances for each
#'   stage along the pseudotime.
#'
#' @import MCMCpack BiocParallel car mvtnorm SummarizedExperiment SingleCellExperiment
#' @importFrom S4Vectors subjectHits queryHits Rle
#' @importFrom methods is
#' @importFrom rtracklayer start offset end
#' @importFrom Matrix cov2cor toeplitz update
#' @importFrom stats pgamma poly qgamma rnorm runif
#' @export
#'
#' @examples
#' library(mist)
#' Dat_path <- system.file("extdata", "small_sampleData_sce.rds", package = "mist")
#' beta_sigma_list_group <- estiParamTwo(
#'     Dat_sce = Dat_path,
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
    # 1. Input Validation
    # Check if Dat_sce is provided
    if (is.null(Dat_sce)) {
        stop("Missing Data Object!",
            call. = TRUE, domain = NULL
        )
    }

    # Check if Dat_name is provided
    if (is.null(Dat_name_g1)) {
        stop("Missing group 1 Data Name!",
            call. = TRUE, domain = NULL
        )
    }
    if (is.null(Dat_name_g2)) {
        stop("Missing group 2 Data Name!",
            call. = TRUE, domain = NULL
        )
    }

    # Check if ptime_name is provided
    if (is.null(ptime_name_g1)) {
        stop("Missing group 1 Pseudotime Name!",
            call. = TRUE, domain = NULL
        )
    }
    if (is.null(ptime_name_g2)) {
        stop("Missing group 2 Pseudotime Name!",
            call. = TRUE, domain = NULL
        )
    }

    Data_sce <- readRDS(file = Dat_sce)
    scDNAm_mat_group1 <- assay(Data_sce, Dat_name_g1)
    scDNAm_mat_group2 <- assay(Data_sce, Dat_name_g2)
    ptime_group1 <- colData(Data_sce)[[ptime_name_g1]]
    ptime_group2 <- colData(Data_sce)[[ptime_name_g2]]

    if (ncol(scDNAm_mat_group1) != length(ptime_group1)) {
        stop("Cell Numbers of Data Matrix and Pseudotime Vector for Group 1 are not the same!", call. = TRUE)
    }

    if (ncol(scDNAm_mat_group2) != length(ptime_group2)) {
        stop("Cell Numbers of Data Matrix and Pseudotime Vector for Group 2 are not the same!", call. = TRUE)
    }

    # 2. Normalize pseudotime to 0-1 for both groups
    ptime_group1 <- ptime_group1[is.finite(ptime_group1) & !is.na(ptime_group1)]
    ptime_group2 <- ptime_group2[is.finite(ptime_group2) & !is.na(ptime_group2)]

    ptime_group1_all <- ptime_group1 / max(ptime_group1)
    ptime_group2_all <- ptime_group2 / max(ptime_group2)
    message("Pseudotime cleaning and normalization to [0, 1] completed.")
    # 3. Remove Genomic Features with only 0/1 values in both groups
    rmRes_group1 <- BiocParallel::bplapply(scDNAm_mat_group1, rmBad, ptime_all = ptime_group1_all)
    rmRes_group2 <- BiocParallel::bplapply(scDNAm_mat_group2, rmBad, ptime_all = ptime_group2_all)

    rmIndex_group1 <- which(unlist(rmRes_group1) == 1)
    rmIndex_group2 <- which(unlist(rmRes_group2) == 1)

    scDNAm_mat_clean_group1 <- scDNAm_mat_group1[!seq_len(nrow(scDNAm_mat_group1)) %in% rmIndex_group1, ]
    scDNAm_mat_clean_group2 <- scDNAm_mat_group2[!seq_len(nrow(scDNAm_mat_group2)) %in% rmIndex_group2, ]
    message("Removal of genomic features with too many 0/1 values completed.")
    # 4. Parameter Estimation using Gibbs Sampling for both groups
    beta_sigma_list_group1 <- bplapply(seq_len(nrow(scDNAm_mat_clean_group1)), run_bayesian_estimation,
        dat_ready = scDNAm_mat_clean_group1,
        ptime_all = ptime_group1_all
    )

    beta_sigma_list_group2 <- bplapply(seq_len(nrow(scDNAm_mat_clean_group2)), run_bayesian_estimation,
        dat_ready = scDNAm_mat_clean_group2,
        ptime_all = ptime_group2_all
    )

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
