#' Plot Methylation Levels vs. Pseudotime for a Specific Gene
#'
#' Generates a scatter plot of methylation levels vs. pseudotime for a specific gene,
#' overlaid with the fitted curve based on the estimated coefficients.
#'
#' @param Dat_sce The updated sce object with A numeric matrix of estimated parameters for all genomic features in the rowData, including:
#'   - \eqn{\beta_0} to \eqn{\beta_4}: Estimated coefficients for the polynomial of degree 4.
#'   - \eqn{\sigma^2_1} to \eqn{\sigma^2_4}: Estimated variances for each stage along the pseudotime.
#' @param Dat_name A character string specifying the name of the assay to extract the methylation data.
#' @param ptime_name A character string specifying the name of the colData to extract the pseudotime vector.
#' @param gene_name A character string specifying the gene name to plot.
#'
#' @return A ggplot2 scatter plot with an overlayed fitted curve.
#' @import rlang
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' library(ggplot2)
#' data <- readRDS(system.file("extdata", "small_sampleData_sce.rds", package = "mist"))
#' Dat_sce_new <- estiParamSingle(
#'     Dat_sce = data,
#'     Dat_name = "Methy_level_group1",
#'     ptime_name = "pseudotime"
#' )
#' plotGene(Dat_sce = Dat_sce_new,
#' Dat_name = "Methy_level_group1",
#' ptime_name = "pseudotime",
#' gene_name = "ENSMUSG00000000037")
plotGene <- function(Dat_sce, Dat_name, ptime_name, gene_name) {
  # Check if ggplot2 is available
  if (!requireNamespace('ggplot2', quietly = TRUE))
    stop("Install 'ggplot2' to use this function.")
  # Check if Dat_sce is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }
  mist_pars_matrix <- rowData(Dat_sce)$mist_pars
  
  # Convert the matrix back to a list
  beta_sigma_list <- split(as.data.frame(mist_pars_matrix), rownames(mist_pars_matrix))
  
  # Convert each row back to a numeric vector
  beta_sigma_list <- lapply(beta_sigma_list, as.numeric)
  # Extract scDNAm matrix and pseudotime
  scDNAm_mat <- assay(Dat_sce, Dat_name)
  ptime <- colData(Dat_sce)[[ptime_name]]

  # Check if the gene exists in the methylation matrix
  if (!gene_name %in% rownames(scDNAm_mat)) {
    stop("The specified gene name is not found in the scDNAm_mat.")
  }

  # Extract methylation levels and pseudotime for the specified gene
  methylation_levels <- scDNAm_mat[gene_name, ]
  pseudotime <- ptime

  # Extract the coefficients for the specified gene
  if (!gene_name %in% names(beta_sigma_list)) {
    stop("The specified gene name is not found in the beta_sigma_list.")
  }
  coefficients <- beta_sigma_list[[gene_name]][seq_len(5)] # Beta_0 to Beta_4

  # Generate the fitted values using the degree-4 polynomial and inverse logit transformation
  fitted_values <- vapply(pseudotime, function(t) {
    z_t <- c(1, t, t^2, t^3, t^4) # Polynomial terms (include intercept)
    linear_predictor <- sum(z_t * coefficients) # Degree-4 polynomial
    1 / (1 + exp(-linear_predictor)) # Inverse logit transformation
  }, numeric(1))

  # Create the scatter plot with the fitted curve
  plot_data <- data.frame(
    Pseudotime = pseudotime,
    Methylation = methylation_levels,
    Fitted = fitted_values
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Pseudotime)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$Methylation), color = "blue", alpha = 0.7, size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$Fitted), color = "red", linewidth = 1) +
    ggplot2::labs(
      title = paste("Methylation Levels vs. Pseudotime for Gene:", gene_name),
      x = "Pseudotime",
      y = "Methylation Level"
    ) +
    ggplot2::theme_classic()
}
