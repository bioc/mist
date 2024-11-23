#' Plot Methylation Levels vs. Pseudotime for a Specific Gene
#'
#' Generates a scatter plot of methylation levels vs. pseudotime for a specific gene,
#' overlaid with the fitted curve based on the estimated coefficients.
#'
#' @param Dat_sce A `SingleCellExperiment` object containing the single-cell DNA methylation level.
#'   Methylation levels should be stored as assays, with genomic feature (gene) names in rownames
#'   and cells in colnames.
#' @param Dat_name A character string specifying the name of the assay to extract the methylation data.
#' @param ptime_name A character string specifying the name of the colData to extract the pseudotime vector.
#' @param beta_sigma_list A named numeric list of estimated coefficients (output of `estiParamSingle`).
#' @param gene_name A character string specifying the gene name to plot.
#'
#' @return A ggplot2 scatter plot with an overlayed fitted curve.
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' plotGene(Dat_sce = Dat_sce,
#' Dat_name = "Methy_level_group1",
#' ptime_name = "pseudotime",
#' beta_sigma_list,
#' gene_name = "ENSMUSG00000000037")
#' }
plotGene <- function(Dat_sce, Dat_name, ptime_name, beta_sigma_list, gene_name) {
  # Check if Dat_sce is a SingleCellExperiment object
  if (!methods::is(Dat_sce, "SingleCellExperiment")) {
    stop("Dat_sce must be a SingleCellExperiment object.",
         call. = TRUE, domain = NULL)
  }

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
  coefficients <- beta_sigma_list[[gene_name]][1:5] # Beta_0 to Beta_4

  # Generate the fitted values using the degree-4 polynomial and inverse logit transformation
  fitted_values <- sapply(pseudotime, function(t) {
    z_t <- c(1, t, t^2, t^3, t^4) # Polynomial terms (include intercept)
    linear_predictor <- sum(z_t * coefficients) # Degree-4 polynomial
    1 / (1 + exp(-linear_predictor)) # Inverse logit transformation
  })

  # Create the scatter plot with the fitted curve
  library(ggplot2)
  plot_data <- data.frame(
    Pseudotime = pseudotime,
    Methylation = methylation_levels,
    Fitted = fitted_values
  )

  ggplot(plot_data, aes(x = Pseudotime)) +
    geom_point(aes(y = Methylation), color = "blue", alpha = 0.7, size = 2) +
    geom_line(aes(y = Fitted), color = "red", linewidth = 1) +
    labs(
      title = paste("Methylation Levels vs. Pseudotime for Gene:", gene_name),
      x = "Pseudotime",
      y = "Methylation Level"
    ) +
    theme_classic()
}
