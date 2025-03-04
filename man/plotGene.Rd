% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGene.R
\name{plotGene}
\alias{plotGene}
\title{Plot Methylation Levels vs. Pseudotime for a Specific Gene}
\usage{
plotGene(Dat_sce, Dat_name, ptime_name, gene_name)
}
\arguments{
\item{Dat_sce}{The updated sce object with A numeric matrix of estimated parameters for all genomic features in the rowData, including:
\itemize{
\item \eqn{\beta_0} to \eqn{\beta_4}: Estimated coefficients for the polynomial of degree 4.
\item \eqn{\sigma^2_1} to \eqn{\sigma^2_4}: Estimated variances for each stage along the pseudotime.
}}

\item{Dat_name}{A character string specifying the name of the assay to extract the methylation data.}

\item{ptime_name}{A character string specifying the name of the colData to extract the pseudotime vector.}

\item{gene_name}{A character string specifying the gene name to plot.}
}
\value{
A ggplot2 scatter plot with an overlayed fitted curve.
}
\description{
Generates a scatter plot of methylation levels vs. pseudotime for a specific gene,
overlaid with the fitted curve based on the estimated coefficients.
}
\examples{
library(SingleCellExperiment)
library(ggplot2)
data <- readRDS(system.file("extdata", "group1_sampleData_sce.rds", package = "mist"))
Dat_sce_new <- estiParam(
    Dat_sce = data,
    Dat_name = "Methy_level_group1",
    ptime_name = "pseudotime"
)
plotGene(Dat_sce = Dat_sce_new,
Dat_name = "Methy_level_group1",
ptime_name = "pseudotime",
gene_name = "ENSMUSG00000000037")
}
