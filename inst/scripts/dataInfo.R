# This script documents how the data in inst/extdata was generated.
#
# Data Source:
# The data originates from GEO dataset GSE121708, described in:
# 1. Argelaguet R, Clark SJ, Mohammed H, Stapel LC, et al. Multi-omics profiling of mouse gastrulation at single-cell resolution. Nature 2019 Dec;576(7787):487-491. PMID: 31827285
# 2. Kapourani CA, Argelaguet R, Sanguinetti G, Vallejos CA. scMET: Bayesian modeling of DNA methylation heterogeneity at single-cell resolution. Genome Biol 2021 Apr 20;22(1):114. PMID: 33879195
#
# Dataset Description:
# - The dataset contains information on cell developmental stages during mouse gastrulation.
# - Raw data was downloaded from GEO and annotated using the mm10 genome reference.
#
# Data Processing Steps:
# 1. Annotation:
#    - Data was annotated using the mm10 reference genome.
# 2. Data Cleaning:
#    - Genes expressed in less than 10% of cells were removed.
#    - This resulted in a scDNAm dataset containing 986 cells and 18,220 genes.
# 3. Pseudotime Inference:
#    - Pseudotime was inferred using the default settings in Monocle3.
# 4. Subsampling for Illustration and Testing:
#    - For illustration purposes (e.g., vignette generation), random samples of 5 genes were selected to create `group1_sampleData_sce.rds` and `group2_sampleData_sce.rds`.
#    - Note: The sampling was random without using a seed, so results may not be reproducible.

