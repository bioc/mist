test_mist <- function() {
  set.seed(123)
  Dat_sce <- readRDS(system.file("extdata", "small_sampleData_sce.rds", package = "mist"))
  # beta_sigma_list <- estiParamSingle(Dat_sce = Dat_sce,
  #                                    Dat_name = 'Methy_level_group1',
  #                                   ptime_name = 'pseudotime')
  # dm_results <- dmSingle(beta_sigma_list)

  Dat_sce <- estiParamSingle(Dat_sce = Dat_sce,
                             Dat_name = "Methy_level_group1",
                             ptime_name = "pseudotime")
  Dat_sce <- dmSingle(Dat_sce)


  #checkEquals(length(dm_results), 5)
  #checkTrue(all(dm_results >= 0))
  checkEquals(nrow(rowData(Dat_sce)$mist_pars), 5)
  checkEquals(length(rowData(Dat_sce)$mist_int), 5)
  checkTrue(all(rowData(Dat_sce)$mist_int >= 0))


}
