test_mist <- function() {
  set.seed(123)
  Dat_sce <- readRDS(system.file("extdata", "small_sampleData_sce.rds", package = "mist"))
  # beta_sigma_list <- estiParamSingle(Dat_sce = Dat_sce,
  #                                    Dat_name = 'Methy_level_group1',
  #                                   ptime_name = 'pseudotime')
  # dm_results <- dmSingle(beta_sigma_list)

  Dat_sce <- estiParamTwo(Dat_sce = Dat_sce,
                                        Dat_name_g1 = 'Methy_level_group1',
                                        Dat_name_g2 = 'Methy_level_group2',
                                        ptime_name_g1 = 'pseudotime',
                                        ptime_name_g2 = 'pseudotime_g2')
  Dat_sce <- dmTwoGroups(Dat_sce)


  #checkEquals(length(dm_results), 5)
  #checkTrue(all(dm_results >= 0))
  checkEquals(nrow(rowData(Dat_sce)$mist_pars_group1), 5)
  checkEquals(nrow(rowData(Dat_sce)$mist_pars_group2), 5)
  checkEquals(length(rowData(Dat_sce)$mist_int_2group), 5)
  checkTrue(all(rowData(Dat_sce)$mist_int_2group >= 0))


}
