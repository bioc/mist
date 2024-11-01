test_mist <- function() {
  set.seed(123)
  Dat_path <- system.file("extdata", "small_sampleData_sce.rds", package = "mist")
  beta_sigma_list <- estiParamSingle(Dat_sce = Dat_path,
                                     Dat_name = 'Methy_level_group1',
                                    ptime_name = 'pseudotime')
  dm_results <- dmSingle(beta_sigma_list)

  beta_sigma_list_group <- estiParamTwo(Dat_sce = Dat_path,
                                        Dat_name_g1 = 'Methy_level_group1',
                                        Dat_name_g2 = 'Methy_level_group2',
                                        ptime_name_g1 = 'pseudotime',
                                        ptime_name_g2 = 'pseudotime_g2')
  dm_results_two <- dmTwoGroups(beta_sigma_list_group)


  checkEquals(length(dm_results), 5)
  checkTrue(all(dm_results >= 0))
  checkTrue(all(names(beta_sigma_list_group) == c("Group1", "Group2")))
  checkEquals(length(dm_results_two), 5)
  checkTrue(all(dm_results_two >= 0))


}
