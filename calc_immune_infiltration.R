calc_immune_infiltration <- function(cohort, cohort_name,verbose=F){

require(dplyr)
require(tidyr)

    list_files_dm1 <- c("g_infiltration","g_tmb","g_msi")
    list_files_dm2 <- c("onco_result_immune_infiltration","onco_result_tmb_annotated","onco_result_msi_annotated")

      ############  Importing the tempus data and create the plotting dataframe
  if (is.character(cohort)) {
    if (verbose) {
      print("Trying to load cohort")
    }
      input_td <- tryCatch({load_tempus_data(cohort, collection=NULL,
                                   list_files=list_files_dm1)},
                 error=function(e){ tryCatch({load_tempus_data(cohort, collection=NULL,
                                                               list_files=list_files_dm2)
                 }, error=function(f){"dud_cohort"})})
      } else if (is.list(cohort)) {
    if (verbose) {
      print("Cohort already loaded, working")
    }
    input_td <- cohort
      }

  ####################


  data_model_check <- calc_data_model(input_td,
    list_tables_dm_1 = list_files_dm1,
    list_tables_dm_2 = list_files_dm2
  )

  stopifnot(data_model_check %in% c("1.0", "2.0"))
  ############ Define DM specific functions ##################
  prep_infiltration_dm1 <- function(input_td) {
    infiltration <- Reduce("left_join", lapply(input_td[list_files_dm1], function(x) unique(x[intersect(c("patient_id", "tmb_v1","percentile","msi_status","est_immune_cells","est_b_cells","est_cd4_cells","est_cd8_cells", "est_mac_cells","est_nk_cells"), colnames(x))])))


    return(infiltration)
  }

  prep_infiltration_dm2 <- function(input_td) {
      infiltration <- Reduce("left_join", lapply(input_td[list_files_dm2], function(x) unique(x[intersect(c("patient_id", "tmb_result_harmonized_forward_xt","msi_status","b_cells_fraction_of_immune_decimal","cd4_cells_fraction_of_immune_decimal","cd8_cells_fraction_of_immune_decimal","macrophages_fraction_of_immune_decimal","nk_cells_fraction_of_immune_decimal","immune_cells_fraction_of_isolate_decimal"), colnames(x))])))
      colnames(infiltration) <- gsub("_cells_fraction_of_immune_decimal|_fraction_of_immune_decimal|_fraction_of_isolate_decimal","_cells", colnames(infiltration))
      colnames(infiltration)[grep("_cells",colnames(infiltration))] <- paste0("est_",colnames(infiltration)[grep("_cells",colnames(infiltration))])
      colnames(infiltration)[grep("tmb_result",colnames(infiltration))] <- "tmb_v1"
      infiltration$percentile <- round(100*get_quantile(infiltration$tmb_v1))

   return(infiltration)
  }

  ##### Calling DM specific functions #####
  if (data_model_check == "1.0") {
    infiltration <- prep_infiltration_dm1(input_td)
  }

  if (data_model_check == "2.0") {
    infiltration <- prep_infiltration_dm2(input_td)
  }

    infiltration$Cohort <- cohort_name

    return(infiltration)
}
