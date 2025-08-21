calc_treatment_journeys <- function(cohort,treatment_lines = 5,num_groups=6,cohort_name=NULL,verbose=F, drop_untreated=TRUE){

    require(ggsankey)
require(dplyr)
require(tidyr)

    list_files_dm1 <- c("regimen")
    list_files_dm2 <- c("onco_regimen")

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
  prep_treatment_journeys_dm1 <- function(input_td) {
    regimen_table <- input_td$regimen
    regimen_summary <- regimen_table %>%
      dplyr::arrange(.data$patient_id, .data$regimen_rank) %>%
      dplyr::select(.data$patient_id,
        Rank = .data$regimen_rank, Drug_Name = .data$regimen_name, Drug_Class = .data$regimen_class,
        Class_Group = .data$regimen_class_group
      )
    return(regimen_summary)
  }

  prep_treatment_journeys_dm2 <- function(input_td) {
    regimen_table <- input_td$onco_regimen
    regimen_summary <- regimen_table %>%
      dplyr::arrange(.data$patient_id, .data$regimen_sequence) %>%
      dplyr::select(.data$patient_id,
        Rank = .data$regimen_sequence, Drug_Name = .data$agents, Drug_Class = .data$therapy_class,
        Class_Group = .data$therapy_class_group
      )
    return(regimen_summary)
  }

  ##### Calling DM specific functions #####
  if (data_model_check == "1.0") {
    regimen_summary <- prep_treatment_journeys_dm1(input_td)
  }

  if (data_model_check == "2.0") {
    regimen_summary <- prep_treatment_journeys_dm2(input_td)
  }
################################################
########## Format for output #################
################################################
  regimen_summary <- regimen_summary %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::mutate(Rank2 = paste0(1:length(.data$Rank), "L")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Rank2) %>%
    dplyr::mutate(Drug_Name2 = forcats::fct_lump(.data$Drug_Name,
      n = num_groups
    )) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Rank2) %>%
    dplyr::mutate(Class_Group2 = forcats::fct_lump(.data$Class_Group,
      n = num_groups
    )) %>%
    dplyr::ungroup()
  regimen_summary <- regimen_summary %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::arrange(.data$patient_id, .data$Rank2) %>%
    dplyr::ungroup() %>%
    data.frame()
  regimen_sankey <- regimen_summary %>%
    dplyr::select(.data$patient_id,
      Treatment_Line = .data$Rank2, Drug = .data$Drug_Name2, Class = .data$Class_Group2
    ) %>%
    dplyr::ungroup()
################################################
########## Finalize output #################
################################################

      regimen_sankey_class <- regimen_sankey %>% tidyr::pivot_wider(
    id_cols = .data$patient_id,
    names_from = .data$Treatment_Line, values_from = .data$Class
  )
  regimen_sankey_drug <- regimen_sankey %>% tidyr::pivot_wider(
    id_cols = .data$patient_id,
    names_from = .data$Treatment_Line, values_from = .data$Drug
  )

    all_lines <- unique(regimen_summary$Rank2)
  treatment_cols <- intersect(
    paste0(1:treatment_lines, "L"),
    all_lines
  )
  regimen_sankey_class_final <- regimen_sankey_class %>% ggsankey::make_long(treatment_cols)
    regimen_sankey_drug_final <- regimen_sankey_drug %>% ggsankey::make_long(treatment_cols)
      if (drop_untreated) {
    regimen_sankey_class_final <- regimen_sankey_class_final %>%
      tidyr::drop_na(.data$node)
    regimen_sankey_drug_final <- regimen_sankey_drug_final %>%
      tidyr::drop_na(.data$node)
  } else {
    regimen_sankey_class_final <- regimen_sankey_class_final %>%
      tidyr::replace_na(list(node = "No F/U"))
    regimen_sankey_drug_final <- regimen_sankey_drug_final %>%
      tidyr::replace_na(list(node = "No F/U"))
  }



    out <- list("Drug"=as.data.frame(regimen_sankey_drug_final),"Class"=as.data.frame(regimen_sankey_class_final))
    return(out)
    }
