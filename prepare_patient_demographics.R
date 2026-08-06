#' Plot oncoprint of patients in Tempus data
#'
#' @param data_cohort Path to patient cohort of interest or the named list object
#' output by load_tempus_data. The path option is slower and not advised if multiple
#' oncoprints will be generated sequentially.
#' @param treatments_of_interest Vector of genes the user wants plotted (MANDATORY if
#' more than 50 genes detected in mmf)
#' @export
#' @examples
#' \dontrun{
#' prepare_patient_demographics(data_cohort,
#'   alteration_func = NULL,
#'   gene_identifier = "gene_canonical_name",
#'   title = NULL,
#'   filter_germline = TRUE,
#'   copy_number_threshold = 8,
#'   patient_order = NULL,
#'   patient_universe = NULL
#' )
#' }
#'
prepare_patient_demographics <- function(data_cohort, treatments_of_interest = c("capmatinib", "tepotinib", "crizotinib")) {
  demographics_list <- list()


  #### Metastatic or Locally advanced status
  metastatic_pts <- tempusr::calc_earliest_metastasis(data_cohort) |>
    select(.data$patient_id,
      stage_date = .data$earliest_mdx_date,
      stage_year = .data$earliest_mdx_year,
      stage_precision = .data$earliest_mdx_precision,
      stage_source = .data$mdx_source
    ) |>
    dplyr::mutate(stage_type = "metastatic")



  locally_advanced_pts <- data_cohort$onco_tumor_characterization |>
    dplyr::filter(
      # in the overall_stage column
      (.data$stage_system == "AJCC Stage" & .data$overall_stage %in% c("Stage 3A", "Stage 3B", "Stage 3C")) |
        # in the n_stage column
        stringr::str_detect(.data$n_stage, "N3"),
      # limit to known dates
      !is.na(.data$tumor_characterization_date),
      # exclude year only dates
      .data$tumor_characterization_date_precision != "YYYY",
      # exclude metastatic patients
      !.data$patient_id %in% metastatic_pts$patient_id
    ) |>
    # get earliest date per patient
    dplyr::group_by(.data$patient_id) |>
    dplyr::slice_min(.data$tumor_characterization_date,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    # select the columns of interest
    dplyr::select(.data$patient_id,
      stage_date = .data$tumor_characterization_date,
      stage_year = .data$tumor_characterization_date_year_indexed,
      stage_precision = .data$tumor_characterization_date_precision
    ) |>
    # denote which table the information is from
    dplyr::mutate(
      stage_date = as.Date(stage_date),
      stage_source = "onco_tumor_characterization",
      stage_type = "locally_advanced"
    )


  stage_info <- bind_rows(metastatic_pts, locally_advanced_pts) %>% dplyr::select(.data$patient_id, Stage = .data$stage_type)
  demographics_list$Stage <- stage_info



  #### Patient metadata

  demographics_list$Stage_Detailed <- dplyr::select(input_td$onco_tumor_characterization, .data$patient_id,.data$stage_system, .data$overall_stage, .data$t_stage, .data$n_stage,.data$m_stage) %>% dplyr::filter(stage_system=="AJCC Stage", !grepl("Stage 1|Stage 2", .data$overall_stage)) %>% dplyr::group_by(.data$patient_id) %>% fill(everything(), .direction = "down") %>% fill(everything(), .direction = "up") %>% dplyr::slice_head(n=1)%>% dplyr::select(-stage_system) %>% unique %>% data.frame
  patient_metadata <- tempusr::prepare_metadata_patient(data_cohort) %>%
    dplyr::select(.data$patient_id, histology = .data$cancer_type, .data$sex, .data$race, .data$ethnicity, .data$age_pdx) %>%
    dplyr::mutate(age_categorical = case_when(age_pdx < 49 ~ "<49", age_pdx > 49 & age_pdx < 60 ~ "50-59", age_pdx > 59 & age_pdx < 70 ~ "60-69", age_pdx > 69 & age_pdx < 79 ~ "70-79", age_pdx > 79 ~ ">80", .default = "Unknown")) %>%
    dplyr::select(-age_pdx)


  demographics_list$Patient_Data <- patient_metadata
  #### Sample metadata
  sample_metadata <- tempusr::prepare_metadata_sample(data_cohort) %>%
    dplyr::select(patient_id, tumor_biospecimen_id = biospecimen_id, assay_dna, assay_rna, biopsy_procedure, biopsy_site, tissue_rollup, metastatic_status, purity = final_tumor_percentage_dna) %>%
    dplyr::mutate(purity_categorical = case_when(purity < 30 ~ "<30%", purity >= 30 & purity <= 70 ~ "30-70%", purity > 70 ~ ">70%", .default = "Unknown")) %>%
    dplyr::select(-purity)

  demographics_list$Sample_Data <- sample_metadata

  #### LOT

  txs <- paste0(treatments_of_interest, collapse = "|")
  dplyr::filter(data_cohort$onco_regimen, grepl(txs, agents)) %>%
    arrange(regimen_sequence) %>%
    group_by(patient_id) %>%
    dplyr::slice_head(n = 1) %>%
    ungroup()

  LOT_summary <- dplyr::filter(data_cohort$onco_regimen, grepl(txs, agents)) %>%
    arrange(regimen_sequence) %>%
    group_by(patient_id) %>%
    dplyr::slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(patient_id, agents, therapy_class_group, regimen_sequence, line_of_therapy_number)

  demographics_list$LOT <- LOT_summary

#### Metastatic site locations

   tx_regimen <- dplyr::filter(data_cohort$onco_regimen, grepl(txs, agents)) %>%
    arrange(regimen_sequence) %>%
    group_by(patient_id) %>%
    dplyr::slice_head(n = 1) %>%
       ungroup() %>% dplyr::select(patient_id,regimen_name=agents, Event_Date=start_date_indexed, Event_End_Date=end_date_indexed)


  met_locs <- data_cohort$onco_metastasis %>% dplyr::select(patient_id, Met_Date=metastasis_date_indexed, site) %>% group_by(patient_id,Met_Date) %>% dplyr::mutate(Combined_Sites=paste0(sort(site), collapse=",")) %>% dplyr::select(-site) %>% unique %>% arrange(patient_id, Met_Date) %>% group_by(patient_id) %>% dplyr::slice_head(n=1) %>% left_join(tx_regimen) %>% dplyr::mutate(Metastatic_Status_at_Tx= case_when( Met_Date <=Event_Date ~"Metastatic before treatment", Met_Date> Event_Date & Met_Date <=Event_End_Date ~ "Metastatic on treatment",Met_Date>Event_End_Date~"Metastatic after treatment end",Met_Date >Event_Date ~ "Metastatic after treatment start",.default="Unknown"), Has_Liver_Mets=grepl("Liver|liver", Combined_Sites), Has_Brain_Mets=grepl("Brain|brain", Combined_Sites)) %>% dplyr::select(patient_id, Combined_Sites,Metastatic_Status_at_Tx,Has_Liver_Mets, Has_Brain_Mets) %>% unique %>% ungroup


  demographics_list$Metastasis <- met_locs





  #### Biopsy timing relative to treatment



  ######## REALLY LONG BIT OF CODE
  metadata_patient <- prepare_metadata_patient(data_cohort)
  metadata_sample <- prepare_metadata_sample(data_cohort)


  event_timelines <- dplyr::select(metadata_patient, patient_id, Event_Date = pdx_date) %>% dplyr::mutate(Event = "Diagnosis", Assay = "Diagnosis")
  data_model <- tempusr:::calc_data_model(data_cohort, "regimen", "onco_regimen")

  if (data_model == "1.0") {
    regimen <- dplyr::select(data_cohort[[grep("regimen", names(data_cohort))]], patient_id, regimen_name, regimen_class, regimen_class_group, regimen_rank, Event_Date = regimen_start_date_indexed, Event_End_Date = regimen_end_date_indexed) %>% dplyr::mutate(Event = "Treatment", Assay = "Treatment")
  } else if (data_model == "2.0") {
    regimen <- dplyr::select(data_cohort[[grep("regimen", names(data_cohort))]], patient_id, regimen_name = agents, regimen_class = therapy_class, regimen_class_group = therapy_class_group, regimen_rank = regimen_sequence, Event_Date = start_date_indexed, Event_End_Date = end_date_indexed) %>% dplyr::mutate(Event = "Treatment", Assay = "Treatment")
  }

  regimen$Event_Date <- as.Date(regimen$Event_Date)
  regimen$Event_End_Date <- as.Date(regimen$Event_End_Date)




  biospecimen <- metadata_sample

  NGS <- rbind(dplyr::select(dplyr::filter(biospecimen, !is.na(assay_dna)), patient_id, tumor_biospecimen_id = sample_id_dna, Event_Date = biopsy_date, Assay = assay_dna), dplyr::select(dplyr::filter(biospecimen, !is.na(assay_rna)), patient_id, tumor_biospecimen_id = sample_id_rna, Event_Date = biopsy_date, Assay = assay_rna)) %>% dplyr::mutate(Event = "Tempus_NGS")



  final_timelines <- bind_rows(event_timelines, NGS, regimen) %>%
    dplyr::filter(!is.na(Event_Date)) %>%
    group_by(patient_id, Event_Date) %>%
    dplyr::mutate(Combined_Event = paste0(sort(unique(Event)), collapse = ","), Combined_Assay = paste0(sort(unique(setdiff(Assay, NA))), collapse = ",")) %>%
      ungroup()

  final_timelines$Combined_Assay <- ifelse(grepl(txs, final_timelines$regimen_name), paste0(final_timelines$Combined_Assay, ":",txs), final_timelines$Combined_Assay)

    final_timelines <- final_timelines %>%
    arrange(patient_id, Event_Date)




  final_timelines <- final_timelines %>% group_by(patient_id) %>%
    dplyr::mutate(Anchor_Date = ifelse(lubridate::is.Date(min(Event_Date[grepl(txs, Combined_Assay)])), min(Event_Date[grepl(txs, Combined_Assay)]), NA), Anchor_Date = as.Date(Anchor_Date)) %>%
    group_by(patient_id) %>%
      dplyr::mutate(Relative_Timing = case_when(Event_Date < Anchor_Date | Event_End_Date < Anchor_Date  ~ "Pre-Treatment" , Event_Date > Anchor_Date & Event_End_Date <=Anchor_Date ~ "On Treatment", Event_Date ==Anchor_Date ~ "On Treatment",Event_Date > Anchor_Date ~ "After Treatment",.default="Unknown")) %>% dplyr::select(-Event,-Assay,-tumor_biospecimen_id,-regimen_name,-regimen_rank,-regimen_class,-regimen_class_group) %>% unique %>% ungroup



  final_timelines <- final_timelines %>% group_by(patient_id) %>% dplyr::mutate(Event_Order=1:length(Combined_Event) ,Anchor_Index = ifelse(is.integer(unique(Event_Order[which(Anchor_Date == Event_Date)])), unique(Event_Order[which(Anchor_Date == Event_Date)])[1], NA)) %>% data.frame


  final_timelines$Interval_to_Anchor <- final_timelines$Event_Date - final_timelines$Anchor_Date
  final_timelines$Final_Order <- final_timelines$Event_Order - final_timelines$Anchor_Index


final_timelines <- dplyr::select(final_timelines, patient_id,  Combined_Event, Combined_Assay, Relative_Timing, Interval_to_Anchor,Relative_Order=Final_Order) %>% unique %>% ungroup()


  demographics_list$Biopsy_Timing <- final_timelines %>% dplyr::filter(grepl("Tempus_NGS", Combined_Event))


#####




  #### Smoking status
  smoking_summary <- data_cohort$onco_smoking_status %>%
    dplyr::select(patient_id, smoking_status_date_indexed, standardized_smoking_status) %>%
    group_by(patient_id) %>%
    arrange(smoking_status_date_indexed) %>%
    dplyr::slice_head(n = 1) %>% dplyr::select(-smoking_status_date_indexed) %>%
    ungroup()


  demographics_list$Smoking <- smoking_summary


#### IHC data
  IHC_summary <- dplyr::select(input_td$onco_result_ihc_passing, patient_id, tumor_biospecimen_id=biospecimen_id, gene_symbol, fraction_positive=tumor_cell_staining_percentage,biopsy_site=collection_anatomical_site) %>% dplyr::mutate(fraction_positive_categorical=case_when(fraction_positive<1 ~ "<1%", fraction_positive>=1 & fraction_positive<49 ~"1-49%", fraction_positive >49~">=50%", .default="Unknown"), gene_symbol =gsub("CD274","PD-L1",gene_symbol)) %>% dplyr::select(-fraction_positive) %>% group_by(patient_id,gene_symbol) %>% dplyr::mutate(Num_Samples_Tested=length(tumor_biospecimen_id)) %>% dplyr::select(-tumor_biospecimen_id) %>% unique %>% dplyr::slice_head(n=1) %>% ungroup



  demographics_list$Pathos_IHC <- IHC_summary


  third_party_IHC <- dplyr::select(input_td$onco_reported_third_party_overview_marker,patient_id, gene_symbol, marker_name, analysis_method, sample_type, marker_result_categorical, numerical_result,numerical_result_type, sample_site, sample_tumor_type) %>% dplyr::filter(grepl("IHC|ISH", analysis_method))


  demographics_list$External_IHC <- third_party_IHC


  return(demographics_list)
  }
