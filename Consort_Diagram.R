library(consort)
consort_df <- clinical_data$onco_patient |>
  select(patient_id) |>
  mutate(metastatic_date_filter = ifelse(patient_id %in% met_pts$patient_id, 1, NA),
         mcrpc_ever_filter = ifelse(patient_id %in% mcrpc_pts$patient_id, 1, NA),
         mcrpc_with_dates_filter = ifelse(patient_id %in% mcrpc_with_dates_pts$patient_id, 1, NA),
         p2a_pts_filter = ifelse(patient_id %in% p2a_pts$patient_id, 1, NA),
         p2a_exclusion = case_when(!patient_id %in% mcrpc_with_dates_pts$patient_id ~ NA,
                                   patient_id %in% {phase2a_patients_recap |> filter(is.na(arpi_flag))}$patient_id ~ "No ARPI w/ Start Date",
                                   patient_id %in% {phase2a_patients_recap |> filter(is.na(arpi_post_mcrpc_flag))}$patient_id ~ "First ARPI NOT post-mCRPC",
                                   patient_id %in% {phase2a_patients_recap |> filter(is.na(lot_flag))}$patient_id ~ "No LoT Data",
                                   !patient_id %in% p2a_pts$patient_id ~ "No Derivable Anchor Date"
                                   ),
         p2a_exclusion = factor(p2a_exclusion,
                                levels = c("No ARPI w/ Start Date",
                                           "First ARPI NOT post-mCRPC",
                                           "No LoT Data",
                                           "No Derivable Anchor Date"
                                           )),
         biomarker_filter = ifelse(patient_id %in% biomarker_pts$patient_id, 1, NA)
         )

consort_diagram <- consort_plot(
      data = consort_df,
      orders = c(
        patient_id = "Prostate RNAseq",
        metastatic_date_filter = "Has Metastatic Date",
        mcrpc_ever_filter = "Ever mCRPC",
        mcrpc_with_dates_filter = "Has mCRPC Date",
        p2a_exclusion = "Excluded",
        p2a_pts_filter = "Phase 2a Like",
        biomarker_filter = "Biomarker Evaluable"
      ),
      side_box = c("p2a_exclusion")
    )

plot(consort_diagram)
