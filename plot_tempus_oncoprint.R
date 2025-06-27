#' Plot oncoprint of patients in Tempus data
#'
#' @param mmf Data frame of the Tempus mutations, either Master Molecular File or Master Molecular File Filtered (data.frame)
#' @param alteration_func A named list object that assigns color and sizes to different mutation types
#' @param patient_identifier The column in mmf where patient identifiers are located
#' @param gene_identifier The column in mmf where gene identifiers are located
#' @param title Plot title (Optional string)
#' @param filter_germline This is a boolean for whether germline mutations should be removed from consideration. If you've removed the somatic_germline column from the mmf, this should be changed to FALSE
#' @param copy_number_threshold The minimum copy number required to call an increase in copy number an amplification. Anything below this value will be classified as a "weak amplification"
#' @param patient_order The ordering of patients during final plotting. This is helpful if you plan to add additional annotations later.
#' @param patient_universe The total set of patients in the universe (only necessary if using a filtered version of the mmf file as input where some patients may not have alterations.

#'
#'
#'
#' @return Complex Heatmap Object featuring your input mutations
#' @note Any filtering or other grouping you wish to perform need to be done before passing the mmf table.
#'
#'
#' @export
#' @examples
#' \dontrun{
#' pathos_oncoprint(mmf,
#'   alteration_func = NULL,
#'   gene_identifier = "gene_canonical_name",
#'   title = NULL,
#' filter_germline=TRUE,
#' copy_number_threshold=8,
#' patient_order=NULL,
#' patient_universe=NULL
#'
#' )
#' }
#'
plot_tempus_oncoprint <- function(mmf, alteration_func = NULL, patient_identifier = "analysis_id", gene_identifier = "gene_canonical_name", title = NULL, filter_germline = TRUE, copy_number_threshold = 8, patient_order = NULL,patient_universe=NULL) {

    mmf <- as.data.frame(mmf)
  colnames(mmf)[grep(paste0("^", gene_identifier, "$"), colnames(mmf))] <- "Gene"

  if (patient_identifier %in% c("patient_id", "analysis_id")) {
    mmf_SNV <- dplyr::filter(mmf, .data$variant_type_code == "SHRTVRNT", .data$functional_impact != "B") %>%
      dplyr::select(.data$patient_id, .data$analysis_id, .data$Gene, .data$variant_type_code, .data$result, .data$functional_impact, .data$mutation_effect, .data$variant_allele_freq, .data$somatic_germline)
  } else {
    mmf_SNV <- dplyr::filter(mmf, .data$variant_type_code == "SHRTVRNT", .data$functional_impact != "B") %>%
      dplyr::select(.data$patient_id, .data$analysis_id, !!as.name(patient_identifier), .data$Gene, .data$variant_type_code, .data$result, .data$functional_impact, .data$mutation_effect, .data$variant_allele_freq, .data$somatic_germline)
  }
  if (filter_germline) {
    mmf_SNV <- dplyr::filter(mmf_SNV, .data$somatic_germline == "S")
  }

  if (patient_identifier %in% c("patient_id", "analysis_id")) {
    mmf_CNV_ref <- dplyr::filter(mmf, .data$variant_type_code == "CNALTER") %>%
      dplyr::select(.data$patient_id, .data$analysis_id, .data$Gene, .data$variant_type_code, .data$result, .data$functional_impact, .data$somatic_germline, .data$copy_number)
  } else {
    mmf_CNV_ref <- dplyr::filter(mmf, .data$variant_type_code == "CNALTER") %>%
      dplyr::select(.data$patient_id, .data$analysis_id, !!as.name(patient_identifier), .data$Gene, .data$variant_type_code, .data$result, .data$functional_impact, .data$result, .data$somatic_germline, .data$copy_number)
  }

  mmf_CNV_calc <- tryCatch(
    {
      tempusr::calc_cnv(mmf)
    },
    error = function(e) {
      "Could not run TempusR::Calc_CNV"
    }
  )


  if (!is.character(mmf_CNV_calc)) {
    colnames(mmf_CNV_calc)[2] <- "Gene"
    mmf_CNV_ref <- dplyr::left_join(dplyr::select(mmf_CNV_ref, -(.data$copy_number)), mmf_CNV_calc, by = c("analysis_id", "Gene"))
  } else if (!is.character(mmf_CNV_calc)) {
    ## print(mmf_CNV_calc)
    mmf_CNV_ref <- mmf_CNV_ref %>%
      dplyr::group_by(.data$patient_id, .data$analysis_id, .data$Gene) %>%
      dplyr::mutate(copy_number = mean(.data$copy_number)) %>%
      unique() %>%
      dplyr::ungroup()
  }

  mmf_CNV_ref <- mmf_CNV_ref %>% dplyr::mutate(Mut = dplyr::case_when(.data$copy_number == 0 ~ "HOMDEL", .data$copy_number == 1 ~ "SHALLOWDEL", .data$copy_number < copy_number_threshold & .data$copy_number > 2 ~ "WEAKAMP",.data$copy_number >= copy_number_threshold ~ "AMP"))

  mmf_SNV <- mmf_SNV %>% dplyr::mutate(Mut = dplyr::case_when(grepl("stop", .data$mutation_effect) ~ "TRUNC", !grepl("stop", .data$mutation_effect) ~ "MISSENSE"))


 if (patient_identifier %in% c("patient_id", "analysis_id")) {
  mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut))} else{     mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, !!as.name(patient_identifier),.data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id,!!as.name(patient_identifier), .data$Gene, .data$Mut))}

  colnames(mmf_int)[grep(paste0("^", patient_identifier, "$"), colnames(mmf_int))] <- "ID"
  ## str(mmf_int)

  mmf_out <- mmf_int %>%
    dplyr::group_by(.data$ID, .data$Gene) %>%
    dplyr::mutate(Combo_Mut = paste0(sort(unique(.data$Mut)), collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$ID, .data$Gene, .data$Combo_Mut) %>%
    unique() %>%
    data.frame()


  mmf_final <- tidyr::pivot_wider(mmf_out, names_from = .data$ID, id_cols = .data$Gene, values_from = .data$Combo_Mut) %>% data.frame()

  rownames(mmf_final) <- mmf_final[, 1]
    mmf_final[, 1] <- NULL
    saveRDS(mmf_final,"mmf_final.rds")
  mmf_final <- as.matrix(mmf_final)
  ## stop()
  colnames(mmf_final) <- gsub("^X","", gsub("\\.","-",colnames(mmf_final)))
  
    if(!is.null(patient_universe)){
        warning("A patient universe was provided that may not overlap with the patients in the provided mmf file. This could lead to underestimation of the mutation frequency in the patient cohort") 
        missing_patients <- setdiff(colnames(mmf_final),x= patient_universe)
        missing_mat <- matrix(NA, nrow=nrow(mmf_final), ncol=length(missing_patients))
        colnames(missing_mat) <- missing_patients
        rownames(missing_mat) <- rownames(mmf_final)
        mmf_final <- cbind(mmf_final, missing_mat)}
        
  if (is.null(title)) {
    title <- "Oncoprint"
  } else {
    title <- title
  }

  if (is.null(alteration_func)) {
    col <- c("HOMDEL" = "blue", "SHALLOWDEL" = "lightblue", "AMP" = "red", "WEAKAMP" = "lightpink", "MISSENSE" = "darkgreen", "TRUNC" = "black")
    alteration_func <- list(
      background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
      HOMDEL = ComplexHeatmap::alter_graphic("rect", fill = col[["HOMDEL"]]),
      SHALLOWDEL = ComplexHeatmap::alter_graphic("rect", fill = col[["SHALLOWDEL"]]),
      AMP = ComplexHeatmap::alter_graphic("rect", fill = col[["AMP"]]),
      WEAKAMP = ComplexHeatmap::alter_graphic("rect", fill = col[["WEAKAMP"]]),
      MISSENSE = ComplexHeatmap::alter_graphic("rect", height = 0.6, fill = col[["MISSENSE"]]),
      TRUNC = ComplexHeatmap::alter_graphic("rect", height = 0.4, fill = col[["TRUNC"]])
    )
    ht_legend_param <- list(title = "Alterations", at = names(col), labels = c("Deep Deletion", "Shallow Deletion", paste0("Amplification (>=", copy_number_threshold, ")"), paste0("Weak Amplification (3-", copy_number_threshold - 1, ")"), "Missense Mutation", "Truncating Mutation"))

    if(is.null(patient_order)){
        ht <- ComplexHeatmap::oncoPrint(mmf_final, alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = FALSE, remove_empty_rows = TRUE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE)}
    else if( !is.null(patient_order)){


            ht <- ComplexHeatmap::oncoPrint(mmf_final[,patient_order], alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = FALSE, remove_empty_rows = TRUE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE, column_order=patient_order)}
    return(ht)
  } else if (all(unique(unlist(strsplit(mmf_out$Combo_Mut, ";"))) %in% names(alteration_func))) {
    print("This hasn't been implemented yet")
    ## return(ht)
  } else {
    print("There are mutation types in the data that do not have a plot mapping. Please submit a custom plot mapping")
  }
}
