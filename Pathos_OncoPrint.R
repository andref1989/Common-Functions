#' Plot treatment journeys of patients in Tempus data
#'
#' @param mmf Data frame of the Tempus mutations, either Master Molecular File or Master Molecular File Filtered (data.frame)
#' @param alteration_func A named list object that assigns color and sizes to different mutation types
#' @param patient_identifier The column in mmf where patient identifiers are located
#' @param gene_identifier The column in mmf where gene identifiers are located
#' @param title Plot title (Optional string)

#'
#'
#'
#' @return Complex Heatmap Object featuring your input mutations
#' @note Any filtering or other grouping you wish to perform need to be done before passing the regimen table. Facetting has not been tested and likely will not work. Trying checks again and again and again... seriously will this ever work?? Probably not.
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' pathos_oncoprint(mmf,
#'   alteration_func = NULL,
#'   gene_identifier = "gene_canonical_name",
#'   title = NULL
#' )
#' }
#'
pathos_oncoprint <- function(mmf, alteration_func = NULL, patient_identifier = "analysis_id", gene_identifier = "gene_canonical_name", title = NULL, filter_germline = TRUE, copy_number_threshold = 8, patient_order = NULL) {
  colnames(mmf)[grep(paste0("^", gene_identifier, "$"), colnames(mmf))] <- "Gene"

  if (patient_identifier %in% c("patient_id", "analysis_id")) {
    mmf_SNV <- dplyr::filter(mmf, variant_type_code == "SHRTVRNT", functional_impact != "B") %>%
      dplyr::select(patient_id, analysis_id, Gene, variant_type_code, result, functional_impact, mutation_effect, variant_allele_freq, somatic_germline)
  } else {
    mmf_SNV <- dplyr::filter(mmf, variant_type_code == "SHRTVRNT", functional_impact != "B") %>%
      dplyr::select(patient_id, analysis_id, !!as.name(patient_identifier), Gene, variant_type_code, result, functional_impact, mutation_effect, variant_allele_freq, somatic_germline)
  }
  if (filter_germline) {
    mmf_SNV <- dplyr::filter(mmf_SNV, somatic_germline == "S")
  }

  if (patient_identifier %in% c("patient_id", "analysis_id")) {
    mmf_CNV_ref <- dplyr::filter(mmf, variant_type_code == "CNALTER") %>%
      dplyr::select(patient_id, analysis_id, Gene, variant_type_code, result, functional_impact, result, somatic_germline, copy_number)
  } else {
    mmf_CNV_ref <- dplyr::filter(mmf, variant_type_code == "CNALTER") %>%
      dplyr::select(patient_id, analysis_id, !!as.name(patient_identifier), Gene, variant_type_code, result, functional_impact, result, somatic_germline, copy_number)
  }

  mmf_CNV_calc <- tryCatch(
    {
      tempusr::calc_cnv(mmf)
    },
    error = function(e) {
      "Could not run TempusR::Calc_CNV"
    }
  )


  if (class(mmf_CNV_calc) != "character") {
    colnames(mmf_CNV_calc)[2] <- "Gene"
    mmf_CNV_ref <- left_join(dplyr::select(mmf_CNV_ref, -(copy_number)), mmf_CNV_calc, by = c("analysis_id", "Gene"))
  } else if (class(mmf_CNV_calc) == "character") {
    ## print(mmf_CNV_calc)
    mmf_CNV_ref <- mmf_CNV_ref %>%
      group_by(patient_id, analysis_id, Gene) %>%
      mutate(copy_number = mean(copy_number)) %>%
      unique() %>%
      ungroup()
  }

  mmf_CNV_ref <- mmf_CNV_ref %>% mutate(Mut = case_when(copy_number == 0 ~ "HOMDEL", copy_number == 1 ~ "SHALLOWDEL", copy_number < copy_number_threshold & copy_number > 2 ~ "WEAKAMP", copy_number >= copy_number_threshold ~ "AMP"))

  mmf_SNV <- mmf_SNV %>% mutate(Mut = case_when(grepl("stop", mutation_effect) ~ "TRUNC", !grepl("stop", mutation_effect) ~ "MISSENSE"))


 if (patient_identifier %in% c("patient_id", "analysis_id")) {
  mmf_int <- rbind(dplyr::select(mmf_CNV_ref, patient_id, analysis_id, Gene, Mut), dplyr::select(mmf_SNV, patient_id, analysis_id, Gene, Mut))} else{     mmf_int <- rbind(dplyr::select(mmf_CNV_ref, patient_id, analysis_id, !!as.name(patient_identifier),Gene, Mut), dplyr::select(mmf_SNV, patient_id, analysis_id,!!as.name(patient_identifier), Gene, Mut))}

  colnames(mmf_int)[grep(paste0("^", patient_identifier, "$"), colnames(mmf_int))] <- "ID"
  ## str(mmf_int)

  mmf_out <- mmf_int %>%
    group_by(ID, Gene) %>%
    mutate(Combo_Mut = paste0(sort(unique(Mut)), collapse = ";")) %>%
    ungroup() %>%
    dplyr::select(ID, Gene, Combo_Mut) %>%
    unique() %>%
    data.frame()


  mmf_final <- tidyr::pivot_wider(mmf_out, names_from = ID, id_cols = Gene, values_from = Combo_Mut) %>% data.frame()

  rownames(mmf_final) <- mmf_final[, 1]
  mmf_final[, 1] <- NULL
  mmf_final <- as.matrix(mmf_final)
  ## stop()
  colnames(mmf_final) <- gsub("^X","", gsub("\\.","-",colnames(mmf_final)))


  if (is.null(title)) {
    title <- "Oncoprint"
  } else {
    title <- title
  }

  if (is.null(alteration_func)) {
    col <- c("HOMDEL" = "blue", "SHALLOWDEL" = "lightblue", "AMP" = "red", "WEAKAMP" = "lightpink", "MISSENSE" = "darkgreen", "TRUNC" = "black")
    alteration_func <- list(
      background = alter_graphic("rect", fill = "#CCCCCC"),
      HOMDEL = alter_graphic("rect", fill = col[["HOMDEL"]]),
      SHALLOWDEL = alter_graphic("rect", fill = col[["SHALLOWDEL"]]),
      AMP = alter_graphic("rect", fill = col[["AMP"]]),
      WEAKAMP = alter_graphic("rect", fill = col[["WEAKAMP"]]),
      MISSENSE = alter_graphic("rect", height = 0.6, fill = col[["MISSENSE"]]),
      TRUNC = alter_graphic("rect", height = 0.4, fill = col[["TRUNC"]])
    )
    ht_legend_param <- list(title = "Alterations", at = names(col), labels = c("Deep Deletion", "Shallow Deletion", paste0("Amplification (>=", copy_number_threshold, ")"), paste0("Weak Amplification (3-", copy_number_threshold - 1, ")"), "Missense Mutation", "Truncating Mutation"))

    if(is.null(patient_order)){
        ht <- ComplexHeatmap::oncoPrint(mmf_final, alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = TRUE, remove_empty_rows = TRUE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE)}
    else if( !is.null(patient_order)){
            print(patient_order)

            ht <- ComplexHeatmap::oncoPrint(mmf_final[,patient_order], alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = TRUE, remove_empty_rows = TRUE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE, column_order=patient_order)}
    return(ht)
  } else if (all(unique(unlist(strsplit(mmf_out$Combo_Mut, ";"))) %in% names(alteration_func))) {
    print("This hasn't been implemented yet")
    ## return(ht)
  } else {
    print("There are mutation types in the data that do not have a plot mapping. Please submit a custom plot mapping")
  }
}
