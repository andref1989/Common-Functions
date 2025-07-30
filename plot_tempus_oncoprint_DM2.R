#' Plot oncoprint of patients in Tempus data
#'
#' @param data_cohort Path to patient cohort of interest or the named list object output by load_tempus_data. The path option is slower and not advised if multiple oncoprints will be generated sequentially.
#' @param genes_of_interest Vector of genes the user wants plotted (MANDATORY if more than 50 genes detected in mmf)
#' @param patient_identifier The column in mmf where patient identifiers are located
#' @param gene_identifier The column in mmf where gene identifiers are located
#' @param title Plot title (Optional string)
#' @param filter_germline This is a boolean for whether germline mutations should be removed from consideration. If you've removed the somatic_germline column from the mmf, this should be changed to FALSE
#' @param alteration_func A named list object that assigns color and sizes to different mutation types (Optional)
#' @param copy_number_threshold The minimum copy number required to call an increase in copy number an amplification. Anything below this value will be classified as a "weak amplification" (integer)
#' @param patient_order The ordering of patients during final plotting. This is helpful if you plan to add additional annotations later. (vector)
#' @param patient_universe The total set of patients in the universe (only necessary if using a filtered version of the mmf file as input where some patients may not have alterations. (Optional vector)
#' @param assay_blacklist The specific assays you want to exclude from the oncoprint (Optional string/vector) 

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
plot_tempus_oncoprint <- function(data_cohort, genes_of_interest=NULL,patient_identifier = "patient_id", gene_identifier = c("gene_canonical_name","gene_symbol"), title = NULL, filter_germline = TRUE, alteration_func = NULL, copy_number_threshold = 8, patient_order = NULL,patient_universe=NULL,assay_blacklist=NULL) {

##########################################################################################
############  Importing the tempus data and identifying the data model
##########################################################################################
    if(is.character(data_cohort)){

        input_td <- load_tempus_data(data_cohort,collection="molecular")
    } else if(is.list(data_cohort)){

        input_td <- data_cohort}

############  Checking the data model ##################

    data_model_check <- calc_data_model(input_td,c("g_molecular_master_file"),c("onco_result_snv_indel_passing","onco_result_cnv_gene"))
    stopifnot(data_model_check %in% c("1.0","2.0"))



############  Define DM specific functions ##################
    prep_oncoprint_dm1 <- function(input_td){

        mmf <- input_td[grep("molecular_master_file", names(input_td))]
        gene_identifier <- unique(unlist(lapply(mmf, function(x) intersect(gene_identifier,colnames(x)))))

        if(is.null(genes_of_interest) && length(unique(mmf[[1]][,gene_identifier])) > 50){
            stop("There are more than 50 unique genes in this mmf file, did you forget to set any genes of interest?")
            } else if (is.null(genes_of_interest) && length(unique(mmf[[1]][,gene_identifier])) < 50){
                      common_cols <- Reduce("intersect",lapply(mmf, function(x) colnames(x)))
                      mmf <- do.call("rbind", lapply(mmf, function(x) x[,common_cols]))
                      mmf <- as.data.frame(mmf)} else if (!is.null(genes_of_interest) && length(genes_of_interest) <50){

                                                   common_cols <- Reduce("intersect",lapply(mmf, function(x) colnames(x)))
                      mmf <- do.call("rbind", lapply(mmf, function(x) x[which(x[,gene_identifier] %in% genes_of_interest),common_cols]))
                                                   mmf <- as.data.frame(mmf)} else if (!is.null(genes_of_interest && length(genes_of_interest) > 50)) { stop("There are more than 50 specified genes of interest. Consider reducing that number and replotting")
}


        colnames(mmf)[grep(paste0("^", gene_identifier, "$"), colnames(mmf))] <- "Gene"
        if(!is.null(assay_blacklist)){
        mmf <- dplyr::filter(mmf, !grepl(paste0(assay_blacklist,collapse="|"), assay))}


        if (patient_identifier %in% c("patient_id", "analysis_id")) {
      mmf_SNV <- dplyr::filter(mmf, .data$variant_type_code == "SHRTVRNT", .data$functional_impact != "B",!grepl("stream|intron|UTR",mutation_effect)) %>%
      dplyr::select(.data$patient_id, .data$analysis_id, .data$Gene, .data$variant_type_code, .data$result, .data$functional_impact, .data$mutation_effect, .data$variant_allele_freq, .data$somatic_germline)
  } else {
    mmf_SNV <- dplyr::filter(mmf, .data$variant_type_code == "SHRTVRNT", .data$functional_impact != "B",!grepl("stream|intron|UTR",mutation_effect)) %>%
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

    mmf_CNV_ref <- mmf_CNV_ref %>%
      dplyr::group_by(.data$patient_id, .data$analysis_id, .data$Gene) %>%
      dplyr::mutate(copy_number = mean(.data$copy_number)) %>%
      unique() %>%
      dplyr::ungroup()
  }

  mmf_CNV_ref <- mmf_CNV_ref %>% dplyr::mutate(Mut = dplyr::case_when(.data$copy_number == 0 ~ "HOMDEL", .data$copy_number == 1 ~ "SHALLOWDEL", .data$copy_number < copy_number_threshold & .data$copy_number > 2 ~ "WEAKAMP",.data$copy_number >= copy_number_threshold ~ "AMP"))

  mmf_SNV <- mmf_SNV %>% dplyr::mutate(Mut = dplyr::case_when(grepl("stop", .data$mutation_effect) ~ "TRUNC", grepl("missense", .data$mutation_effect) ~ "MISSENSE"),grepl("splice", .data$mutation_effect) ~ "SPLICE",.default="OTHER")

        ## saveRDS(mmf_CNV_ref,"~/CNV_Ref_DM1.rds")
        ## saveRDS(mmf_SNV,"~/SNV_Ref_DM1.rds")

 if (patient_identifier %in% c("patient_id", "analysis_id")) {
  mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut))} else{     mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, !!as.name(patient_identifier),.data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id,!!as.name(patient_identifier), .data$Gene, .data$Mut))}

  colnames(mmf_int)[grep(paste0("^", patient_identifier, "$"), colnames(mmf_int))] <- "ID"

        return(mmf_int)
        }


    prep_oncoprint_dm2 <- function(input_td){
        mmf <- input_td[grep("onco_result_cnv_gene|onco_result_snv", names(input_td))]

        gene_identifier <- unique(unlist(lapply(mmf, function(x) intersect(gene_identifier,colnames(x)))))
        mmf_SNV <- mmf[grep("snv", names(mmf))]
        mmf_CNV <- mmf[grep("cnv_gene", names(mmf))]


                    if(is.null(genes_of_interest) && (length(unique(mmf_SNV[[1]][,gene_identifier])) > 50||length(unique(mmf_CNV[[1]][,gene_identifier])) > 50)){
            stop("There are more than 50 unique genes in this mmf file, did you forget to set any genes of interest?")
                    } else if (is.null(genes_of_interest) && (length(unique(mmf[[1]][,gene_identifier])) < 50 && length(unique(mmf_CNV[[1]][,gene_identifier])) < 50)){
                common_cols1 <- Reduce("intersect",lapply(mmf_SNV, function(x) colnames(x)))
                common_cols2 <- Reduce("intersect",lapply(mmf_CNV, function(x) colnames(x)))
                mmf_SNV <- do.call("rbind", lapply(mmf_SNV, function(x) x[,common_cols1]))
                mmf_CNV <- do.call("rbind", lapply(mmf_CNV, function(x) x[,common_cols2]))
                mmf_SNV <- as.data.frame(mmf_SNV)
                mmf_CNV <- as.data.frame(mmf_CNV)

                    } else if (!is.null(genes_of_interest) && length(genes_of_interest) <50){


                common_cols1 <- Reduce("intersect",lapply(mmf_SNV, function(x) colnames(x)))
                common_cols2 <- Reduce("intersect",lapply(mmf_CNV, function(x) colnames(x)))
                mmf_SNV <- do.call("rbind", lapply(mmf_SNV, function(x) x[,common_cols1]))
                mmf_CNV <- do.call("rbind", lapply(mmf_CNV, function(x) x[,common_cols2]))
                mmf_SNV <- as.data.frame(mmf_SNV)
                mmf_CNV <- as.data.frame(mmf_CNV)



                        mmf_SNV <- dplyr::filter(mmf_SNV, !!as.name(gene_identifier) %in% genes_of_interest)
                        mmf_CNV <- dplyr::filter(mmf_CNV, !!as.name(gene_identifier) %in% genes_of_interest)

            } else if (!is.null(genes_of_interest && length(genes_of_interest) > 50)) { stop("There are more than 50 specified genes of interest. Consider reducing that number and replotting")
}


            colnames(mmf_SNV)[grep(paste0("^", gene_identifier, "$"), colnames(mmf_SNV))] <- "Gene"

            colnames(mmf_CNV)[grep(paste0("^", gene_identifier, "$"), colnames(mmf_CNV))] <- "Gene"

            if(!is.null(assay_blacklist)){
            mmf_CNV <- dplyr::filter(mmf_CNV, !grepl(paste0(assay_blacklist,collapse="|"), assay))
            mmf_SNV <- dplyr::filter(mmf_SNV, !grepl(paste0(assay_blacklist,collapse="|"), assay))}


        if (patient_identifier %in% c("patient_id", "analysis_id")) {
      mmf_SNV <- dplyr::filter(mmf_SNV, .data$variant_classification != "B",!grepl("stream|intron|UTR",variant_molecular_consequence)) %>%
          dplyr::select(.data$patient_id, analysis_id=.data$tumor_biospecimen_id, .data$Gene, functional_impact=.data$variant_classification,mutation_effect=.data$variant_molecular_consequence, variant_allele_freq= .data$tumor_variant_allele_frequency, somatic_germline=dplyr::any_of(c("variant_original_source","variant_origin")))
  } else {
    mmf_SNV <- dplyr::filter(mmf_SNV, .data$variant_classification != "B",!grepl("stream|intron|UTR",variant_molecular_consequence)) %>%
      dplyr::select(.data$patient_id, analysis_id=.data$tumor_biospecimen_id, !!as.name(patient_identifier), .data$Gene, functional_impact=.data$variant_classification, mutation_effect=.data$variant_molecular_consequence, variant_allele_freq=.data$tumor_variant_allele_frequency, somatic_germline=dplyr::any_of(c("variant_original_source","variant_origin")))
  }
  if (filter_germline) {
    mmf_SNV <- dplyr::filter(mmf_SNV, .data$somatic_germline %in% c("S","somatic"))
  }

      if (patient_identifier %in% c("patient_id", "analysis_id")) {
    mmf_CNV_ref <- dplyr::select(mmf_CNV,.data$patient_id, analysis_id=.data$tumor_biospecimen_id, .data$Gene, .data$copy_number)
  } else {
    mmf_CNV_ref <- dplyr::select(mmf_CNV,.data$patient_id, analysis_id=.data$tumor_biospecimen_id, !!as.name(patient_identifier), .data$Gene, .data$copy_number)
  }

  mmf_CNV_calc <- tryCatch(
    {
      tempusr::calc_cnv(mmf_CNV)
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

            mmf_SNV <- mmf_SNV %>% dplyr::mutate(Mut = dplyr::case_when(grepl("stop", .data$mutation_effect) ~ "TRUNC", grepl("missense", .data$mutation_effect) ~ "MISSENSE",grepl("splice", .data$mutation_effect) ~ "SPLICE",.default="OTHER"))
 if (patient_identifier %in% c("patient_id", "analysis_id")) {
            mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id, .data$Gene, .data$Mut))} else{     mmf_int <- rbind(dplyr::select(mmf_CNV_ref, .data$patient_id, .data$analysis_id, !!as.name(patient_identifier),.data$Gene, .data$Mut), dplyr::select(mmf_SNV, .data$patient_id, .data$analysis_id,!!as.name(patient_identifier), .data$Gene, .data$Mut))}

  colnames(mmf_int)[grep(paste0("^", patient_identifier, "$"), colnames(mmf_int))] <- "ID"

        return(mmf_int)}

    
##########################################################################################
#### Calling formatting functions for oncoprint.
##########################################################################################

    if(data_model_check == "1.0"){

        mmf_int <- prep_oncoprint_dm1(input_td)

    }
    else if (data_model_check =="2.0")  {
        mmf_int <- prep_oncoprint_dm2(input_td)}




##########################################################################################
####### Finalizing formatting of data frames
##########################################################################################
 
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

  mmf_final <- as.matrix(mmf_final)
  colnames(mmf_final) <- gsub("^X","", gsub("\\.","-",colnames(mmf_final)))



##########################################################################################
######### Finalizing patient universe and plotting
##########################################################################################
        if(!is.null(patient_universe)){
        warning("A patient universe was provided that may not overlap with the patients in the provided mmf file. This could lead to underestimation of the mutation frequency in the patient cohort") 
        missing_patients <- setdiff(colnames(mmf_final),x= patient_universe)
        missing_mat <- matrix("", nrow=nrow(mmf_final), ncol=length(missing_patients))
        colnames(missing_mat) <- missing_patients
        rownames(missing_mat) <- rownames(mmf_final)
        mmf_final <- cbind(mmf_final, missing_mat)
        } else{

            mmf_locs <- input_td[grep(paste0(c("g_molecular_master_file","g_molecular_master_file_filtered","onco_result_cnv_chr_segment","onco_result_cnv_gene","onco_result_snv_indel_passing_filtered"),collapse="|"),names(input_td))]

            patient_universe <- unique(unlist(lapply(mmf_locs, function(x) x[,patient_identifier])))
            missing_patients <- setdiff(colnames(mmf_final),x= patient_universe)
            missing_mat <- matrix("", nrow=nrow(mmf_final), ncol=length(missing_patients))
            colnames(missing_mat) <- missing_patients
            rownames(missing_mat) <- rownames(mmf_final)
            mmf_final <- cbind(mmf_final, missing_mat)}




    if (is.null(title)) {
    title <- "Oncoprint"
  } else {
    title <- title
  }

  if (is.null(alteration_func)) {
    col <- c("HOMDEL" = "blue", "SHALLOWDEL" = "lightblue", "AMP" = "red", "WEAKAMP" = "lightpink", "MISSENSE" = "darkgreen", "TRUNC" = "black","SPLICE"="lightgreen","OTHER"="gray50")

    alteration_func <- list(
      background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
      HOMDEL = ComplexHeatmap::alter_graphic("rect", fill = col[["HOMDEL"]]),
      SHALLOWDEL = ComplexHeatmap::alter_graphic("rect", fill = col[["SHALLOWDEL"]]),
      AMP = ComplexHeatmap::alter_graphic("rect", fill = col[["AMP"]]),
      WEAKAMP = ComplexHeatmap::alter_graphic("rect", fill = col[["WEAKAMP"]]),
      MISSENSE = ComplexHeatmap::alter_graphic("rect", height = 0.6, fill = col[["MISSENSE"]]),
      TRUNC = ComplexHeatmap::alter_graphic("rect", height = 0.4, fill = col[["TRUNC"]]),
      SPLICE = ComplexHeatmap::alter_graphic("rect", height = 0.3,fill = col[["SPLICE"]]),
      OTHER = alter_graphic("rect", height = 0.2,width=0.2, fill = col[["OTHER"]])

    )
          mmf_final[is.na(mmf_final)] <- ""
    ht_legend_param <- list(title = "Alterations", at = names(col), labels = c("Deep Deletion", "Shallow Deletion", paste0("Amplification (>=", copy_number_threshold, ")"), paste0("Weak Amplification (3-", copy_number_threshold - 1, ")"), "Missense Mutation", "Truncating Mutation","Splice Region","Other"))

    if(is.null(patient_order)){
        ht <- suppressMessages(ComplexHeatmap::oncoPrint(mmf_final, alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = FALSE, remove_empty_rows = FALSE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE))}
    else if( !is.null(patient_order)){


            ht <- suppressMessages(ComplexHeatmap::oncoPrint(mmf_final[,patient_order], alter_fun = alteration_func, col = col, column_title = title, remove_empty_columns = FALSE, remove_empty_rows = FALSE, heatmap_legend_param = ht_legend_param, alter_fun_is_vectorized = TRUE, column_order=patient_order))}
    return(ht)
  } else if (all(unique(unlist(strsplit(mmf_out$Combo_Mut, ";"))) %in% names(alteration_func))) {
    print("This hasn't been implemented yet")
    ## return(ht)
  } else {
    print("There are mutation types in the data that do not have a plot mapping. Please submit a custom plot mapping")
  }
}
