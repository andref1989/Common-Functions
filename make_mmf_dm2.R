make_mmf_dm2 <- function(input_td,patient_identifier="patient_id",gene_identifier="gene_symbol",assay_blacklist="xF",filter_germline=T, copy_number_threshold=8){
        mmf <- input_td[grep("onco_result_cnv_gene|onco_result_snv", names(input_td))]

        gene_identifier <- unique(unlist(lapply(mmf, function(x) intersect(gene_identifier,colnames(x)))))
        mmf_SNV <- mmf[grep("snv", names(mmf))]
        mmf_CNV <- mmf[grep("cnv_gene", names(mmf))]

        common_cols1 <- Reduce("intersect",lapply(mmf_SNV, function(x) colnames(x)))
        common_cols2 <- Reduce("intersect",lapply(mmf_CNV, function(x) colnames(x)))
        mmf_SNV <- do.call("rbind", lapply(mmf_SNV, function(x) x[,common_cols1]))
        mmf_CNV <- do.call("rbind", lapply(mmf_CNV, function(x) x[,common_cols2]))
        mmf_SNV <- as.data.frame(mmf_SNV)
        mmf_CNV <- as.data.frame(mmf_CNV)
        str(mmf_SNV)
        str(mmf_CNV)
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
        str(mmf_SNV)

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
