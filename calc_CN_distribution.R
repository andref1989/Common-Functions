calc_CN_distribution <- function(cohort_path,cohort_name, cohort_blacklist="xF"){

        if(any(grepl("parquet",list.files(cohort_path,"parquet",recursive=T)))){
            print("Parquet")
            td <- tryCatch({load_tempus_data(cohort_path, collection=NULL,
                               list_files=c("onco_result_cnv_gene"))
    }, error=function(f){"dud_cohort"})} else{
  td <- tryCatch({load_tempus_data(cohort_path, collection=NULL,
                                   list_files=c("g_molecular_master_file"))},
                 error=function(e){ tryCatch({load_tempus_data(cohort_path, collection=NULL, list_files=c("onco_result_cnv_gene"))
                 }, error=function(f){"dud_cohort"})})}


    td <- lapply(td, function(x) dplyr::filter(x, !assay %in% assay_blacklist))
    list_files_dm1 <- c("g_molecular_master_file")
    list_files_dm2 <- c("onco_result_cnv_gene")
    data_model <- calc_data_model(td, list_files_dm1, list_files_dm2)


    cnv <- calc_cnv(td)

    if(data_model =="1.0"){ measured <- dplyr::select(cnv, identifier=analysis_id,gene_canonical_name) %>% unique %>% data.frame
    } else if(data_model=="2.0"){ measured <- dplyr::select(cnv, identifier=tumor_biospecimen_id, gene_canonical_name) %>% unique %>% data.frame
    }



    ## measured <- as.data.frame(table(measured$gene_canonical_name))
    measured <- measured %>% group_by(gene_canonical_name) %>% dplyr::summarize(Num_Samples=n())
    saveRDS(measured,"~/test_measured.rds")

    cnv_out <- as.data.frame(table(dplyr::select(cnv,gene_canonical_name, copy_number))) %>% left_join(measured) %>% dplyr::filter(Freq >0)
    cnv_out$Fraction <- cnv_out$Freq/cnv_out$Num_Samples

    return(cnv_out)
}


calc_CN_distribution_v2 <- function(cohort_path,cohort_name,assay_blacklist=NULL){

  if(any(grepl("parquet",list.files(cohort_path,"parquet",recursive=T)))){
    print("Parquet")
    td <- tryCatch({load_tempus_data(cohort_path, collection=NULL,
                                     list_files=c("onco_result_cnv_gene"))
    }, error=function(f){"dud_cohort"})} else{
     td <-  tryCatch({vroom::vroom(paste0(cohort_path,"/Data/Group_Level_Molecular/g_molecular_master_file.csv"),delim=",", show_col_types=F) %>% dplyr::filter(variant_type_code=="CNALTER") %>%
        dplyr::select(patient_id, analysis_id, gene_canonical_name, variant_type_code, position_1,position_2, copy_number)},error=function(f) {tryCatch({load_tempus_data(cohort_path, collection=NULL, list_files=c("onco_result_cnv_gene"))
        }, error=function(e){"dud_cohort"})})

     if(is.data.frame(td)){ td <- list("g_molecular_master_file"=mmf)} else if(is.list(td)) { td <- td} else{ td <- "dud_cohort"}
  if(!is.null(assay_blacklist)){
  td <- lapply(td, function(x) dplyr::filter(x, !assay %in% assay_blacklist))}
  list_files_dm1 <- c("g_molecular_master_file")
  list_files_dm2 <- c("onco_result_cnv_gene")
  data_model <- calc_data_model(td, list_files_dm1, list_files_dm2)


  cnv <- calc_cnv(td)

  if(data_model =="1.0"){ measured <- dplyr::select(cnv, identifier=analysis_id,gene_canonical_name) %>% unique %>% data.frame
  } else if(data_model=="2.0"){ measured <- dplyr::select(cnv, identifier=tumor_biospecimen_id, gene_canonical_name) %>% unique %>% data.frame
  }



  measured <- measured %>% group_by(gene_canonical_name) %>% dplyr::summarize(Num_Samples=n())

  cnv_out <- as.data.frame(table(dplyr::select(cnv,gene_canonical_name, copy_number))) %>% left_join(measured) %>% dplyr::filter(Freq >0)
  cnv_out$Fraction <- round(cnv_out$Freq/cnv_out$Num_Samples,3)
  cnv_out$Cohort <- cohort_name
  return(cnv_out)
}
