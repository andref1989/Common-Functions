calc_CI <- function(data_cohort,amplification_threshold=5, deletion_threshold=1,assay_blacklist="xF"){
    ## if(!file.exists("~/CNA_Data.rds")){
    list_files_dm1 <- c("g_molecular_master_file")
    list_files_dm2 <- c("onco_result_cnv_gene")

    if(is.character(data_cohort)){


        data_cohort <- tryCatch({load_tempus_data(data_cohort, collection=NULL,list_files=list_files_dm2)}, error=function(e){ load_tempus_data(data_cohort, collection=NULL,list_files=list_files_dm1)})
        } else{ data_cohort <- data_cohort}

    gene_anno <- data_cohort[[1]] %>% dplyr::select(gene_symbol, gene_start_pos, gene_end_pos) %>% unique %>% dplyr::mutate(width=gene_end_pos-gene_start_pos)

    CNA_data <- tempusr::calc_cnv(data_cohort,assay_blacklist=assay_blacklist) %>% left_join(dplyr::select(gene_anno,gene_canonical_name=gene_symbol, width))
    ## saveRDS(CNA_data,"~/CNA_Data.rds")
## } else{
    ## CNA_data <- readRDS("~/CNA_Data.rds")
## }

    CNA_data <- CNA_data %>% group_by(tumor_biospecimen_id) %>% dplyr::mutate(Total_Sequenced=sum(width)) %>% group_by(tumor_biospecimen_id, gene_canonical_name) %>% dplyr::mutate(Altered=ifelse(copy_number>=amplification_threshold|copy_number<=deletion_threshold,TRUE, FALSE)) %>% group_by(tumor_biospecimen_id) %>% dplyr::mutate(Total_Width=sum(width),Altered_Width=as.numeric(Altered)*width) %>% dplyr::summarize(Total_Altered=sum(Altered_Width), Total_Sequenced=Total_Width) %>% unique %>% ungroup

    CNA_data$Altered_Fraction <- CNA_data$Total_Altered/CNA_data$Total_Sequenced

    return(CNA_data)

    }
