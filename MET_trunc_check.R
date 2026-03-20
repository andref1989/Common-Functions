library(dplyr)
library(reshape2)
library(Hmisc)
##library(diptest)
source("~/Andre_F_functions_git/Andre_F_functions.R")

devtools::load_all("~/Pathos_Projects/tempusR/")

cohorts <- system("pathostk cohort list",intern=TRUE)
cohorts2 <- lapply(cohorts, function(x) unlist(strsplit(x,"\t"))[-1])
cohort_df <- as.data.frame(do.call("rbind", cohorts2[-1]))
colnames(cohort_df) <- cohorts2[[1]]
cohort_df <- cohort_df[which(as.integer(cohort_df[,2]) <=2e4),]
##nrow(cohort_df)
for(i in 1:nrow(cohort_df)){
  if(!file.exists(paste0("/home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))){
    system(paste0("mkdir -p /home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))
    system(paste0("gcsfuse -o ro --only-dir ",cohort_df[i,1]," --implicit-dirs pathos-data /home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))
  }
}

pathos_dir <- "~/MDM2_TP53_eval/"

if(file.exists("~/characterization_check_CN.rds")){ characterization_check <- readRDS("~/characterization_check_CN.rds")} else{
  characterization_check <- data.frame(stringsAsFactors = F)}
for(i in setdiff(cohort_df[,1],characterization_check$Cohort)){
  print(paste0("Starting ",i))
  if(!file.exists(paste0(pathos_dir,i,"/Data"))){
    indir <- paste0(pathos_dir,i)
    cohort_path <- system(paste0("find ",indir," -wholename \'*/Data\' "),intern=TRUE)
  } else{ cohort_path <- paste0(pathos_dir,i)}

  int <- tryCatch({ calc_CN_distribution(cohort_path,i)}, error=function(e) {"dud cohort"})

  cohort <- gsub("/","_",i)
  saveRDS(int, paste0("~/Cohort_Characterization/copy_number/",cohort,"_characterization.rds"))

  int_df <- data.frame("Cohort"=i,"Status"="Finished")
  characterization_check <- rbind(characterization_check,int_df)
  saveRDS(characterization_check,"~/characterization_check_CN.rds")
  gc()
  print(paste0("Finished ",i))
}


extract_met_status <- function(cohort_path, cohort_name,assay_blacklist=NULL){
  if(any(grepl("parquet",list.files(cohort_path,"parquet",recursive=T)))){
    print("Parquet")
    td <- tryCatch({load_tempus_data(cohort_path, collection=NULL,
                                     list_files=c("onco_result_cnv_gene"))
    }, error=function(f){"dud_cohort"})} else{
      td <- tryCatch({load_tempus_data(cohort_path, collection=NULL,
                                       list_files=c("g_molecular_master_file"))},
                     error=function(e){ tryCatch({load_tempus_data(cohort_path, collection=NULL, list_files=c("onco_result_cnv_gene"))
                     }, error=function(f){"dud_cohort"})})
      td <- lapply(td, function(x) dplyr::filter(x,variant_type_code=="CNALTER"))
    }
  if(!is.null(assay_blacklist)){
    td <- lapply(td, function(x) dplyr::filter(x, !assay %in% assay_blacklist))}
  list_files_dm1 <- c("g_molecular_master_file")
  list_files_dm2 <- c("onco_result_cnv_gene")
  data_model <- calc_data_model(td, list_files_dm1, list_files_dm2)

}
