library(dplyr)
library(reshape2)
library(Hmisc)
library(diptest)
source("~/Andre_F_functions_git/Andre_F_functions.R")

devtools::load_all("~/Pathos_Projects/tempusR/")

cohorts <- system("pathostk cohort list",intern=TRUE)
cohorts2 <- lapply(cohorts, function(x) unlist(strsplit(x,"\t"))[-1])
cohort_df <- as.data.frame(do.call("rbind", cohorts2[-1]))
colnames(cohort_df) <- cohorts2[[1]]
cohort_df <- cohort_df[which(as.integer(cohort_df[,2]) <=7e3),]
##nrow(cohort_df)
for(i in 1:nrow(cohort_df)){
  if(!file.exists(paste0("/home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))){
    system(paste0("mkdir -p /home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))
    system(paste0("gcsfuse -o ro --only-dir ",cohort_df[i,1]," --implicit-dirs pathos-data /home/forbesa/MDM2_TP53_eval/",cohort_df[i,1]))
  }
}



pathos_dir <- "~/MDM2_TP53_eval/"

if(file.exists("~/characterization_check.rds")){ characterization_check <- readRDS("~/characterization_check.rds")} else{
  characterization_check <- data.frame(stringsAsFactors = F)}
for(i in setdiff(cohort_df[,1],characterization_check$Cohort)){
  print(paste0("Starting ",i))
  int <- tryCatch({ characterize_cohort(paste0(pathos_dir,i),i)}, error=function(e) {"dud cohort"})
  cohort <- gsub("/","_",i)
  source("~/Andre_F_functions_git/calc_treatment_journeys.R")

  int2 <- tryCatch( { calc_treatment_journeys(paste0(pathos_dir,i),num_groups=10,cohort_name=cohort)}, error=function(e) {"dud_cohort"})
  if(is.list(int2)){
    lapply(names(int2),function(x) saveRDS(int2[[x]], paste0("~/Cohort_Characterization/treatments/",cohort,"_",x,"_characterization.rds")))}
  if(is.list(int)){
  lapply(names(int),function(x) saveRDS(int[[x]], paste0("~/Cohort_Characterization/",x,"/",cohort,"_characterization.rds")))}
  if(file.exists(paste0("~/Cohort_Characterization/expression/",cohort,"_characterization.rds"))){
    int_df <- data.frame("Cohort"=i,"Status"="Finished")}
  characterization_check <- rbind(characterization_check,int_df)
  saveRDS(characterization_check,"~/characterization_check.rds")
  gc()
  print(paste0("Finished ",i))

}
