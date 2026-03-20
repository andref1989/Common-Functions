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

for(i in cohort_df[,1]){
  print(paste0("Starting ",i))
  cohort <- gsub("/","_",i)
  source("~/Andre_F_functions_git/calc_immune_infiltration.R")

  int3 <- tryCatch( { calc_immune_infiltration(paste0(pathos_dir,i),cohort_name=cohort)}, error=function(e) {"dud_cohort"})
  if(is.data.frame(int3)){
   saveRDS(int3, paste0("~/Cohort_Characterization/immune_infiltration/",cohort,"_immune_characterization.rds"))}

  gc()
  print(paste0("Finished ",i))

}
