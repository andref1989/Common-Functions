identify_modules_from_KD <- function(directory="result_012",Network=NULL,gene="MDM2",to_file=TRUE,test=FALSE){
    library(googleCloudStorageR)
    library(gargle)
    output_list <- list()
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

#### set inputs ####

# set gene of interest - code requires both hgnc_symbol and ensembl_gene_id


gene_of_interest_hgnc_symbol <- gene

gene_of_interest_ensembl_gene_id <- tempusr::gene_annotation |>
  dplyr::filter(hgnc_symbol == gene_of_interest_hgnc_symbol) |>
  dplyr::pull(ensembl_gene_id)

# set p44 result folder to use


    p44_result_version <- directory
    sub_directory <- sprintf("project_0044/%s",p44_result_version)
    ## str(p44_result_version)
    ## str(sub_directory)

#### get data frame of all WGCN ####

# list all files
##networks <- list.files
networks <- googleCloudStorageR::gcs_list_objects(bucket = "pathos-research-results",
                                                  prefix = sub_directory)

# determine name of the network
    ## str(networks)

networks <- networks |>
  dplyr::mutate(network = stringr::str_split_i(name,
                                               "\\/",
                                               3)) |>
  dplyr::select(network) |>
  dplyr::distinct()

#### exclude CCLE networks

    if(is.null(Network)){
networks <- networks |>
  dplyr::filter(network != "CCLE")

} else{ networks <- dplyr::filter(networks, grepl(Network, network)|network==Network)}

#### determine which networks have a BN ####
    networks <- dplyr::filter(networks, grepl("[0-9]{2}Q[0-9]",network))

##### For testing only ######################
    if(test){
    networks <- as.data.frame(networks[c(1,sample(1:nrow(networks),1)),])
    colnames(networks)[1] <- "network"}
#############################################

## for each network
print("Working on the following networks:")
print(networks)
for (index in c(1:nrow(networks))) {

  # get path to BN

    
  path_BN <- googleCloudStorageR::gcs_list_objects(bucket = "pathos-research-results",
                                                   prefix = sprintf("project_0044/%s/%s/",
                                                                    p44_result_version,
                                                                    networks$network[index])) %>%
    dplyr::filter(stringr::str_detect(name,
                                      "BN_digraph_pruned_formatted"))

### determine if BN exists


  if (nrow(path_BN) >= 1) {

    has_BN <- TRUE

  } else {

    has_BN <- FALSE

  }
    networks$has_BN[index] <- has_BN


}

#### determine which BNs contain gene of interest and, if so, run KDA on modules ####

### function to parse txt from GCS


parse_txt <- function(object) {

  data <- httr::content(object,
                        as = "raw")

  data <- readr::read_tsv(data,
                          col_names = FALSE,
                          show_col_types = FALSE)

  colnames(data) <- c("source",
                      "target")

  return(data)

}

# add column to store whether or not gene of interest is in BN

networks$gene_in_BN <- NA

# create data frame to save KDA results

results_kda <- data.frame()

# for each network
##if(!file.exists(paste0(gene,"_",directory,"_results_kda.csv"))){
for (index in c(1:nrow(networks))) {

  # if there is a BN

  if (networks$has_BN[index]) {

    # load edges

    bn_edges <- pathosr::load_bn_data(sprintf("project_0044/%s/%s/tables/BayesianNetwork",
                                              p44_result_version,
                                              networks$network[index]))

    bn_edges <- bn_edges$BN_digraph_pruned_formatted

    # determine if gene of interest is in BN

    gene_in_BN <- any(gene_of_interest_ensembl_gene_id %in% union(bn_edges$source,
                                                                  bn_edges$target))

    # if the gene is in the BN

    if (gene_in_BN) {

      # print status

      print(sprintf("Running KDA for %s",
                    networks$network[index]))

      # load module gene membership for KDA

      module_gene <- googleCloudStorageR::gcs_get_object(bucket = "pathos-research-results",
                                                         object = sprintf("project_0044/%s/%s/tables/KDA/KDAInputFile.txt",
                                                                          p44_result_version,
                                                                          networks$network[index]),
                                                         parseFunction = parse_txt)

      # format module gene membership for KDA

      module_gene_formatted <- split(x = module_gene$source,
                                     f = module_gene$target)

 ## run KDA
if(!file.exists(paste0(directory,"_results_kda.csv"))){
      results_kda_temp <- pathosr::calc_key_drivers(network = as.matrix(bn_edges),
                                                    signatures = module_gene_formatted)

      # format and save KDA results

      results_kda <- dplyr::bind_rows(results_kda,
                                      results_kda_temp |>
                                        dplyr::mutate(network = networks$network[index]))

} else{ print("Reading in existing KDA file")
    results_kda <- read.csv(paste0(directory,"_results_kda.csv"))}

    # otherwise

  }} else {

    gene_in_BN <- NA

  }

  # save whether gene of interest was in BN

  networks$gene_in_BN[index] <- gene_in_BN

}



# save KDA results
if(to_file==TRUE & !file.exists(paste0(directory,"_results_kda.csv"))){
write.csv(x = results_kda,
          file = paste0(directory,"_results_kda.csv"),
          row.names = FALSE)
} else { output_list[["KDA"]] <- results_kda}


#### for all networks where BN contains gene of interest calculate module enrichment for genes associated with worse OS ####

# create data frame to store OS results
    print(networks)
    print(table(networks$gene_in_BN))
    networks <- dplyr::filter(networks, !is.na(gene_in_BN))
    if(any(networks$gene_in_BN) ){
results_os <- data.frame()
module_os_summary_all <- data.frame()
# for each network

for (index in c(1:nrow(networks))) {

        data_wgcn <- pathosr::load_wgcn_data(sprintf("project_0044/%s/%s",
                                                 p44_result_version,
                                                 networks$network[index]),
                                         list("os_signatures",
                                              "wgcna_gene_info"))

##      str(data_wgcn)
      if(!any(is.null(c(data_wgcn[c("os_signatures","wgcna_gene_info")])))){
    # format the OS signatures

    gene_info <- data_wgcn$os_signatures |>
      dplyr::mutate(effect = dplyr::case_when(estimate > 1 & adj.p < 0.05 ~ "worse_OS",
                                              estimate < 1 & adj.p < 0.05 ~ "better_OS",
                                              TRUE ~ "no_effect")) |>
      dplyr::select(ensembl_gene_id,
                    effect)


          ### Summarize modules
          module_os_summary <- left_join(dplyr::select(data_wgcn$os_signatures, Gene=ensembl_gene_id, hgnc_symbol,estimate),
                                         dplyr::select(data_wgcn$wgcna_gene_info, Gene=ensembl_gene_id,module_color,module_label),by="Gene") %>% group_by(module_color) %>% summarize_all(mean) %>% dplyr::select(module=module_color,Mean_Estimate=estimate) %>% dplyr::mutate(Mean_Effect = dplyr::case_when(Mean_Estimate > 1 ~ "worse_OS",
                                              Mean_Estimate < 1 ~ "better_OS",
                                              TRUE ~ "no_effect"))

          module_os_summary$network <- networks[index,"network"]
          module_os_summary_all <- rbind(module_os_summary_all,module_os_summary)
          ##str(module_os_summary)
if(to_file==TRUE){
          write.csv(module_os_summary_all,paste0(gene,"_",directory,"_module_os_summary.csv"),row.names=F)} else{ output_list[["OS_Summary"]] <- module_os_summary}



  # if gene of interest was in the BN

  if (!is.na(networks$gene_in_BN[index]) & networks$gene_in_BN[index]) {

    # load WGCN data

    data_wgcn <- pathosr::load_wgcn_data(sprintf("project_0044/%s/%s",
                                                 p44_result_version,
                                                 networks$network[index]),
                                         list("os_signatures",
                                              "wgcna_gene_info"))

##      str(data_wgcn)
      if(!any(is.null(c(data_wgcn[c("os_signatures","wgcna_gene_info")])))){
    # format the OS signatures

    gene_info <- data_wgcn$os_signatures |>
      dplyr::mutate(effect = dplyr::case_when(estimate > 1 & adj.p < 0.05 ~ "worse_OS",
                                              estimate < 1 & adj.p < 0.05 ~ "better_OS",
                                              TRUE ~ "no_effect")) |>
      dplyr::select(ensembl_gene_id,
                    effect)


          ### Summarize modules
          module_os_summary <- left_join(dplyr::select(data_wgcn$os_signatures, Gene=ensembl_gene_id, hgnc_symbol,estimate),
                                         dplyr::select(data_wgcn$wgcna_gene_info, Gene=ensembl_gene_id,module_color,module_label),by="Gene") %>% group_by(module_color) %>% summarize_all(mean) %>% dplyr::select(module=module_color,Mean_Estimate=estimate) %>% dplyr::mutate(Mean_Effect = dplyr::case_when(Mean_Estimate > 1 ~ "worse_OS",
                                              Mean_Estimate < 1 ~ "better_OS",
                                              TRUE ~ "no_effect"))

          module_os_summary$network <- networks[index,"network"]
          module_os_summary_all <- rbind(module_os_summary_all,module_os_summary)
          ##str(module_os_summary)
if(to_file==TRUE & !file.exists(paste0(gene,"_",directory,"_module_os_summary.csv"))){
          write.csv(module_os_summary_all,paste0(gene,"_",directory,"_module_os_summary.csv"),row.names=F)} else{ output_list[["OS_Summary"]] <- module_os_summary}
    # for each module

    list_modules <- unique(data_wgcn$wgcna_gene_info$module_color)

    for (index_module in list_modules) {

      # get the contingency table

      contingency <- data_wgcn$wgcna_gene_info |>
        dplyr::mutate(in_module = dplyr::case_when(module_color == index_module ~ "in_module",
                                                   TRUE ~ "not_in_module")) |>
        dplyr::select(ensembl_gene_id,
                      in_module) |>
        dplyr::inner_join(gene_info,
                          by = "ensembl_gene_id") |>
        dplyr::count(in_module,
                     effect,
                     name = "num_genes") |>
        dplyr::full_join(expand.grid(in_module = c("in_module",
                                                   "not_in_module"),
                                     effect = c("worse_OS",
                                                "better_OS",
                                                "no_effect")),
                         by = c("in_module",
                                "effect")) |>
        dplyr::mutate(num_genes = tidyr::replace_na(num_genes,
                                                    0)) |>
        tidyr::pivot_wider(names_from = "effect",
                           values_from = "num_genes") |>
        dplyr::arrange(dplyr::desc(in_module)) |>
        tibble::column_to_rownames("in_module") 

        contingency_worse <- dplyr::select(contingency,no_effect,worse_OS)
        contingency_better <- dplyr::select(contingency,no_effect,better_OS)

##        str(contingency_better)
##        str(contingency_worse)
      # calculate fishers exact

        results_os_temp_worse <- fisher.test(contingency_worse,alternative = "greater")
        results_os_temp_better <- fisher.test(contingency_better,alternative = "greater")

        pval_out <- c(results_os_temp_worse$p.value,results_os_temp_better$p.value)
        estimate_out <- c(results_os_temp_worse$estimate,results_os_temp_better$estimate)
      # save

      results_os <- dplyr::bind_rows(results_os,
                                     data.frame(network = networks$network[index],
                                                module = index_module,
                                                odds_ratio = estimate_out,
                                                p_value = pval_out,os_effect=c("worse","better")))}} 

   else { print("No OS signatures file found")}

      }
# adjust the p-value

results_os <- results_os |> 
  dplyr::group_by(network) |>
  dplyr::mutate(p_value_adj = p.adjust(p_value,
                                       "fdr")) |>
  dplyr::ungroup()

# save OS results
if(to_file==TRUE){
write.csv(x = results_os,
          file = paste0(gene,"_",directory,"_results_os.csv"),
          row.names = FALSE)
} else{ output_list[["OS"]] <- results_os}
#### format results for BD ####

### combine KDA (only gene of interest) and OS results

results_bd <- results_kda |>
  dplyr::filter(node == gene_of_interest_ensembl_gene_id) |>
  dplyr::select(network,
                module = signature,
                kda_odds_ratio = OR,
                kda_p_value_adj = padj) |>
  dplyr::full_join(results_os |>
                     dplyr::select(network,
                                   module,
                                   os_odds_ratio = odds_ratio,
                                   os_p_value_adj = p_value_adj),
                   by = c("network",
                          "module"))

# identify positive results
##    str(results_bd)
##    saveRDS(results_bd,"results_bd.rds")

results_bd <- results_bd |>
  dplyr::filter(kda_p_value_adj < 0.05,
                os_p_value_adj < 0.05) |>
  dplyr::group_by(network) |>
  dplyr::mutate(modules = paste0(module,
                                    collapse = ", "), kda_odds_ratios=paste0(round(kda_odds_ratio,3),collapse="~"), os_odds_ratios=paste0(round(os_odds_ratio,3),collapse="~")) |>
  dplyr::mutate(result = sprintf("%s is a key driver of the %s module(s), which is(are) also enriched for genes whose expression is associated with OS.",
                                 gene_of_interest_hgnc_symbol,
                                 modules)) |>
  dplyr::select(network,
                result,kda_odds_ratios,os_odds_ratios) %>% unique

# annotate non-positive results

##str(results_bd)
results_bd <- networks |>
  dplyr::mutate(result_no_BN = sprintf("There is no Bayesian network (yet) we can use to determine key drivers."),
                result_gene_not_in_BN = sprintf("%s is not in the Bayesian network (so we don't know what it may be a key driver of).",
                                                gene_of_interest_hgnc_symbol),
                result_none = sprintf("%s is not a key driver of any modules enriched for gene associated with worse outcomes.",
                                      gene_of_interest_hgnc_symbol)) |>
  dplyr::full_join(results_bd,
                   by = c("network")) |>
  dplyr::mutate(result = dplyr::case_when(!is.na(result) ~ result,
                                          has_BN == FALSE ~ result_no_BN,
                                          gene_in_BN == FALSE ~ result_gene_not_in_BN,
                                          TRUE ~ result_none)) |>
  dplyr::select(network,
                result,kda_odds_ratios,os_odds_ratios)

# save results
if(to_file==TRUE){
write.csv(x = results_bd,
          file = paste0(gene,"_",directory,"_results_bd.csv"),
          row.names = FALSE)
} else{ output_list[["BD"]] <- results_bd}
}}}
 else{ print("Gene wasn't present in any BN so not running anything... sorry")}}
