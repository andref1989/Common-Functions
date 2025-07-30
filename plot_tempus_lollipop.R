#' Plot lollipop plot of patients in Tempus data
#'
#' @param data_cohort Path to patient cohort of interest or the named list object output by load_tempus_data. The path option is slower and not advised if multiple oncoprints will be generated sequentially.
#' @param target_gene The gene of interest (string)
#' @param output_html The path to where you wish to store the final image. If none provided, it will default to file "(target_gene)_mutations.html" in the working directory. HTML output is currently required but will likely be optional in the future. (string)
#' @param variant_type_blacklist Mutation types to be excluded from the plots (typically those with benign functional impact) (string or vector) Default: "B"
#' @param title Plot title (Optional string)
#' @param assay_blacklist A list of assay versions you wish to exclude from consideration (xT v1,v2 etc) (Optional string)
#' @param filter_germline This is a boolean for whether you wish to exclude germline alterations from plotting.  (logical) Default: TRUE
#' @param return_df This is a boolean for whether you wish to return MAF file used to generate the figure.  (logical) Default: FALSE

#' @return Currently returns an interactive html file with your mutations in it. Optionally you can return a mutation file with the changes to the input required for plotting.
#' @note Any filtering or other grouping you wish to perform need to be done before passing the input cohort.
#' @note To embed an interactive version of the resulting html in RMarkdown reports, use the following:
#' xfun::file_string("path_to_output_html") or htmltools::includeHTML("path_to_output_html")
#'
#' @export
#' @examples
#' \dontrun{
#' plot_tempus_lollipop(data_cohort,target_gene="TP53",
#'   output_html="test.html",variant_type_blacklist=c("B","US"),
#'   gene_column = "gene_symbol",
#'   title = NULL, filter_germline=TRUE)
#'
#' }
#'

plot_tempus_lollipop <- function(data_cohort,target_gene,output_html=NULL,title=NULL,variant_type_blacklist="B",assay_blacklist=NULL,filter_germline=TRUE,return_df=FALSE){

##########################################################################################
############  Importing the tempus data and identifying the data model
##########################################################################################

  list_files_dm1 <- c("g_molecular_master_file")
  list_files_dm2 <- c("onco_result_snv_indel_passing",
                      "onco_result_cnv_gene")

    if(is.character(data_cohort)){

      # detect if the supplied cohort is DM 1 or DM 2
      files <- list.files(data_cohort,
                          recursive = TRUE)
      dm1_file_exists <- any(grepl(list_files_dm1, files))
      dm2_file_exists <- any(grepl(paste0(list_files_dm2, collapse = "|"), files))

      if(dm1_file_exists){

        input_td <- load_tempus_data(data_cohort,
                                     collection = NULL,
                                     list_files = list_files_dm1)
      } else if(dm2_file_exists){

        input_td <- load_tempus_data(data_cohort,
                                     collection = NULL,
                                     list_files = list_files_dm2)
      } else {     stop("The file path does not have the necessary files for either DM1 or DM2.")

        }

    } else if(is.list(data_cohort)){

        input_td <- data_cohort}

############  Checking the data model ##################

    data_model_check <- calc_data_model(input_td,
                                        list_tables_dm_1 = list_files_dm1,
                                        list_tables_dm_2 = list_files_dm2)

    stopifnot(data_model_check %in% c("1.0","2.0"))

############  Define DM specific functions ##################

    prep_lollipop_dm1 <- function(input_td){

        tempus_mmf <- input_td[["g_molecular_master_file"]]
        if(filter_germline){ tempus_mmf <- dplyr::filter(somatic_germline !="G")}
    if(!is.null(variant_type_blacklist)){
        tempus_mmf <- dplyr::filter(tempus_mmf, !functional_impact %in% variant_type_blacklist)
        if(nrow(tempus_mmf) <1){ print("Not enough mutations available after removing blacklisted types"); stop()}
}

    if(!is.null(assay_blacklist)){
        mutations <-  dplyr::filter(tempus_mmf, !assay %in% assay_blacklist,
                                    gene_canonical_name==target_gene,variant_type=="Short Variant")
    } else { mutations <-  dplyr::filter(tempus_mmf,
                                         gene_canonical_name==target_gene,variant_type=="Short Variant")
                                         }
    if(nrow(mutations) <1){ print(paste0("Not enough mutations available after removing ",paste0(assay_blacklist,collapse=","), "assays")); stop()}

    AA_translator <- data.frame(c("Ala","Arg","Asn","Asp","Cys","Gln",
                                  "Glu","Gly","His","Ile","Leu","Lys",
                                  "Met","Phe","Pro","Pyl","Ser","Sec",
                                  "Thr","Trp","Tyr","Val"),
                                c("A","R","N","D","C","Q","E","G","H",
                                  "I","L","K","M","F","P","O","S","U",
                                  "T","W","Y","V"),
                                c("Alanine","Arginine","Asparagine",
                                  "Aspartic acide","Cysteine","Glutamine",
                                  "Glumatic acid","Glycine","Histidine",
                                  "Isoleucine","Leucine","Lysine","Methionine",
                                  "Phenylalanine","Proline","Pyrolysine",
                                  "Serine","Selenocysteine","Threonine",
                                  "Tryptophan","Tyrosine","Valine"))



    colnames(AA_translator) <- c("Abbv","Letter","Full")
    rownames(AA_translator) <- AA_translator[,1]



    mutations$AA1 <- AA_translator[unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[1]))),"Letter"]

    mutations$AA2 <- AA_translator[unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[2]))),"Letter"]

####For translating multi-amino acid changes
    aa1_vec <- unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[1])))
    aa2_vec <- unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[2])))

    index1 <- which(nchar(aa1_vec) >3)
    index2 <- which(nchar(aa2_vec) >3)


    if(any(index1)){
        mutations$AA1[index1] <- unlist(lapply(index1, function(x) paste0(AA_translator[unlist(strsplit(gsub("(.{3})","\\1 ", aa1_vec[x])," ")),"Letter"], collapse="-")))}
    if(any(index2)){
        mutations$AA2[index2] <- unlist(lapply(index2, function(x) paste0(AA_translator[unlist(strsplit(gsub("(.{3})","\\1 ", aa2_vec[x])," ")),"Letter"], collapse="-")))}

####################################


    mutations$AA_Pos <- unlist(lapply(mutations$amino_acid_change, function(x) gsub("[A-Z]|[a-z]|\\.","",x)))

        mutations$AA_change <- paste0("p.",mutations$AA1,mutations$AA_Pos,mutations$AA2)

        mutations <- mutations[order(mutations$AA_Pos),]


        mutations$Mut_Class <- ifelse(mutations$mutation_effect=="synonymous_variant", "Silent",ifelse(mutations$mutation_effect=="missense_variant","Missense_Mutation","Unknown"))
        mutations$Gene <- mutations$gene_canonical_name
##    mutations$Mut_Class <- ifelse(mutations$AA1==mutations$AA2, "Silent", mutations$Mut_Class)

        write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)
        mutation_dat <- readMAF("mutations.txt", gene.symbol.col = "Gene", variant.class.col = "Mut_Class", protein.change.col = "AA_change",sep='\t')
        return(mutation_data)
    }
##### DM2 Function ####

    prep_lollipop_dm2 <- function(input_td) {

        tempus_mmf <- input_td[["onco_result_snv_indel_passing"]] %>% dplyr::filter(!grepl("UTR|intron", variant_molecular_consequence))
        if(nrow(tempus_mmf) <1){ stop("No coding mutations detected")}
        if(filter_germline){ tempus_mmf <- dplyr::filter(tempus_mmf, variant_origin !="germline") }

        if(!is.null(variant_type_blacklist)){

        tempus_mmf <- dplyr::filter(tempus_mmf, !variant_classification %in% variant_type_blacklist)
        if(nrow(tempus_mmf) <1){ stop("Not enough mutations available after removing blacklisted types")}
        }

    if(!is.null(assay_blacklist)){
        mutations <-  dplyr::filter(tempus_mmf, !assay %in% assay_blacklist,
                                    gene_symbol==target_gene,variant_type_detailed %in% c("snv","mnp"))

    } else { mutations <-  dplyr::filter(tempus_mmf,
                                         gene_symbol==target_gene,
                                         variant_type_detailed %in% c("snv","mnp"))}
 
    if(nrow(mutations) <1){ stop(paste0("Not enough mutations available after removing ",paste0(assay_blacklist,collapse=","), "assays"))}

            mutations$AA_change <- mutations$p_var
        mutations$AA_Pos <- unlist(lapply(mutations$p_var, function(x) gsub("[A-Z]|[a-z]|\\.","",x)))

        mutations <- mutations[order(mutations$AA_Pos),]
        mutations$Gene <- mutations$gene_symbol

    mutations$Mut_Class <- ifelse(mutations$variant_molecular_consequence=="synonymous_variant", "Silent",ifelse(mutations$variant_molecular_consequence=="missense_variant","Missense_Mutation",ifelse(grepl("splice", mutations$variant_molecular_consequence),"Splice","Unknown")))
##    mutations$Mut_Class <- ifelse(mutations$AA1==mutations$AA2, "Silent", mutations$Mut_Class)

        write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)
        mutation_dat <- readMAF("mutations.txt", gene.symbol.col = "Gene", variant.class.col = "Mut_Class", protein.change.col = "AA_change",sep='\t')
        return(mutation_dat)
    }

        ## mutations$AA1 <- unlist(lapply(mutations$p_var, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[1])))
        ## mutations$AA2 <- unlist(lapply(mutations$p_var, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[2])))



##################################################################



##########################################################################################
#### Calling formatting functions for oncoprint.
##########################################################################################

    if(data_model_check == "1.0"){

        mmf_int <- prep_lollipop_dm1(input_td)

    }
    else if (data_model_check =="2.0")  {
        mmf_int <- prep_lollipop_dm2(input_td)}



##    str(mutations)



######## Finalize plotting ###############


    plot_options <- g3Lollipop.options()
    if(!is.null(title)){ plot_options$titleText <- title} else{ plot_options$titleText <- paste0(target_gene, " mutation locations")}
    mutation_fig <- g3Lollipop(mmf_int,gene.symbol=target_gene,output.filename=title,gene.symbol.col="Gene", protein.change.col = "AA_change",btn.style="blue",plot.options=plot_options)
    if(!is.null(output_html)){
        htmlwidgets::saveWidget(mutation_fig,output_html)} else {
                                                             htmlwidgets::saveWidget(mutation_fig,paste0(target_gene,"_mutations.html"))}

    if(return_df){return(mutations)}
    }

