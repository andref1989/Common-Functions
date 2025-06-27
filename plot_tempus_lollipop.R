#' Plot Lollipop plots of mutations in a given gene in the Tempus data.
#'
#' @param tempus_mmf Data frame of sequential events (data.frame)
#' @param target_gene The number of sequential events to include in the plot for the cohort (Integer)
#' @param output_html The number of unique drug names/features/events to include in the plot. The plot will only include top "n" features(Integer)
#' @param title Plot title (Optional string)
#' @param variant_type_blacklist Whether to show the color legend for the plot or not.(Boolean)
#' @param assay_subset The column containing the sequential order of events (per individual data required) (string)
#' @param return_df  FALSE
#'
#' @return optional formatted data frame in MAF format of the mutations being plotted
#' @return Interactiv HTML file containing the final lollipop plots for your input mmf. This HTML can be placed directly into reports using the following code snippet. htmltools::includeHTML("mutations.html")
#' @note Any filtering or other grouping you wish to perform need to be done before passing the dataframe. Facetting has not been tested and likely will not work.
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' plot_tempus_lollipop(tempus_mmf,
#'   target_gene="TP53",
#'   output_html="test_mutations.html",
#'   title=NULL,
#'   title = "Test Plot",
#'   variant_type_blacklist=NULL,
#'   assay_subset=NULL,
#'   return_df=FALSE)
#'   }

plot_tempus_lollipop <- function(tempus_mmf,target_gene,output_html=NULL,title=NULL,variant_type_blacklist=NULL,assay_subset=NULL,return_df=FALSE){
    require(g3viz)
    require(dplyr)

    if(!is.null(variant_type_blacklist)){
        tempus_mmf <- dplyr::filter(tempus_mmf, !functional_impact %in% variant_type_blacklist)
    } else{ tempus_mmf <- dplyr::filter(tempus_mmf, !functional_impact %in% c("B","LB" )) }
        if(nrow(tempus_mmf) <1){ print("Not enough mutations available after removing blacklisted types"); stop()}


    if(!is.null(assay_subset)){
        mutations <-  dplyr::filter(tempus_mmf, assay %in% assay_subset,
                                    gene_canonical_name==target_gene,variant_type=="Short Variant")
    } else { mutations <-  dplyr::filter(tempus_mmf,
                                         gene_canonical_name==target_gene,
                                         variant_type=="Short Variant")}
    if(nrow(mutations) <1){ print(paste0("Not enough mutations available after subsetting to ",paste0(assay_subset,collapse=","))); stop()}

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

    mutations$AA_Pos <- unlist(lapply(mutations$amino_acid_change, function(x) gsub("[A-Z]|[a-z]|\\.","",x)))

    mutations$AA_change <- paste0("p.",mutations$AA1,mutations$AA_Pos,mutations$AA2)
    mutations <- mutations[order(mutations$AA_Pos),]


#### Default mutation mapping from Tempus:
    uniq_muts <- sort(unique(mutations$mutation_effect))
    short_uniq_muts <- unlist(lapply(uniq_muts, function(x) unlist(strsplit(x,"&"))[1]))

    mut_translator <- list("3_prime_UTR_variant"="3'UTR","5_prime_UTR_variant"="5'UTR","disruptive_inframe_deletion"="In_Frame_Del","disruptive_inframe_insertion"="In_Frame_Ins","frameshift_variant"="Frame_Shift","inframe_deletion"="In_Frame_Del","inframe_insertion"="In_Frame_Ins","initiator_codon_variant"="Start_Codon_SNP","missense_variant"="Missense_Mutation","splice_acceptor_variant"="Splice_Site","splice_donor_variant"="Splice_Site","splice_region_variant"="Splice_Region","start_lost"="Start_Codon_Del","stop_gained"="Nonsense_Mutation","stop_lost"="Nonstop_Mutation", "synonymous_variant"="Silent","upstream_gene_variant"="5'Flank", "downstream_gene_variant"="3'Flank")

    new_muts <- unlist(mut_translator[short_uniq_muts])
    names(new_muts) <- uniq_muts



    mutations$Mut_Class <- new_muts[mutations$mutation_effect]
    mutations$Mut_Class[is.na(mutations$Mut_Class)] <- "Unknown"

    write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)
##    str(mutations)


    mutation_dat <- readMAF("mutations.txt", gene.symbol.col = "gene_canonical_name", variant.class.col = "Mut_Class", protein.change.col = "AA_change",sep='\t')

    plot_options <- g3Lollipop.options()
    if(!is.null(title)){ plot_options$titleText <- title} else{ title=paste0(target_gene," SNVs")}

    mutation_fig <- suppressMessages(g3Lollipop(mutation_dat,gene.symbol=target_gene,gene.symbol.col="gene_canonical_name", protein.change.col = "AA_change",btn.style="blue",plot.options=plot_options))
    if(!is.null(output_html)){
        htmlwidgets::saveWidget(mutation_fig,output_html)} else {
                                                             htmlwidgets::saveWidget(mutation_fig,paste0(target_gene,"_mutations.html"))}

    if(return_df){return(mutations)}
}
