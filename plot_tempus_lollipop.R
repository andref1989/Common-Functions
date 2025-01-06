plot_lollipop <- function(tempus_mmf,target_gene,output_html,title=NULL,variant_type_blacklist=NULL,assay_subset=NULL,return_df=FALSE){
    require(g3viz)
    require(dplyr)

    if(!is.null(variant_type_blacklist)){
        tempus_mmf <- dplyr::filter(tempus_mmf, functional_impact %in% variant_type_blacklist)
        if(nrow(tempus_mmf) <1){ print("Not enough mutations available after removing blacklisted types"); stop()}
}

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

    mutations$Mut_Class <- ifelse(mutations$mutation_effect=="synonymous_variant", "Silent",ifelse(mutations$mutation_effect=="missense_variant","Missense_Mutation","Other"))

    write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)
##    str(mutations)


    mutation_dat <- readMAF("mutations.txt", gene.symbol.col = "gene_canonical_name", variant.class.col = "Mut_Class", protein.change.col = "AA_change",sep='\t')


    plot_options <- g3Lollipop.options()
    if(!is.null(title)){ plot_options$titleText <- title}
    mutation_fig <- g3Lollipop(mutation_dat,gene.symbol=target_gene,output.filename=title,gene.symbol.col="gene_canonical_name", protein.change.col = "AA_change",btn.style="blue",plot.options=plot_options)
    if(!is.null){
        htmlwidgets::saveWidget(mutation_fig,output_html)} else {
                                                             htmlwidgets::saveWidget(mutation_fig,paste0(target_gene,"_mutations.html"))}

    if(return_df){return(mutations)}
    }
