plot_CCLE_lollipop <- function(CCLE_mut_df,variant_blacklist=c("Silent","Splice_Site"),target_gene,output_html=NULL,title=NULL, return_df=FALSE,gene_column="Hugo_Symbol",protein_change_column="Protein_Change",codon_column="Codon_Change"){
    require(g3viz)
    require(dplyr)

    if(!is.null(variant_blacklist)){
        mutations <- dplyr::filter(CCLE_mut_df, !Variant_Classification %in% variant_blacklist)
        if(nrow(mutations) <1){ print("Not enough mutations available after removing blacklisted types"); stop()}
}

    mutations <- dplyr::filter(mutations, get(gene_column)==target_gene)

    ## mutations$AA_Pos <- unlist(lapply(mutations[,protein_change_column], function(x) gsub("[A-Z]|[a-z]|\\.","",x)))

    ## mutations$AA1 <- unlist(lapply(mutations[,protein_change_column], function(x) gsub("p.|[0-9]{1,5}","", unlist(strsplit(x,">"))[1])))
    ## mutations$AA2 <- unlist(lapply(mutations[,protein_change_column], function(x) unlist(strsplit(x,">"))[2]))

    ## mutations$AA_change <- paste0("p.",mutations$AA1,mutations$AA_Pos,mutations$AA2)

    ## if(c(any(mutations$AA1 == "") ,any(mutations$AA2 == ""))){

    ##     index1 <- c(which(mutations$AA1==""), which(mutations$Protein_Change==""))
    ##     index2 <- c(which(mutations$AA2==""), which(mutations$Protein_Change==""))
    ##     if(!is.null(codon_column)){
    ##         codon1 <- unlist(lapply(index1,function(x) unlist(strsplit(mutations[x,codon_column],">"))[1]))
    ##         codon1 <- unlist(lapply(index1,function(x) unlist(strsplit(mutations[x,codon_column],">"))[1]))



    ##     }
    ## mutations <- mutations[order(mutations$AA_Pos),]

    codon_translator <- data.frame("Codon"=c("TAA","TAG","TGA","GCT","GCC","GCA","GCG","TGT","TGC","GAT","GAC","GAA","GAG","TTT","TTC","GGT","GGC","GGA","GGG","CAT","CAC","ATT","ATC","ATA","AAA","AAG","TTA","TTG","CTT","CTC","CTA","CTG","ATG","AAT","AAC","CCT","CCC","CCA","CCG","CAA","CAG","CGT","CGC","CGA","CGG","AGA","AGG","TCT","TCC","TCA","TCG","AGT","AGC","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG","TGG","TAT","TAC","del","fs","DEL","FS"),"Amino_Acid"=c("*","*","*","A","A","A","A","C","C","D","D","E","E","F","F","G","G","G","G","H","H","I","I","I","K","K","L","L","L","L","L","L","M","N","N","P","P","P","P","Q","Q","R","R","R","R","R","R","S","S","S","S","S","S","T","T","T","T","V","V","V","V","W","Y","Y","del","fs","del","fs"))
    rownames(codon_translator) <- codon_translator[,1]
################################################################################################################

    Mutations_Pos <- str_extract( mutations[,protein_change_column], "\\d+")


    index <- which(mutations[,protein_change_column] =="")
    if(any(index) && !is.null(codon_column)){
        codon_list <- mutations[index,codon_column]
        codon_loc <- unlist(lapply(codon_list, function(x) as.numeric(str_extract(x,"\\d+"))%/%3))+1
        AA1 <- lapply(codon_list, function(x) unlist(str_extract_all(unlist(strsplit(unlist(strsplit(x,">"))[1],"\\)"))[2],"[ATCGatcg]{3}|fs|del")))
        AA1 <- unlist(lapply(AA1, function(x) paste0(setdiff(codon_translator[toupper(x),2],NA),collapse="")))

        AA2 <- lapply(codon_list, function(x) unlist(str_extract_all(unlist(strsplit(x,">"))[2],"[ATCGatcg]{3}|fs|del")))
        AA2 <- unlist(lapply(AA2, function(x) paste0(setdiff(codon_translator[toupper(x),2],NA),collapse="")))
        AA2[is.na(AA2)] <- ""
        prot_change <- paste0("p.",AA1,codon_loc,AA2)
        mutations[index,protein_change_column] <- prot_change

        }

    mutations$Mut_Class <- mutations$Variant_Classification


    write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)



    mutation_dat <- readMAF("mutations.txt", gene.symbol.col = gene_column, variant.class.col = "Mut_Class", protein.change.col = protein_change_column,sep='\t')


    plot_options <- g3Lollipop.options()
    if(!is.null(title)){ plot_options$titleText <- title} else{ plot_options$titleText <- paste0(target_gene, " mutation locations")}
    mutation_fig <- g3Lollipop(mutation_dat,gene.symbol=target_gene,output.filename=title,gene.symbol.col=gene_column, protein.change.col = protein_change_column,btn.style="blue",plot.options=plot_options)
    if(!is.null(output_html)){
        htmlwidgets::saveWidget(mutation_fig,output_html)} else {
                                                             htmlwidgets::saveWidget(mutation_fig,paste0(target_gene,"_mutations.html"))}

    if(return_df){return(mutations)}
    }
