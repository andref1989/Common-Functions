translate_ENSMBL_to_HGNC <- function(ENSEMBL_IDs,split=TRUE){

    require(tempusr)
    if(split==TRUE){ENSEMBL_IDs <- unlist(lapply(ENSEMBL_IDs, function(x) unlist(strsplit(x,"\\."))[1]))} else { ENSEMBL_IDs <- ENSEMBL_IDs}
    gene_df <- tempusr::gene_annotation
    rownames(gene_df) <- gene_df$ensembl_gene_id
    final_symbols <- gene_df[ENSEMBL_IDs,"hgnc_symbol"]
    return(final_symbols)}
