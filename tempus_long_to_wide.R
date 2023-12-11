tempus_rna_expression_to_mat <- function(rna_df,gene_identifier="Gene", gene_val_column="log2_gene_tpm_corrected", patient_identifier="patient_id"){
    require(reshape2)
    require(dplyr)
    gene_check <- any(grepl("^Gene$|^gene$", colnames(rna_df)))
    if(gene_check ==TRUE){ print("Not translating Ensembl genes to Symbols")} else{
                                                                                print("Translating gene codes to symbols")
                                                                                rna_df$Gene <- translate_ENSMBL_to_HGNC(rna_df$gene_code)}

    rna_mat <- reshape2::dcast(rna_df, as.formula(paste0(gene_identifier,"~", patient_identifier)),value.var=gene_val_column, fun.aggregate=mean)
    if(length(unique(rna_mat[,gene_identifier]))==nrow(rna_mat)){
        rownames(rna_mat) <- rna_mat[,gene_identifier]
        rna_mat[,gene_identifier] <- NULL
    } else{ print("Cannot set rownames to Gene symbols")}

    return(rna_mat)}
