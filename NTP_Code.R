get_most_variable <- function(mat, margin=c('row','column'),cutoff=0,quantile=NULL){
    if(any(is.na(mat)) == TRUE){
        mat[is.na(mat)] <- 0
        print("There are NA values in the matrix. Working anyway")}

    else{ print("Working")}

    if(is.null(quantile)){
        cutoff=0 }
    else{ cutoff = quantile}

    print(cutoff)
    if( margin =="row"){
        vec <- apply(mat, 1, var)
        vec <- sort(vec,decreasing=TRUE)

        variable <- which(vec > quantile(vec,cutoff,na.rm=TRUE))


        #out <- mat[variable,]
    }
    else if( margin== "column"){
        vec <- apply(mat,2,var)
        variable <- which(vec > quantile(vec,cutoff))}
    return(variable)}


nearest_template_prediction_v3 <- function(reference_df,alt_df, tsne_df, cluster_column="Cluster",is_ranked=FALSE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,verbose=TRUE,random_sample=100,variable_length=TRUE,zscore_cutoff=1.4,templates=NULL, template_to_1=FALSE,sample_col="Sample") {
    require(lsa)
    require(dplyr)

    if(length(method)>1){ method <- "simple"} else{ method <- method}


    predicted_template_df <- data.frame(stringsAsFactors=F)

    template_df <- unique(tsne_df[unique(c(cluster_column))])
##    str(template_df)
    template_df <- apply(template_df,2, as.character)
##    gene_set <- unique(unlist(lapply(unique(tsne_df[,TF_column]), function(x) unlist(strsplit(x,"-")))))
    gene_set <- rownames(alt_df)


    if(is.null(templates)==TRUE){
    templates <- get_template_genes_v2(reference_df, tsne_df,cluster_column,quantile_cutoff=quantile_cutoff,method=method, template_length=template_length, is_ranked=is_ranked,value_returned="Templates",variable_length=variable_length, zscore_cutoff=zscore_cutoff)} else { templates <- templates}

    if(!all(unique(unlist(templates)) %in% union(rownames(reference_df), rownames(alt_df)))){ print("There are template genes missing from the query dataframe. Dropping missing genes")
    templates <- lapply(templates, function(x) intersect(x, intersect(rownames(reference_df),rownames(alt_df))))} else{ templates <- templates}
    ##    names(templates) <- sort(unique(template_df$Disease))
    print("Finished templates")
    ## str(templates)
###########
    samples <- lapply(sort(unique(template_df[,cluster_column])), function(x) dplyr::filter(tsne_df, get(cluster_column)==x)$Sample)
    names(samples) <- sort(unique(template_df[,cluster_column]))
    samples <- lapply(samples, function(x) intersect(x, colnames(reference_df)))

##    str(samples)

    sample_template_list <- lapply(names(samples), function(x) rowMeans(reference_df[templates[[x]],samples[[x]]],na.rm=TRUE))
##    str(sample_template_list)
    if(template_to_1==TRUE){
    for(x in 1:length(sample_template_list)){ int_ref <- sample_template_list[[x]]; int_ref[1:length(int_ref)] <- 1; sample_template_list[[x]] <- int_ref}} else { sample_template_list <- sample_template_list}

        names(sample_template_list) <- sort(unique(template_df[,cluster_column]))
    ##str(sample_template_list)
##    str(template_df)
##    str(sample_template_list)

#######################
#######################

    for(i in colnames(alt_df)){
##        str(i)
##        str(alt_df)

        distance_list <- unlist(lapply(sort(unique(template_df[,cluster_column])), function(x) 1-cosine(alt_df[templates[[x]],i],sample_template_list[[x]])))
        names(distance_list) <- sort(unique(template_df[,cluster_column]))

 ##       str(distance_list)

        random_list <- list()
        for(k in names(templates)){
  ##          print(k)
            random_list[[k]] <- unlist(lapply(1:random_sample, function(x) 1-cosine(alt_df[sample(gene_set, length(templates[[k]])),i],sample_template_list[[k]])))
  ##          str(random_list)

        }
        #lapply(random_list,function(x) print(summary(x)))

        pval_list <- unlist(lapply(names(distance_list), function(x) 1-pnorm(distance_list[[x]],mean(random_list[[x]]),sd(random_list[[x]]),lower.tail=F)))
        names(pval_list) <- names(distance_list)
        ###print(distance_list)
        ###print("Pvalue=")
        ###print(pval_list)
        candidate_template_df <- data.frame(pval_list,distance_list,names(distance_list),i,stringsAsFactors=F)
        colnames(candidate_template_df) <- c("P_val","Distance","Template","Cell_Line")
        candidate_template_df$Q_val <- p.adjust(candidate_template_df$P_val,method="BH")
        candidate_template_df$Signif <- candidate_template_df$Q_val <= 0.1
        sub_df <- dplyr::filter(candidate_template_df, Signif==TRUE)
##        str(candidate_template_df)
##        str(sub_df)
        if(nrow(sub_df) >=1){
            ##print("Working 1")
            candidate_template_df$Prediction <- sub_df[which(sub_df$Distance == min(sub_df$Distance)),"Template"]} else if (nrow(sub_df) == 0){
                ##print("Working 2")
                candidate_template_df$Prediction <- candidate_template_df[which(candidate_template_df$Distance == min(candidate_template_df$Distance,na.rm=TRUE)),"Template"]}
##        str(candidate_template_df)

##        str(random_list)
        predicted_template_df <- rbind(predicted_template_df,candidate_template_df)
                                        #str(sub_df)
        if(verbose==TRUE){
        print(paste0("Finished ", grep(paste0("^",i,"$"), colnames(alt_df)), " of ",ncol(alt_df)))}
    }

    output <- list(templates,predicted_template_df,sample_template_list)
    names(output) <- c("Templates","Predictions","Reference_Templates")
    print("Finished")
        return(output)
}


parse_NTP <- function(NTP_list, metadata_df,signif_cutoff=0.1,is_patient=FALSE){
    require(dplyr)

    if(is.data.frame(NTP_list)){

        NTP_list <- list(NTP_list)

        names(NTP_list) <- "NTP" } else{
            NTP_list <- NTP_list
    }

    for(i in names(NTP_list)){
        sub <- NTP_list[[i]]
        ##str(sub)
        sub$Signif <- sub$Q_val <= signif_cutoff
        sub_df <- dplyr::filter(sub, Signif==TRUE)
        for(k in unique(sub$Cell_Line)){
##            str(k)
            index <- which(sub$Cell_Line==k)
            int <- dplyr::filter(sub_df, Cell_Line==k)
##            str(int)
            if(nrow(int) >=1){
            ##print(int[which(int$Distance == min(int$Distance)),"Template"] == unique(sub[index,"Prediction"]))
            sub[index,"Prediction"] <- int[which(int$Distance == min(int$Distance)),"Template"]} else if (nrow(int) == 0){
#                print("Changing")
                sub[index,"Prediction"] <- sub[index,"Prediction"]




                                                                                               }}

        if(is_patient==FALSE){
        sub$Subtype <- metadata_df[sub$Cell_Line,"TCGA"]
        sub$Detailed_Cancer <- metadata_df[sub$Cell_Line,"lineage_subtype"]
        sub$Detailed_Cancer_Ext <- metadata_df[sub$Cell_Line,"disease_subtype"]
        sub$Primary <- metadata_df[sub$Cell_Line,"primary_or_metastasis"]
        sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff
        } else if(is_patient==TRUE){
            sub$Subtype <- metadata_df[sub$Cell_Line,"Subtype"]
            sub$Detailed_Cancer <- "Unknown"
            sub$Detailed_Cancer_Ext <- "Unknown"
            sub$Primary <- "Unknown"
            sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff

                                                                         } else if(is_patient =="skip"){
                                                                             sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
                                                                             sub$Q_val_cutoff <- signif_cutoff
                                                                             }

        NTP_list[[i]] <- sub
    }
    if(is_patient =="skip"){ print("Skipping accuracy assessment")
            output <- do.call("rbind",lapply(NTP_list, function(x) unique(x[c(4,6:ncol(x))])))
            final <- list(NTP_list,output)
            return(final)
            stop()} else{
##    str(NTP_list)
##str(NTP_list)

#print(colnames(output))

    output <- output%>% group_by(Subtype, Method, Primary) %>% mutate(Accuracy_fraction=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method) %>% mutate(Accuracy_fraction_overall=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method,Signif) %>% mutate(Accuracy_fraction_signif=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output$Primary <- gsub("^$","Unknown",output$Primary)
##output <- as.data.frame(output, stringsAsFactors=F)
##str(output)
##str(NTP_list)
final <- list(NTP_list,output)
##str(final)
    return(list(NTP_list,output))
        }}


get_template_genes_v2 <- function(expression_df, tsne_df, cluster_column="Cluster",is_ranked=TRUE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,value_returned=c("Templates","Matrix"),variable_length=TRUE,zscore_cutoff=1.8,abs_val=FALSE,sig_limit=1000){

    if(length(method)>1){ method <- "simple"} else{ method <- method}
    if(length(value_returned)>1){ method <- "Templates"} else { method <- method}
    require(dplyr)
    source("~/Andre_F_functions.R")

    sample_column <- grep("sample", colnames(tsne_df),ignore.case=TRUE, value=TRUE)
    tsne_df$Sample <- tsne_df[,sample_column]
    tsne_df <- tsne_df[which(tsne_df$Sample %in% colnames(expression_df)),]


    sub_df <- expression_df[,intersect(tsne_df$Sample, colnames(expression_df))]
    sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]

    if(method == "simple"){
        df <- zscore_df(do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample]))),"row")
##        str(df)

        df_abs <- do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample])))
        colnames(df) <- unique(tsne_df[,cluster_column])
        colnames(df_abs) <- unique(tsne_df[,cluster_column])
        if(is_ranked == TRUE){
            print("Working 1")
            if(variable_length==FALSE){
                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])} else if(variable_length==TRUE){

                    if(abs_val ==TRUE){ df <- abs(df)*-1


                                    } else{ df <- df}
                    templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]<=zscore_cutoff),x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))}

            names(templates) <- names(templates_abs)



        } else if (is_ranked==FALSE){
            print("Working 2")
            if(variable_length==FALSE){

                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                                if(abs_val ==TRUE){ df <- abs(df)} else{ df <- df}
                                          templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))




                                      }
        names(templates) <- colnames(df)

                                  }}

    else if(method == "randomforest"){
                                          require(randomForest)
        rownames(sub_df) <- gsub("-","_",rownames(sub_df))
        ##str(sub_df)
##        sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]
    transposed_df <- as.data.frame(t(sub_df))
    transposed_df[,cluster_column] <- as.factor(tsne_df[colnames(expression_df),cluster_column])
    uniq_classes <- unique(tsne_df[,cluster_column])
##    print(index)
##   print(cluster_column)


##    str(transposed_df)
    print("Starting Random Forest")
    rf <- randomForest(as.formula(paste0(cluster_column,"~.")),data=transposed_df,importance=TRUE)
        importance_df <- as.data.frame(importance(rf))[setdiff(colnames(importance(rf)),uniq_classes)]
##        str(importance_df)

    importance_df <- importance_df[order(-importance_df[,2]),]
    final_importance_matrix <- data.frame(importance(rf)[rownames(importance_df),uniq_classes],stringsAsFactors=F)


    colnames(final_importance_matrix) <- uniq_classes
##        str(final_importance_matrix)
##        str(colnames(final_importance_matrix))

    templates_abs <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, sort_order=TRUE))[1:sig_limit])

    templates_rel <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(zscore_df(final_importance_matrix),x,sort_order=TRUE))[1:sig_limit])
#    str(templates_rel)
        ##str(final_importance_matrix)
    names(templates_abs) <- colnames(final_importance_matrix)
    names(templates_rel) <- colnames(final_importance_matrix)

  ##      str(templates_abs)
  ##      str(templates_rel)
     if(variable_length==FALSE){
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                zscore_importance_matrix <- zscore_df(final_importance_matrix,"row",central_tendency="median")
                templates_rel <- lapply(colnames(zscore_importance_matrix), function(x) names(sort(zscore_importance_matrix[which(zscore_importance_matrix[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                names(templates_rel) <- colnames(zscore_importance_matrix)
##                str(templates_rel)
##                str(templates_abs)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))      }


#    templates <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, decreasing=TRUE))[1:template_length])
    templates <- lapply(templates, function(x) gsub("_","-",x))
        names(templates) <- names(templates_abs)


}
#    str(templates)
        templates <- lapply(templates, function(x) x[!is.na(x)])
    if(value_returned=="Templates"){ return(templates)} else if(value_returned=="Matrix"){ return(final_importance_matrix)}

}
