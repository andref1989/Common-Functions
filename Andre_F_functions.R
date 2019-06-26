ANF_qqplot <- function (pval_df,TF_list, regulating_TFs,Gene){
    pval_df$observed <- -log10(pval_df[,1])
    pval_df$expected <- -log10(runif(length(pval_df[,1])))
    pval_df$TF <- "No"
    pval_df$ID <- rownames(pval_df)
    TFs <- TF_list
    pval_df$TF[which(rownames(pval_df) %in% TFs)] <- "Yes"
    pval_df$TF[which(rownames(pval_df) %in% regulating_TFs)] <- paste0(Gene,"\nregulating")
    qplot <- ggplot(pval_df, aes(sample=observed))+stat_qq()
    qplot_df <- ggplot_build(qplot)$data[[1]]
    qplot_df$TF <- pval_df$TF[order(pval_df$observed)]
    qplot_df$ID <- pval_df$ID[order(pval_df$observed)]
    qplot_df$theoretical <- sort(-log10(runif(length(qplot_df$theoretical))))
    p <- ggplot(qplot_df,aes(theoretical,sample,label=ID))+geom_point(aes(color=qplot_df$TF))+geom_text(nudge_x=-0.2,size=1.5,aes(color=qplot_df$TF),check_overlap=TRUE)+geom_abline(linetype="dotted",slope=1,intercept=0)+geom_hline(aes(yintercept=-log10(0.05/length(qplot_df$TF))), color="coral", size=2)
    
    return(p)}

ascat_to_granges <- function(ascat_summary){ ascat_summary <- read.table(ascat_summary,header=T,stringsAsFactors=F,sep=',');ascat_summary[,1] <- gsub(23,"X",ascat_summary[,1]);ascat_summary[,1] <- gsub(24,"Y",ascat_summary[,1]);gr <- GRanges(seqnames=Rle(ascat_summary[,1]), ranges= IRanges(start=ascat_summary[,2], end=ascat_summary[,3]),strand="*",CN=ascat_summary[,6], Sample_ID=ascat_summary[,8]); gr <- gr.chr(gr);                               return(gr)}

bed_to_granges <- function(bed){ bed <- read.table(bed, stringsAsFactors=F, sep="\t")
                                 gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*")
                                 return(gr)}

bed_to_granges_enhancer <- function(bed){ bed <- read.table(bed, stringsAsFactors=F, sep="\t")
                                          gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",Target=bed$V4) ;return(gr)}
bed_to_granges_enhancer_extra <- function(bed){ bed <- read.table(bed, stringsAsFactors=F, sep="\t"); gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",Target=bed$V4, sample=bed$V5) ;return(gr)}


bed_to_granges_extra <- function(bed){ bed <- read.table(bed, stringsAsFactors=F, sep="\t")
                                 gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",sample=bed$V7) ;return(gr)}

gr_to_bed <- function(gr,outfile, metadata=FALSE){

                                   Chrom <- as.character(gr@seqnames)
                                   Start <- gr@ranges@start
                                   End <- gr@ranges@start+gr@ranges@width

                                   if(metadata==FALSE){
                                       df <- data.frame(Chrom,Start,End)}
                                   else if(metadata==TRUE){
                                       metadata <- as.data.frame(gr@elementMetadata@listData)
                                       df <- data.frame(Chrom,Start,End,metadata)
                                       }                                           
                                  
                                   write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=Tb)
                                   print(head(df))}

links_to_bed <- function(link_names, peak_file,outfile, metadata=FALSE){
    split1 <- link_names
    split1 <- unlist(strsplit(split1,"~"))
    peaks <- split1[seq(1,length(split1),2)]
    targets <- split1[seq(2,length(split1),2)]
    index <- which(peaks %in% names(peak_file))
    true_peaks <- peaks[index]
    true_targets <- targets[index]
    str(true_peaks)
    out_gr <- peak_file[true_peaks]
    Chrom <- as.character(out_gr@seqnames)
    Start <- out_gr@ranges@start
    End <- out_gr@ranges@start+out_gr@ranges@width
    if(metadata==FALSE){
        df <- data.frame(Chrom,Start,End)}
    else if(metadata==TRUE){
        metadata <- as.data.frame(out_gr@elementMetadata@listData)
        df <- data.frame(Chrom,Start,End,metadata)
    }
    df$Identifier <- true_peaks
    df$Target <- true_targets
    df <- df[order(df[,1],df[,2]),]
    write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=F)
    print(head(df))}

    
mutations_gr_to_bed <- function(gr,outfile){
    Chrom <- as.character(gr@seqnames)
    Start <- gr@ranges@start
    End <- gr@ranges@start+gr@ranges@width
    Ref <- gr$Ref
    Alt <- gr$Alt
    Gene <- gr$Gene
    df <- data.frame(Chrom,Start,End,Ref,Alt)

    write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=F)}    
                                   

    
Binarize_CN_table_subtype <- function(CN_table){
    CN_table$Binarized_Subtype <- 0
    CN_table[CN_table[,1] >6,"Binarized_Subtype"] <- 1
    return(CN_table)}


CN_to_granges <- function(CN_df){ CN_df <- read.table(CN_df,header=F,stringsAsFactors=F,sep='\t');CN_df[,1] <- gsub(23,"X",CN_df[,1]);CN_df[,1] <- gsub(24,"Y",CN_df[,1]); CN_df[,6] <- gsub("a|a_2|b|c","",CN_df[,6]);gr <- GRanges(seqnames=Rle(CN_df[,1]), ranges= IRanges(start=CN_df[,2], end=CN_df[,3]),strand="*",CN=CN_df[,5], Sample_ID=CN_df[,6]); gr <- gr.chr(gr);                               return(gr)}



Copy_Number_Eval <- function(CN_gr,interval_gr,expression_matrix,Gene){
    CNA_table <- as.data.frame(expression_matrix[,Gene])
    samples_intersecting_Gene <- intersect_with_metadata(CN_gr,interval_gr)
    seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
    
    rownames(CNA_table) <- rownames(expression_matrix)
    CNA_table$Copy_Number <- NA
    colnames(CNA_table) <- c(Gene,"Copy_Number")
    CNA_table$ID <-  rownames(CNA_table)
    for(i in rownames(CNA_table)){ index <- which(samples_intersecting_Gene$Sample_ID == i); CN <- mean(samples_intersecting_Gene$CN[index]);CNA_table[i,]$Copy_Number <- CN}
    CNA_table$Subtype <- NA; for(i in rownames(CNA_table)){sample_subtype <- expression_matrix[i,"Status"];CNA_table[i,]$Subtype <- sample_subtype}
return(CNA_table)}

Cross_validation <- function(matrix,num_of_folds,featurelist,Gene){
    fold_size <- round(dim(matrix)[1]/num_of_folds)
    out_list <- list()
    blacklist <- c()

    for(i in 1:num_of_folds){
        test <- sample(setdiff(rownames(matrix),blacklist),fold_size)
        blacklist <- c(blacklist,test)
        train <- setdiff(rownames(matrix),test)
        assign(paste0("model_",i),multiple_regression(matrix[train,],featurelist,Gene))
        test_result <- predict.glm(get(paste0("model_",i)),matrix[test,])
        df_out <- data.frame(matrix[test,Gene],test_result)
        colnames(df_out) <- c("Actual","Prediction")
        df_out$Features <- paste0(featurelist,collapse=",")
        out_list <- append(out_list,list(df_out))
        print(paste0("Finished Processing Fold ",i))
    }
    return(out_list)}


Condense_cross_validation_output <- function(Cross_validation_output,feature_info){
    final_df <- data.frame()
    cor_val_vec <- c()
    for(i in 1:length(Cross_validation_output)){
        sub <- Cross_validation_output[[i]]
        cor_val <- cor(sub[,1],sub[,2],use='complete.obs')^2
        Fold <- paste0("Fold_",i)

        out_df <- data.frame(cor_val,Fold,unique(sub$Features))
        final_df <- rbind(final_df,out_df)
    }
    colnames(final_df) <- c("R_Squared","Fold","Features")
    final_df$Feature_Type <- feature_info
   return(final_df) 
}



Plot_cross_validation_output <- function(Cross_validation_output,featurelist){
    final_df <- data.frame()
    cor_val_vec <- c()
    for(i in 1:length(Cross_validation_output)){
        sub <- Cross_validation_output[[i]]
        cor_val <- cor(sub[,1],sub[,2],use='complete.obs')^2
        Fold <- paste0("Fold_",i)

        out_df <- data.frame(cor_val,Fold,sub$Features)
        final_df <- rbind(final_df,out_df)
    }
    colnames(final_df) <- c("R_Squared","Fold","Features")
    final_df$Features <- paste0(featurelist,collapse="-")
    p <- ggplot(final_df,aes_string(x="Fold",y="R_Squared"))+geom_point()+geom_hline(yintercept=mean(final_df[,1]),size=2, alpha=0.5,color="coral", linetype=2)+labs(title=paste0("Model Performance in ",length(Cross_validation_output),"-fold cross-validation"),subtitle=final_df$Features)
    print(p)
    dev.off()

}

Plot_condensed_cross_validation_output <- function(condensed_output,fig_title){
    f_vec <- c()
    for(i in 1:length(condensed_output$Features)){
        f_len <- length(unlist(strsplit(as.character(condensed_output$Features[i]),",")))
        f_vec <- c(f_vec,f_len)}
    condensed_output$Model <- as.factor(f_vec)

    p <- ggplot(condensed_output, aes(x=Model,y=R_Squared,fill=Feature_Type))+geom_violin(trim=TRUE,draw_quantiles=c(0.25,0.5,0.75))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(title=fig_title,x="Number of model features")+geom_smooth(data=condensed_output,aes(x=as.numeric(Model),y=R_Squared))
    print(p)
    dev.off()
}

get_targeting_enhancers <- function(gr,metadata_column,gene){ out <- gr[which(gr@elementMetadata[,metadata_column] == gene)]; return(out)}

ID_enhancer_mutations <- function(enhancer_gr, Mutation_calls){
    infile <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    infile_gr <- GRanges(seqnames=Rle(c(infile$Chrom)), ranges=IRanges(start=infile$start, end=infile$end), Gene=infile$Gene, Ref=infile$Ref, Alt=infile$Alt, Sample=infile$Sample)
    mutations_gr <- infile_gr[findOverlaps(infile_gr,enhancer_gr)@from]
    return(mutations_gr)
}

ID_gene_mutations <- function(gene_of_interest, Mutation_calls){
    infile <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    infile_gr <- GRanges(seqnames=Rle(c(infile$Chrom)), ranges=IRanges(start=infile$start, end=infile$end), Gene=infile$Gene, Ref=infile$Ref, Alt=infile$Alt, Sample=infile$Sample)
    mutations_gr <- infile_gr[which(infile_gr$Gene == gene_of_interest),]
    return(mutations_gr)
}

intersect_with_metadata <- function(gr1,gr2, ignore.strand=FALSE){
    if(ignore.strand == FALSE){ out <- gr1[findOverlaps(gr1,gr2)@from,]} else if( ignore.strand ==TRUE){ out <- gr1[findOverlaps(gr1,gr2,ignore.strand=TRUE)@from,]}
    ; return(out)}

setdiff_with_metadata <- function(gr1,gr2){ int <- findOverlaps(gr1,gr2); out <- gr1[setdiff(1:length(gr1),unique(int@from)),]; return(out)}

Collapse_gr_keep_Gene <- function(gr){
    reduced <- reduce(gr)
    final_gr <- reduced
    for (i in 1:length(reduced)){
        overlap_enhancer = findOverlaps(reduced[i], gr)
        row_names = vector()
        for( j in 1:length(overlap_enhancer)){
            names = gr[overlap_enhancer[j]@to]$Gene
            row_names = c(names, row_names)
            
        }

        final_gr$nearby[i] = row_names
    }
    return(final_gr)}

make_enhancer_mutation_matrix <- function(targeting_enhancer_gr, Mutation_calls){ 
    infile_df <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile_df) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    rows <- length(unique(infile_df$Sample))
    columns <- length(targeting_enhancer_gr)
    infile_df$Chrom <- gsub("23","X",infile_df$Chrom)
    infile_df$Chrom <- gsub("24","Y",infile_df$Chrom)
    infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom)), ranges=IRanges(start=infile_df$start, end=infile_df$end), Gene=infile_df$Gene, Ref=infile_df$Ref, Alt=infile_df$Alt, Sample=infile_df$Sample)
    mutation_matrix <- matrix(data=0,nrow=rows,ncol=columns)
    rownames(mutation_matrix) <- unique(infile_df$Sample)
    colnames(mutation_matrix) <- names(targeting_enhancer_gr)
    #colnames(mutation_matrix) <- paste0(seqnames(targeting_enhancer_gr),":" ,ranges(targeting_enhancer_gr)@start,"-",ranges(targeting_enhancer_gr)@start+ranges(targeting_enhancer_gr)@width)
    
    seqlevelsStyle(infile_gr) <- seqlevelsStyle(targeting_enhancer_gr)
    
    
    #targeting_enhancer_gr <- reduce(targeting_enhancer_gr)
    for(i in 1:length(targeting_enhancer_gr)){ 
                                               mutations_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr[i])@from]
                                               if (length(mutations_gr) > 0){ out <- mutations_gr
                                                                              print(paste0(length(out)," mutations in ",i,"th enhancer"))
                                                                              mutation_matrix[out$Sample,i] <- 1}
                                               else if( length(mutations_gr) <= 0){print(paste0("No intersecting mutations in ", i,"th enhancer"))}
                                               
                                           }
                                              
return(mutation_matrix)

}

make_enhancer_sv_matrix <- function(targeting_enhancer_gr, Mutation_calls, complex_flag=FALSE,blacklist=NULL){
    infile_df <- data.frame(Mutation_calls$ID,Mutation_calls$variant_type,Mutation_calls$chr_from,Mutation_calls$chr_from_bkpt,Mutation_calls$chr_from_range,Mutation_calls$chr_from_strand,Mutation_calls$chr_to,Mutation_calls$chr_to_bkpt,Mutation_calls$chr_to_range,Mutation_calls$chr_to_strand, Mutation_calls$Score,Mutation_calls$variant_type ,stringsAsFactors=F)
    colnames(infile_df) <- c("ID","Type","Chrom1","Breakpoint1","Range1","Strand1","Chrom2","Breakpoint2","Range2","Strand2", "Score","Variant_Type")

    if(complex_flag== TRUE){
        infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom1)), ranges=IRanges(start=infile_df$Breakpoint1,end=(as.numeric(infile_df$Breakpoint1)+as.numeric(infile_df$Range1))),strand=infile_df$Strand1, Chrom2=infile_df$Chrom2, Breakpoint2=infile_df$Breakpoint2, Breakpoint_range2=infile_df$Range2,Score = infile_df$Score, strand2 = infile_df$Strand2 ,Sample=infile_df$ID, Type = infile_df$Variant_Type);    print(table(infile_gr$Variant_Type))
        

    }

    else {
        infile_df <- infile_df[which(infile_df$Chrom1 == infile_df$Chrom2),]
        intrachrom_var <- infile_df[which(infile_df$Variant_Type != "inversion"),]
        infile_df <- intrachrom_var        
        infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom1)), ranges=IRanges(start=infile_df$Breakpoint1, end=infile_df$Breakpoint2),strand=infile_df$Strand1, Chrom2=infile_df$Chrom2, Breakpoint2=infile_df$Breakpoint2, Breakpoint_range2=infile_df$Range2,Score = infile_df$Score, strand2 = infile_df$Strand2 ,Sample=infile_df$ID, Type = infile_df$Variant_Type)  }
    table(infile_df$Variant_Type)
    rows <- length(unique(infile_gr$Sample))
    columns <- length(targeting_enhancer_gr)
    mutation_matrix <- matrix(data=0,nrow=rows,ncol=columns)
    seqlevelsStyle(infile_gr) <- seqlevelsStyle(targeting_enhancer_gr)
    rownames(mutation_matrix) <- unique(infile_gr$Sample)
    colnames(mutation_matrix) <- names(targeting_enhancer_gr)
    
    if(is.null(blacklist) == FALSE){
        if(class(blacklist) == "GRanges"){
            print("Applying Blacklist")
            infile_gr <- setdiff_with_metadata(infile_gr,blacklist)
       }
                        else{ print("Blacklist needs to be in GenomicRanges format")}
    }
    else if(is.null(blacklist) == TRUE){ infile_gr <- infile_gr}
    temp_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr)@from]
    print(table(temp_gr$Type))
    
    
        for(i in 1:length(targeting_enhancer_gr)){ 
            mutations_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr[i])@from]

            
                                               if (length(mutations_gr) > 0){ out <- mutations_gr
                                                                              print(paste0(length(out)," mutations in ",i,"th enhancer"))
                                                                              pos <- out[which(out$Type != "deletion"),]
                                                                              neg <- out[which(out$Type == "deletion"),]
                                                                              mutation_matrix[neg$Sample,i] <- -1

                                                                              mutation_matrix[pos$Sample,i] <- 1
                                                                              
                                                                              
                                                                          }
                                               else if( length(mutations_gr) <= 0){print(paste0("No intersecting mutations in ", i,"th enhancer"))}
                                               
                                           }
                                              
return(mutation_matrix)

}

make_CDS_CN_matrix <- function(CN_gr, Gene, CDS_gr,padding=0,normalized=TRUE){
    gene_val <- length(grep(paste0("^",Gene,"$"),CDS_gr$Gene))
    if(gene_val > 0){
        interval_gr <- CDS_gr[grep(paste0("^",Gene,"$"), CDS_gr$Gene)]+padding
        names(interval_gr) <- paste0("exon_",1:length(interval_gr),"_CDS_CN")
        seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
        print(interval_gr)
        
                      
        CN_matrix <- matrix(0,ncol=length(interval_gr) ,nrow=length(unique(CN_gr$Sample_ID)))
        
        colnames(CN_matrix) <- names(interval_gr)
        rownames(CN_matrix) <- unique(CN_gr$Sample_ID)
        for(i in colnames(CN_matrix)){
            CN_intersect <- intersect_with_metadata(CN_gr,interval_gr[i])
            if(normalized==TRUE){
                for(j in unique(CN_intersect$Sample_ID)){
                    Sample_mean_CN <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                    Sample_sd_CN <- sd(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                    CN_val <- (mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN)-Sample_mean_CN)/Sample_sd_CN
                    CN_matrix[j,i] <- CN_val
                }}
            else if(normalized==FALSE){
                for(j in unique(CN_intersect$Sample_ID)){
                    CN_val <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN,na.rm=TRUE)
                    CN_matrix[j,i] <- CN_val
                }}
            print(paste0("Finished processing ",i))}

            

        } else {print("Cannot find this gene in database")}

    return(CN_matrix)
   
}
               
make_enhancer_CN_matrix <- function(CN_gr,interval_gr, blacklist=NULL,normalized=TRUE){

    seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
    if(is.null(blacklist) ==FALSE){
        print("Applying blacklist")
        seqlevelsStyle(blacklist) <- seqlevelsStyle(interval_gr)
        CN_gr <- setdiff_with_metadata(CN_gr,blacklist)
    } else{ print("No blacklist")}
    names(interval_gr) <- paste0(names(interval_gr),"_CN")
    print(interval_gr)
    CN_matrix <- matrix(0,ncol=length(interval_gr) ,nrow=length(unique(CN_gr$Sample_ID)))
        
        colnames(CN_matrix) <-names(interval_gr)
        rownames(CN_matrix) <- unique(CN_gr$Sample_ID)
        for(i in colnames(CN_matrix)){
            CN_intersect <- intersect_with_metadata(CN_gr,interval_gr[i])
            if(normalized == TRUE){
                for(j in unique(CN_intersect$Sample_ID)){
                Sample_mean_CN <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                
                Sample_sd_CN <- sd(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                
                CN_val <- (mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN)-Sample_mean_CN)/Sample_sd_CN
                
                CN_matrix[j,i] <- CN_val
            }}
            else if(normalized == FALSE){
                for(j in unique(CN_intersect$Sample_ID)){
                    CN_val <- mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN,na.rm=TRUE)
                    CN_matrix[j,i] <- CN_val
                }}
            print(paste0("Finished processing ",i))
        }

    return(CN_matrix)
   
}
                                                            
                             



calc_model_pval <- function(model){
    f <- summary(model)}
make_l1_enhancer_matrix <- function(reg_element_matrix,expression_matrix,Gene){
    reg_df <- as.data.frame(reg_element_matrix)
    reg_sub <- data.frame(reg_df[intersect(rownames(expression_matrix), rownames(reg_df)),])
    rownames(reg_sub) <- intersect(rownames(expression_matrix), rownames(reg_df))
    colnames(reg_sub) <- colnames(reg_df)
    reg_sub[,Gene] <- NA
    for(i in rownames(reg_sub)){ reg_sub[i,Gene] <- expression_matrix[i,Gene]}
    return(reg_sub)}



make_plot_df <- function(lm_input,Gene){
    out_df <- data.frame(lm_input$model[Gene],lm_input$fitted.values)
    colnames(out_df) <- c("Actual","Estimate")
    out_df$ID <- names(lm_input$fitted.values)
    subtype_column <- grep("Subtype|Status",colnames(lm_input$model))
    if(length(subtype_column) == 0){ out_df$Subtype <- "None"}
    else {out_df$Subtype <- lm_input$model[,subtype_column]}
    final_plot <- ggplot(out_df, aes(x=Actual,y=Estimate,label=ID))+geom_point(aes(color=out_df$Subtype))+geom_text(nudge_y=-0.1,aes(color=out_df$Subtype),size=1.5,check_overlap=TRUE)+geom_abline(linetype="dotted",slope=1,intercept=0)
    return(final_plot)

}

model_AUC <- function(matrix,feature_list,Gene){
    feature_list <- intersect(feature_list, colnames(matrix))
    out_df_final <- data.frame(0,0,0,0)
    colnames(out_df_final) <- c("Actual","Estimate","Feature_Count","Features")
    rownames(out_df_final) <- "Empty"

    for(i in 1:length(feature_list)){
        feature_list_sub <- sample(feature_list,i)
        formula_in <- paste0(feature_list_sub,collapse="+")
        final_formula <- paste0(Gene,"~",formula_in)
        lm_out <- glm(as.formula(final_formula),data=matrix)
        features <- paste(feature_list_sub,collapse="-")
        out_df <- data.frame(lm_out$model[Gene],lm_out$fitted.values,i,features)
        colnames(out_df) <- c("Actual","Estimate","Feature_Count","Features")
        rownames(out_df) <- paste0(names(lm_out$fitted.values),"_",i,"_features")
        out_df_final <- rbind(out_df_final,out_df)
    }
    out_df_final <- out_df_final[2:dim(out_df_final)[1],]
    out_df_final$Status <- NA
    subtype_col <- grep("Status|Subtype",colnames(matrix))
    for(i in rownames(matrix)){index <- grep(paste0(i,"_[0-9]*"),rownames(out_df_final))
                               out_df_final[index,"Status"] <- matrix[i,subtype_col]}

    return(out_df_final)}

        
feature_selection_experiment <- function(featurelist,matrix,Gene){
    feature_df_list <- list()
    featurelist <- intersect(featurelist,colnames(matrix))
    for (i in c(seq(1,length(featurelist),10),length(featurelist))){
        sample_size <- i
        readout_vec <- c()
        test_features_vec <- c()
        for(j in 1:100){ test_features <- sample(featurelist,sample_size)
                         model <- multiple_regression(matrix,test_features,Gene)
                         readout <- cor(model$model[Gene],model$fitted.values,use='complete.obs')^2
                         readout_vec <- c(readout_vec,readout)
                         test_features_vec <- c(test_features_vec,paste0(test_features,collapse="-"))}
        assign(paste0(i,"_feature_run"), data.frame(readout_vec,test_features_vec))
        feature_df_list <- append(feature_df_list,list(get(paste0(i,"_feature_run"))))
        print(paste0("Finished Processing run with ",i," Features out of a possible ",length(featurelist)))

    }
return(feature_df_list)
}

Collapse_feature_selection <- function(feature_selection_output,feature_description){
    out_df <- data.frame()
    for (i in 1:length(feature_selection_output)){
        sub_df <- feature_selection_output[[i]]
        sub_df$Feature_Count <- length(unlist(strsplit(as.character(sub_df[1,2]),"-")))
        sub_df$Description <- feature_description
        out_df <- rbind(out_df,sub_df)}
    colnames(out_df) <- c("R_Squared","Features","Feature_Count","Description")
    return(out_df)}
        


extract_best_features <- function(feature_selection_output, feature_description){
    feature_df <- data.frame()
    for(i in 1:length(feature_selection_output)){
        sub <- feature_selection_output[[i]]
        index <- which(sub[,1] == max(sub[,1],na.rm=TRUE))
        features <- as.character(sub[index,2])
        features <- paste0(features,collapse=",")
        out_df <- data.frame(features,feature_description)
        feature_df <- rbind(feature_df,out_df)}
    return(feature_df)}

Group_cross_validation <- function(best_features_output,matrix,num_folds=10,Gene,feature_info){
    in_df <- best_features_output
    final_df <- data.frame()
    for(i in 1:length(in_df[,1])){
        featurelist <- unique(unlist(strsplit(as.character(in_df[i,1]),"-|,")))
        cross_val <- Cross_validation(matrix,num_folds,featurelist,Gene)
        condensed <- Condense_cross_validation_output(cross_val,feature_info)
        final_df <- rbind(final_df,condensed)}
    return(final_df)}
    
        
glmnet_regression <- function(featurelist, matrix, Gene, family="gaussian",alpha){
    featurelist_sub <- intersect(featurelist,colnames(matrix))
    check_len <- length(featurelist)-length(featurelist_sub)
    if(check_len > 0){ print("Some of your provided features are *not* in the provided matrix")}
    featurelist <- featurelist_sub
    predictor <- as.matrix(matrix[,featurelist])
    
    response <- matrix[,Gene]

    fit <- glmnet(x=predictor,y=response,family=family,alpha=alpha)
    return(fit)}

cv_glmnet_regression <- function(featurelist, matrix, Gene, family="gaussian",alpha,nfolds=10, measure='mse'){
    featurelist_sub <- intersect(featurelist,colnames(matrix))
    check_len <- length(featurelist)-length(featurelist_sub)
    if(check_len > 0){ print("Some of your provided features are *not* in the provided matrix")}
    featurelist <- featurelist_sub
    predictor <- as.matrix(matrix[,featurelist])
    response <- matrix[,Gene]
    fit <- cv.glmnet(x=predictor,y=response,family=family,alpha=alpha,type.measure=measure,nfolds=nfolds)
    return(fit)}
run_best_glmnet <- function(cv_glmnet_object, matrix, Gene){

    cols <- coef(cv_glmnet_object)@Dimnames[[1]]
    sub_cols <- intersect(cols, colnames(matrix))
    best_lambda <- cv_glmnet_object$lambda.min
    index <- which(cv_glmnet_object$lambda == best_lambda)
    num_pos <- cv_glmnet_object$nzero[index]
    newx <- as.matrix(matrix[,sub_cols])
    fit <- predict(cv_glmnet_object,newx=newx, s= best_lambda)
    readout <- round(cor(matrix[,Gene], fit)^2,3)
    
    return(list(readout,num_pos))

}
plot_glmnet_regression <- function(glmnet_object,cross_validation_status=FALSE, outfile, title){
    if(cross_validation_status == FALSE){
        pdf(outfile)
        plot(glmnet_object, xvar='lambda')
        plot(glmnet_object, xvar='dev')
        dev.off()
    }
    if(cross_validation_status == TRUE){
        pdf(outfile)
        best_lambda <- glmnet_object$lambda.min
        best_index <- which(glmnet_object$glmnet.fit$lambda == glmnet_object$lambda.min)
        best_model <- glmnet_object$glmnet.fit$dev.ratio[best_index]
       
        plot(glmnet_object, ylim=c(0,10));abline(v= log(best_lambda),col="red", lty=2,main= title)
        plot(glmnet_object$glmnet.fit, xvar='lambda');abline(v= log(best_lambda),col="red", lty=2,main= title)
        plot(glmnet_object$glmnet.fit, xvar='dev', xlim= c(0,1)); abline(v= best_model,col="red", lty=2,main= title)
        dev.off()
    }
    return(best_model)
}

plot_case_study <-function(matrix,Gene,title,excluded_columns){
    muts <-matrix[setdiff(colnames(matrix),union(Gene,excluded_columns))]
    mut_rows <- c()
    for(i in rownames(muts)){ mut_row <- muts[i,]; count_element <- as.integer(mut_row !=0); mut_rows <- c(mut_rows, sum(count_element))}
    matrix$Num_mutations <- as.character(mut_rows)
    mut_status <- as.integer(mut_rows !=0)
    mut_status <- gsub(0,"No",mut_status)
    mut_status <- gsub(1,"Yes",mut_status)
 
    matrix$Mutation_Status <- mut_status
 
    matrix$Subtype <- expression_df[rownames(matrix),"Status"]
     

    color_change <- list("LumA"="Blue","LumB"="Green","Basal"="Red","Her2"="Purple","Normal"="White","Unknown"="Black")
    print(paste0("Plotting ", Gene))
    

    matrix_clean <- muts
    test <- names(table(unlist(unique(apply(matrix_clean,2,unique)))))
    print(test)

    check_cols <- names(table(unlist(unique(apply(matrix_clean,2,unique)))))
    color_pal <- list("-1"="Red","0"="White","1"="Green")


    #colfunc2 <- colorRamp2(breaks= sort(as.integer(check_cols)),colors = unlist(color_pal[check_cols]))



   
                                        #input_list <- sort(unique(matrix$Num_mutations))
    input_list <- names(which((table(matrix$Num_mutations) >1) == TRUE))
    
    comparison_list <- list(); for(i in 2:length(input_list)){ int <- list(c("0",input_list[i])); comparison_list <- append(comparison_list,int)}
    heights <- max(matrix[,Gene]) + seq(length.out= length(comparison_list))

                                        #matrix$Num_mutations <- mut_rows
    matrix[is.na(matrix$Subtype),"Subtype"] <- "Unknown"
    if((sum(abs(as.matrix(matrix_clean))) >0) && (dim(as.matrix(matrix_clean))[2] >1)){
        #colfunc2 <- colorRampPalette(colors= unlist(color_pal[check_cols]),bias=10.01, space='rgb')
        #heatmap.2(as.matrix(matrix_clean), col=colfunc2(length(unlist(color_pal[check_cols]))+10),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
                                        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
 #       heatmap.2(as.matrix(matrix_clean), col=colfunc2,RowSideColors=unlist(color_change[matrix$Subtype]), tracecol='white',main=title,cexRow=0.25,cexCol=1.25) 
        

       # heatmap.2(as.matrix(matrix_clean), col=colfunc2,RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25) 
        
        p <- ggplot(matrix,aes(x=Num_mutations,y=get(Gene),color=Num_mutations))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Number of Mutated elements",y=paste0(Gene," Expression"))+guides(color=FALSE); p <- p+geom_signif(test='t.test',comparisons=comparison_list,map_signif_level=FALSE,step_increase=0.1,tip_length=0,margin_top=0.4); print(p)
                                                                                                                                                                                          if(length(table(matrix$Num_mutations)) > 2){
            q <- ggplot(matrix,aes(x=Mutation_Status,y=get(Gene),color=Mutation_Status))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Has a mutated element",y=paste0(Gene," Expression")); print(q)} else{print("Not showing T/F for this gene")}
    } else {print("There are no mutations in the gene of interest")}


}

plot_case_study_redux <- function(matrix,Gene,Title,excluded_columns){
    muts <-matrix[setdiff(colnames(matrix),union(Gene,excluded_columns))]
    featurelist <- setdiff(colnames(matrix),excluded_columns)
    lm_out <- multiple_regression(matrix,featurelist,Gene)
    coeffs <- round(summary(lm_out)$coefficients[,1],digits=4)
    coeff_pvals <- round(summary(lm_out)$coefficients[,4],digits=4)

    print(coeffs)
    print(coeff_pvals)
    
    print(length(abs(coeffs)>0))
    coeff_names <- names(coeffs[which(abs(coeffs) >0)])
    submat <- muts[union(setdiff(coeff_names,"(Intercept)"),grep("_CDS_CN",colnames(matrix),value=TRUE))]
    mut_rows <- c()
    for(i in rownames(muts)){ mut_row <- muts[i,]; count_element <- as.integer(mut_row !=0); mut_rows <- c(mut_rows, sum(count_element))}
    matrix$Num_mutations <- as.character(mut_rows)
    mut_status <- as.integer(mut_rows !=0)
    mut_status <- gsub(0,"No",mut_status)
    mut_status <- gsub(1,"Yes",mut_status)

    matrix$Mutation_Status <- mut_status

    matrix$Subtype <- expression_df[rownames(matrix),"Status"]
    color_change <- list("LumA"="Blue","LumB"="Green","Basal"="Red","Her2"="Purple","Normal"="White","Unknown"="Black")
    print(paste0("Plotting ", Gene))
    
    matrix_clean <- muts
    check_cols <- names(table(unlist(applyhe(matrix_clean,2,unique))))
    print(check_cols)
    check_cols[which(as.numeric(check_cols) >=2)] <- "2"
    check_cols[which(as.numeric(check_cols) <=-2)] <- "-2"
    check_cols <- unique(check_cols)
    print("New")
    print(check_cols)
    color_pal <- list("-2"="Black","-1"="Red","0"="White","1"="Green","2"="Blue")

   if(length(check_cols) <=4){
                               breaks =seq(as.numeric(check_cols[1]),as.numeric(check_cols[length(check_cols)]),length.out=length(check_cols)+1)} else if(length(check_cols) ==5){ breaks=c(-2,-1,-0.5,0.5,1,2)}

#    breaks <-c(-2,-1,-0.5,0.5,1,2)
    print(breaks)
    print(unlist(color_pal[check_cols]))
    num_breaks=5

    expression_levels <- cut(expression_df[rownames(matrix),Gene],breaks=num_breaks)
    col_levels <- brewer.pal(num_breaks,"RdBu")
    col_list <- list();for(i in names(table(expression_levels))){col_list[[i]] <- col_levels[which(names(table(expression_levels))==i)]}
    rowCols <- unlist(col_list[expression_levels])






   
                                        #input_list <- sort(unique(matrix$Num_mutations))
    input_list <- names(which((table(matrix$Num_mutations) >1) == TRUE))
    
    comparison_list <- list(); for(i in 2:length(input_list)){ int <- list(c("0",input_list[i])); comparison_list <- append(comparison_list,int)}
    heights <- max(matrix[,Gene]) + seq(length.out= length(comparison_list))

                                        #matrix$Num_mutations <- mut_rows
    matrix[is.na(matrix$Subtype),"Subtype"] <- "Unknown"
    if((sum(abs(as.matrix(matrix_clean))) >0) && (dim(as.matrix(matrix_clean))[2] >1)){
       #colfunc2 <- colorRampPalette(colors= unlist(color_pal[check_cols]),bias=10.01, space='rgb')
        #heatmap.2(as.matrix(matrix_clean), col=colfunc2(length(unlist(color_pal[check_cols]))+10),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
                                        #heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)

        
        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]),trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256), keysize=1, key.title=NA, key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=rowCols,trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(submat), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]),trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(submat), col=unlist(color_pal[check_cols]),RowSideColors=rowCols,trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')

        
        p <- ggplot(matrix,aes(x=Num_mutations,y=get(Gene),color=Num_mutations))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Number of Mutated elements",y=paste0(Gene," Expression"))+guides(color=FALSE); p <- p+geom_signif(test='t.test',comparisons=comparison_list,map_signif_level=FALSE,step_increase=0.1,tip_length=0,margin_top=0.4); print(p)
      if(length(table(matrix$Num_mutations)) > 2){                                                                                                                                                                      
            q <- ggplot(matrix,aes(x=Mutation_Status,y=get(Gene)))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Has a mutated element",y=paste0(Gene," Expression"))+geom_signif(test='t.test', comparisons=list(c("No","Yes")),map_signif_level=FALSE,tip_length=0); print(q)} else{print("Not showing T/F for this gene")}}             else {print("There are no mutations in the gene of interest")}


}


    
multiple_regression <- function(matrix,feature_list,Gene){
    feature_list <- intersect(feature_list, colnames(matrix))

    formula_in <- paste0(feature_list,collapse="+")
    final_formula <- paste0(Gene,"~",formula_in)
    lm_out <- glm(as.formula(final_formula), data=matrix,maxit=100)
    return(lm_out)}

make_conting_table_model_AUC <- function(model_AUC_output,num_steps=20){
    conting_table_final <- data.frame(0,0,0,0,0,0,0,"Empty")
    colnames(conting_table_final) <- c("TP","FP","TN","FN","Feature_Count","TPR","FPR","Features")
    for(i  in unique(model_AUC_output$Feature_Count)){
        in_df <- model_AUC_output[which(model_AUC_output$Feature_Count == i),]
        step_counter <- seq.int(min(in_df$Estimate),max(in_df$Estimate),length.out=num_steps)
       
        TP_vec <- c()
        FP_vec <- c()
        TN_vec <- c()
        FN_vec <- c()
       
        readout_list <- list("11"="TP","10"="FP","01"="FN","00"="TN")

        for (j in step_counter){
            Prediction <- as.integer(in_df$Estimate > j)
            Actual <- as.integer(in_df$Status %in% c("LumA","LumB"))
            Comparison <- as.integer(Prediction == Actual)
    
          


            readout <- paste0(Prediction,Actual)
            readout_values <- table(unlist(readout_list[readout]))
            TP_vec <- as.vector(c(TP_vec,readout_values["TP"]))
            FP_vec <- as.vector(c(FP_vec,readout_values["FP"]))
            TN_vec <- as.vector(c(TN_vec,readout_values["TN"]))
            FN_vec <- as.vector(c(FN_vec,readout_values["FN"]))
             
            }
        conting_table <- data.frame(TP_vec,FP_vec,TN_vec,FN_vec,i)
        colnames(conting_table) <- c("TP","FP","TN","FN","Feature_Count")
        conting_table[is.na(conting_table)] = 0
        conting_table$TPR <- conting_table$TP/(conting_table$TP+conting_table$FN)
        conting_table$FPR <- 1 -(conting_table$TN/(conting_table$TN+conting_table$FP))
        conting_table$Features <-  as.character(unique(in_df$Features))

                      
        
        conting_table_final <- rbind(conting_table_final,conting_table)
        
    }
    conting_table_final <- conting_table_final[2:dim(conting_table_final)[1],]
    conting_table_final$Features <- as.character(conting_table_final$Features)

    return(conting_table_final)}
 
plot_conting_table <- function(conting_table_AUC, outfile){
    pdf(outfile)
    if(length(unique(conting_table_AUC$Feature_Count)) > 5 ){
        subsample <- sample(unique(conting_table_AUC$Feature_Count),5)}
    else if(length(unique(conting_table_AUC$Feature_Count)) < 5){
        subsample <- sample(unique(conting_table_AUC$Feature_Count),length(unique(conting_table_AUC$Feature_Count)))}
    
    conting_table_sub <- conting_table_AUC[which(conting_table_AUC$Feature_Count %in% subsample),]
    AUC_vec <- c()
    for(i in conting_table_sub$Feature_Count){sub <- conting_table_sub[which(conting_table_sub$Feature_Count == i),];     AUC <- integrate.xy(sub$FPR,sub$TPR); AUC_vec <- c(AUC_vec,AUC)}
    AUC_vec <- round(AUC_vec,digits=3);str(AUC_vec);str(conting_table_sub)
    conting_table_sub$Feature_Count <- paste0(conting_table_sub$Feature_Count,":AUC=",AUC_vec)
    p <- ggplot(conting_table_sub,aes(x=FPR,y=TPR,color=Feature_Count))+geom_line()+geom_abline(slope=1,intercept=0,linetype="dotted")
    q <- p+labs(color="Feature Count",caption= gsub(".pdf","",outfile))
    print(q)

    dev.off()
    out <- gsub(".pdf",".txt",outfile)
    write.table(conting_table_sub,out,sep='\t',quote=F,row.names=F)
    return(conting_table_sub)
}

Plot_vector_as_bar_chart <- function(vector,outfile,x_field,y_field,negative_axis=FALSE,return_df=FALSE){
    df <- data.frame(vector)
    df$Name <- names(vector)
    colnames(df) <- c("Value","Name")

    str(df)
    vector_mean <- mean(abs(df$Value))
    print(vector_mean)
    if(negative_axis ==TRUE){
        p <- ggplot(df,aes(x=Name,y=Value,label=Name))+labs(x=x_field,y=y_field)+geom_col(alpha=0.65)+geom_text(nudge_y=0.01,size=3,check_overlap=TRUE)+theme(axis.text.x=element_text(face='bold',size=3,angle=45))+geom_hline(yintercept=vector_mean,linetype='dotted',col="salmon",size=2)+geom_hline(yintercept=-1*vector_mean,linetype='dotted',col="salmon",size=2)}
    else if(negative_axis==FALSE){
        p <- ggplot(df,aes(x=Name,y=Value,label=Name))+labs(x=x_field,y=y_field)+geom_col(alpha=0.65)+geom_text(nudge_y=0.01,size=3,check_overlap=TRUE)+theme(axis.text.x=element_text(face='bold',size=3,angle=45))+geom_hline(yintercept=vector_mean,linetype='dotted',col="salmon",size=2)}
    pdf(outfile)
    print(p)
    dev.off()
    if(return_df == TRUE){ return(df)} else{print("Finished")}
}


prepare_matrix <- function(matrix, gene) {
    pheno_input <- matrix[,gene]
    cleaned_matrix <- matrix
    cleaned_matrix[,gene] <- NULL
    return(list(pheno_input,cleaned_matrix))
}

pval_calculator <- function(pheno_input, x_input ){

    
    X_mx <- as.matrix(cbind(1,x_input))
    n_samples <- length(x_input)[1]
    
    
    MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input

    y_hat <- X_mx %*% MLE_beta


    SSM <- sum((y_hat - mean(pheno_input))^2)
    SSE <- sum((pheno_input - y_hat)^2)

    df_M <- 2
    df_E <- n_samples - 3
    MSM <- SSM / df_M
    MSE <- SSE / df_E
    

    Fstatistic <- MSM / MSE
    pval <- pf(Fstatistic,df1 = 2, n_samples-3,lower.tail = FALSE)



    return(pval)
}

significant_pval_ID <- function(pval_df){ cutoff <- 0.05/length(pval_df[,1]); significant_ids <- which(pval_df[,1] <= cutoff); significant <- pval_df[significant_ids,]; names(significant) <- rownames(pval_df)[significant_ids]; return(significant)}

table_to_granges <- function(table){ gr <- GRanges(seqnames=Rle(table[,1]), ranges= IRanges(start=table[,2], end=table[,3]),strand="*"); return(gr)}

tidy_matrix <- function (matrix){
    matrix <- as.data.frame(lapply(matrix,function(y) as.numeric(gsub(-Inf,-10000,y))))
    matrix[is.na(matrix)] <- 0
    return(matrix)

}

violin_plots <- function(matrix, list_of_features, group_column=NULL,outfile){
    pdf(outfile)
    featurelist <- sort(intersect(list_of_features,colnames(matrix)))
    str(featurelist)
    str(list_of_features)
    print(setdiff(list_of_features,featurelist))
    if(is.null(group_column) == FALSE){group_column <- intersect(group_column, colnames(matrix))}
    else {group_column <- group_column}

    submatrix <- matrix[,c(featurelist,group_column)]
    for (i in featurelist){ p <- ggplot(submatrix, aes(x=get(group_column),y=get(i), fill=get(group_column)))+geom_violin(trim=FALSE)+labs(fill=group_column)+xlab(group_column)+ylab(i)
                            q <- p+stat_summary(fun.data="mean_sdl", geom="crossbar",width =0.04)+geom_signif(test='t.test',comparisons = list(c("LumA","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+1.5)+geom_signif(test='t.test',comparisons = list(c("LumB","Normal")), map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+1) + geom_signif(test='t.test',comparisons= list(c("Her2","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2) + geom_signif(test='t.test',comparisons = list(c("Basal","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2.5)+ geom_signif(test='t.test',comparisons = list(c("Yes","No")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2.5)
                            print(q)
                            print(i)}
    dev.off()}


    
get_CN_expression_correlation <- function(genelist,matrix){
    cor_vec <- c()
    for(i in genelist){ expression <- matrix[i]
                        CN <- get(paste0("CNA_table_",i))$Copy_Number
                        cor_val <- cor(expression,CN,use='complete.obs')
                        cor_vec <- c(cor_vec,cor_val)}
    names(cor_vec) <- genelist
    return(cor_vec)}


Copy_Number_Eval_all <- function(genelist, expression_matrix, Gene, CN_gr, CDS_gr,padding=0 ){
    CNA_table_final <- as.data.frame(expression_matrix[,Gene])
    colnames(CNA_table_final) <- Gene
    for(i in genelist){
        gene_val <- length(grep(paste0("^",i,"$"),CDS_gr$Gene))
        if(gene_val > 0){ interval_gr <- CDS_gr[grep(paste0("^",i,"$"), CDS_gr$Gene)[1]]+padding
                          samples_intersecting_Gene <- intersect_with_metadata(CN_gr,interval_gr)
                          seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
                          rownames(CNA_table) <- rownames(expression_matrix)
                          CNA_table$Copy_Number <- NA
                          colnames(CNA_table) <- c(Gene,"Copy_Number")
                          CNA_table$ID <-  rownames(CNA_table)
                          for(j in rownames(CNA_table)){ index <- which(samples_intersecting_Gene$Sample_ID == j)
                                                         CN <- mean(samples_intersecting_Gene$CN[index])
                                                         CNA_table[j,"Copy_Number"] <- CN }
                          CNA_table$Subtype <- NA
                          for(k in rownames(CNA_table)){
                              sample_subtype <- expression_matrix[k,"Status"]
                              CNA_table[k,]$Subtype <- sample_subtype}
                          colnames(CNA_table) <- c(Gene,paste0(i,"_Copy_Number"),"ID","Subtype")
                          CNA_table_final <- cbind(CNA_table_final, CNA_table[,paste0(i,"_Copy_Number")])
                          col_names <- colnames(CNA_table_final)
                          index <- grep("CNA_table",col_names)
                          col_names[index] <- paste0(i,"_CN")
                          colnames(CNA_table_final) <- col_names
                          print(paste0("Finished ",i))
                          
                      }
        else if(gene_val ==0){
                                        #print(paste0("No CDS for ",i))
        }
        
    }
    CNA_table_final$ID <-CNA_table$ID
    CNA_table_final$Subtype <- CNA_table$Subtype
return(CNA_table_final)
}


get_regulator_TFs <- function(graph,Gene){
    index <- which(names(V(graph)) == Gene)
    regulators <- names(which(graph[,index]!=0))
    return(regulators)}

get_regulatory_targets <- function(graph,TF){
    index <- which(names(V(graph)) == TF)
    targets <- names(which(graph[index,]!=0))
    return(targets)}

delete.na <- function(DF, n=0) {
    DF[rowSums(is.na(DF)) <= n,]
}
split_vector <- function(vector, split_size){ split(vector, ceiling(seq_along(vector)/split_size))}

RIGHT <- function(x,n){
    substring(x,nchar(x)-n+1)
}
get_hub_TFS <- function(graph,q=0.75){
    outdegree <- igraph::degree(graph,mode='out')
    outdegree <- outdegree[which(outdegree >0)]
    hubs <- sort(names(outdegree[which(outdegree >= quantile(outdegree,q))]))
    return(hubs)}
    
get_common_hub_TFs <- function(hublist, graphs_to_exclude) {
    sub_list <- hublist[setdiff(names(hublist),graphs_to_exclude)]
    common_vec <- sub_list[[1]]
    for(i in 2:length(sub_list)){
        common_vec <- intersect(common_vec,sub_list[[i]])}
    blacklist <- unique(unlist(hublist[graphs_to_exclude]))
    common_vec <- setdiff(common_vec,blacklist)
    return(common_vec)}


importJaspar <- function(file=myloc) {
    vec <- readLines(file)
    vec <- gsub("\t"," ",vec)
    vec <- gsub("\\[|\\]", "", vec)
    start <- grep(">", vec); end <- grep(">", vec) - 1
    pos <- data.frame(start=start, end=c(end[-1], length(vec)))
    pwm <- sapply(seq(along=pos[,1]), function(x) vec[pos[x,1]:pos[x,2]])
    pwm <- sapply(seq(along=pwm), function(x) strsplit(pwm[[x]], " {1,}"))
    pwm <- lapply(seq(along=start), function(x) matrix(as.numeric(t(as.data.frame(pwm[(pos[x,1]+1):pos[x,2]]))[,-1]), nrow=4, dimnames=list(c("A", "C", "G", "T"), NULL)))
    names(pwm) <- gsub(">", "", vec[start])
    return(pwm)
}

Jaspar2meme <- function(encode_motifs_file, outdir){
    in_motifs <- importJaspar(encode_motifs_file)
    out_motifs <- list(); for(i in 1:length(in_motifs)){ motif <- t(in_motifs[[i]])
                                                         motif <- motif/rowSums(motif)
                                                         motif_name <- gsub(" ","~",names(in_motifs[i]))
                                                         outlist <- list(motif)
                                                         outname <- paste0(outdir, motif_name,"_MEME.txt")
                                                         
                                                         writeLines("MEME version 4", outname)
                                                         write("ALPHABET= ACGT",outname,append=T)
                                                         write("stands: +-",outname,append=T)
                                                         write(paste0("MOTIF ", motif_name), outname,append=T )
                                                         write(paste0("letter-probability matrix: alength=4 w= ",dim(motif)[1]),outname,append=T)
                                                         vec <- c()
                                                         for(j in 1:dim(motif)[1]){
                                                             out_vec <- paste(motif[j,],collapse=" ")
                                                             vec <- c(vec,out_vec)}
                                                         write(vec, outname, append=TRUE)}
                                                             
}

FIMO_to_granges <- function(FIMO_bed){
    bed <- read.table(FIMO_bed, header=T, stringsAsFactors=F,sep='\t')
    gr <- GRanges(seqnames=Rle(bed[,1]), ranges= IRanges(start=bed[,2], end=bed[,3]),strand=bed[,4], PWM_Score=bed[,5],P_value=bed[,6], sequence=bed[,7], motif=bed[,8])
                                 return(gr)}

Collapse_gr_keep_Gene <- function(gr){
    reduced <- reduce(gr)
    final_gr <- reduced
    for (i in 1:length(reduced)){
        overlap_enhancer = findOverlaps(reduced[i], gr)
        row_names = vector()
        for( j in 1:length(overlap_enhancer)){
            names = gr[overlap_enhancer[j]@to]$Gene
            row_names = c(names, row_names)
            
        }

        final_gr$nearby[i] = row_names
    }
    return(final_gr)}


qqplot <- function(pval_vector,distribution){
    if(distribution %in% c("unif","pois","norm")){
        print("Working with provided distribution")
        pval_df <- data.frame(sort(-log10(pval_vector)),stringsAsFactors=F)
        colnames(pval_df) <- "observed"

        pval_df$expected <- sort(-log10(runif(length(pval_vector))))}
        str(pval_df$expected)
        p <- ggplot(pval_df,aes(expected,observed))+geom_point(size=2)+geom_abline(slope=1,intercept=0,size=2,color="red",linetype=2, alpha=0.5)+labs(title=paste0("QQplot with ",distribution," distribution"))
        return(p)}

lseq <- function(from=1, to=100000, length.out=19){ exp(seq(log(from), log(to), length.out=length.out))}

get_nearest_gr <- function(gr1,gr2){
    int <- distanceToNearest(gr1,gr2)
    sub <- gr1[int@from]
    sub$Distance <- int@elementMetadata@listData$distance
    sub$closest <- int@to
    return(sub)}


pintersect_with_metadata <- function(gr1,gr2){
    hits <- findOverlaps(gr,gr2)
    gr.over <- pintersect(gr[queryHits(hits)],gr2[subjectHits(hits)])
    gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x) sum(width(x)))
    gr$overlap<- 0
    gr$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
    return(gr)}

empty_gr <- function(){
    gr <- GRanges(seqnames='chr1',ranges=IRanges(start=1,end=2))
    return(gr)}

vector_smartmerge <- function(df1,df2,margin){
    if(margin == 'row'){
        row_names <- union(rownames(df1),rownames(df2))
        df_int <- cbind(df1,0)
        col_index <- dim(df_int)[2]
        df_int[rownames(df2),col_index] <- df2[,1]
        colnames(df_int)[col_index] <- colnames(df2)
    }
    else{print("Haven't worked this out yet")}
             return(df_int)}   

run_tsne <- function(df,sample_margin="column",perplexity=NULL){
    if(sample_margin == "column"){ df <- as.data.frame(t(df),stringsAsFactors=F)
                                   rownames(df) <- gsub("\\.","-",rownames(df))
                               }
    else if (sample_margin == "row") { df <- df}
    if(is.null(perplexity) == TRUE){ hyperparam = round((nrow(df)-1)/3)-1}
    else{ hyperparam=perplexity}
    out <- Rtsne(df,check_duplicates=F,perplexity=hyperparam)$Y
    out <- as.data.frame(out,stringsAsFactors=F)
    out$Sample <- rownames(df)
    rownames(out) <- out[,3]
    colnames(out) <- c("Dim1","Dim2","Sample")
    return(out)
}

get_metadata_tsne <- function(df, metadata_file,column=2){
    df$Subtype <- metadata_file[rownames(df),column]
    return(df)}
    
plot_tsne <- function(df,title,color="Subtype",label="Subtype"){
    p <- ggplot(df, aes_string("Dim1","Dim2", color=color,label=label))+geom_point(size=3, alpha=0.5)+geom_text(check_overlap=TRUE, nudge_y=0.5,size=5,show.legend=F)+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)
    return(p)
}

get_drivers_from_clusters <- function(df,tsne_df){
    tsne_df$Cluster_drivers <- "None"
    tsne_df$Top_Driver <- "None"
    colnames(df) <- gsub("\\.","-",colnames(df))
    pval_all <- data.frame(stringsAsFactors=F)
    print("Initializing")

    for(i in unique(tsne_df$Cluster)){
        sample_list <- tsne_df[which(tsne_df$Cluster == i),"Sample"]
        mat_sub <- df[sample_list]
        other_mat <- df[setdiff(colnames(df),sample_list)]
        pval_df <- data.frame(stringsAsFactors=F)
 #       str(sample_list)

        for(j in rownames(other_mat)){
            pval <- wilcox.test(unlist(mat_sub[j,]),unlist(other_mat[j,]))$p.value

            cluster_avg <- mean(unlist(mat_sub[j,]))
            alt_avg <- mean(unlist(other_mat[j,]))
            sample_num <- length(unlist(mat_sub[j,]))
            cluster_id <- i
         
            sign <- mean(unlist(mat_sub[j,])) > mean(unlist(other_mat[j,]))
            pval_df <- rbind(pval_df,cbind(cluster_avg,alt_avg,sample_num,cluster_id,pval,j,sign),stringsAsFactors=F)}
        pval_df$qval <- p.adjust(pval_df$pval,method="BH")

        print(table(pval_df$sign))
#        pval_df <- pval_df[which(pval_df$sign =="TRUE"),]
        pval_df <- pval_df[order(pval_df$qval),]
        top10 <- pval_df[1:10,"j"]
        top_driver <- top10[1]
        #str(top_driver)
        #print(top10)
        tsne_df[sample_list,"Cluster_drivers"] <- paste0(top10,collapse="-")
        tsne_df[sample_list,"Top_Driver"] <- top_driver

        pval_df$Cluster <- paste0("Cluster-",i)
        pval_all <- rbind(pval_all,pval_df,stringsAsFactors=F)
       # str(pval_all)
    }
    return(list(pval_all,tsne_df))
}

make_clusters_from_metric_list <- function(metric_list,cutoff_quantile=0.95,metadata_file){
    metric_clusters <- list()
    for(i in 1:length(metric_list)){
        print(paste0("Working ",names(metric_list)[i]))
        df <- metric_list[[i]]
        tsne_out <- run_tsne(df)
        print("tSNE done")
        clusters_tsne <- hclust(dist(tsne_out[,c(1,2)]))
        members <- cutree(clusters_tsne, h=quantile(clusters_tsne$height,cutoff_quantile))
        print("Clusters identified")
        cluster_out <- as.data.frame(cbind(members,names(members)),stringsAsFactors=F)
        tsne_out <- merge(tsne_out,cluster_out,by.x="Sample",by.y="V2")
        print("Merging clusters to tSNE")
        rownames(tsne_out) <- tsne_out$Sample
        tsne_out <- get_metadata_tsne(tsne_out, metadata_file)
        colnames(tsne_out)[4] <- "Cluster"
        #print(head(tsne_out))
        tsne_out <- get_drivers_from_clusters(df,tsne_out)
        #str(tsne_out)
        metric_clusters <- c(metric_clusters, list(tsne_out)); print(paste0("Finished ",i))}
    names(metric_clusters) <- names(metric_list)
    return(metric_clusters)
}





make_TF_target_mat <- function(TF,indir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/patientEdgeList/",outdir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/"){
    source("/home/anf2034/Andre_F_functions.R")
    target_list <- list()
    print(paste0("Working on ", TF))

file_list <- list.files(indir)
file_list <- file_list[grep("_graph",file_list)]
#file_list <- file_list[1:10]

for(i in file_list){graph <- readRDS(paste0(indir,i))
                    targets <- get_regulatory_targets(graph,TF)
                    target_list <- c(target_list,list(targets))
                    print(paste0("Finished ", grep(i,file_list)," of ",length(file_list)))

                }
saveRDS(target_list, paste0(outdir,TF,"_target_list.rds"))


names(target_list) <- gsub("_graph.rds","", file_list)

union_targets <- c();  for(i in target_list){ union_targets <- unique(c(union_targets,i))}
TF_matrix <- matrix(0,length(target_list),length(union_targets))
rownames(TF_matrix) <- names(target_list)
colnames(TF_matrix) <- union_targets
str(TF_matrix)
for(i in 1:length(target_list)){ targets <- target_list[[i]]; TF_matrix[i, targets] <- 1}

target_summary <- colSums(TF_matrix)/nrow(TF_matrix)

print(paste0("Variance for ", TF,"= ",round(var(target_summary),digits=2)))

saveRDS(TF_matrix,paste0(outdir,TF,"_target_matrix.rds"))
saveRDS(target_summary,paste0(outdir,TF,"_target_summary.rds"))


}

make_var_estimate_test <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/var_estimate.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")
    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),1]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        str(in_mat)
        sub_mat <- in_mat[in_samples,]
        inter_summary <- round(colSums(in_mat)/nrow(in_mat),digits=2)
        intra_summary <- round(colSums(sub_mat)/nrow(sub_mat),digits=2)

        inter_summary <- 100*table(inter_summary)/sum(table(inter_summary))
        intra_summary <- 100*table(intra_summary)/sum(table(intra_summary))

        subtype_df <- rbind(subtype_df,cbind(inter_summary["1"],intra_summary["1"],sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_variability.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_variability.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}
        
make_var_estimate <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/TF_list_Erica_networks.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")
    
    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),1]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]
        inter_summary <- round(colSums(in_mat)/nrow(in_mat),digits=2)
        intra_summary <- round(colSums(sub_mat)/nrow(sub_mat),digits=2)

        inter_summary <- 100*table(inter_summary)/sum(table(inter_summary))
        intra_summary <- 100*table(intra_summary)/sum(table(intra_summary))

        subtype_df <- rbind(subtype_df,cbind(inter_summary["1"],intra_summary["1"],sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_variability.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_variability.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}





make_dist_estimate <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/TF_list_Erica_networks.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")
    
    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),1]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]

        
        intra_dist <- sd(dist(sub_mat))/mean(dist(sub_mat))
        inter_dist <- sd(dist(in_mat))/mean(dist(in_mat))


        subtype_df <- rbind(subtype_df,cbind(inter_dist,intra_dist,sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}

make_dist_estimate_test <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/var_estimate.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")
    
    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),1]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]

        intra_dist <- sd(dist(sub_mat))/mean(dist(sub_mat))
        inter_dist <- sd(dist(in_mat))/mean(dist(in_mat))
        

        subtype_df <- rbind(subtype_df,cbind(inter_dist,intra_dist,sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}

gr.mid <- function (x) {
    start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
    return(x)
}

make_graph_from_el <- function(infile, outdir){
    el <- as.matrix(read.table(infile, header=F, sep='\t', stringsAsFactors=F))
    filename <- paste0("TCGA",unlist(strsplit(infile,"TCGA"))[2])
    graph <- graph_from_edgelist(el)
    saveRDS(graph, paste0(outdir,gsub("edgelist.txt","graph.rds",filename)))}

    
rewiring_auc <- function(dir){
    file_list <- list.files(dir)
    file_list <- file_list[grep("_matrix.rds", file_list)]
    out_list <- list()
    for(i in file_list){
        in_mat <- readRDS(paste0(dir,i))
        mat_index_order <- names(sort(rowSums(in_mat)))
        target_all_len <- c()
        target_len <- c()
        for(j in 1:length(mat_index_order)){
            target_vec <- names(which(in_mat[mat_index_order[j],] == 1))
            target_len <- union(target_len,target_vec)
            target_all_len <- c(target_all_len,length(target_len))
            
                      
        }
        target_all_len <- as.data.frame(target_all_len,stringsAsFactors=F)
        colnames(target_all_len) <- "Union_Targets"
            
        target_all_len$Num_Samples <- 1:dim(target_all_len)[1]
        
        out_list <- c(out_list,list(target_all_len))
        print(paste0("Finished working ", unlist(strsplit(i,"_"))[1]))}
    names(out_list) <- gsub("_target_matrix.rds","",file_list)
    return(out_list)}
                                          
add_graphs <- function(graph1,graph2){
    if((is.weighted(graph1) & is.weighted(graph2)) == FALSE){
        print("One of the graphs is not weighted, please recheck")} else{
            graph_out <- graph1 + graph2
            E(graph_out)$weight_1[is.na(E(graph_out)$weight_1)] <- 0
            E(graph_out)$weight_2[is.na(E(graph_out)$weight_2)] <- 0
            E(graph_out)$weight <- E(graph_out)$weight_1+E(graph_out)$weight_2
            return(graph_out)}}

subtract_graphs <- function(graph1, graph2){
    if((is.weighted(graph1) & is.weighted(graph2)) == FALSE){
        print("One of the graphs is not weighted, please recheck")} else{
            graph_out <- union(graph1,graph2)
            E(graph_out)$weight_1[is.na(E(graph_out)$weight_1)] <- 0
            E(graph_out)$weight_2[is.na(E(graph_out)$weight_2)] <- 0
            E(graph_out)$weight <- E(graph_out)$weight_1 - E(graph_out)$weight_2
            return(graph_out)}}
              
read_peak_gene_el <- function(edgelist,cutoff=0.5){


    
    edgelist <- read.table(edgelist, sep=',',stringsAsFactors=F,header=F)
    split1 <- unlist(strsplit(edgelist[,1],"-(?=[^-]+$)", perl=TRUE))
    edgelist$Peak <- split1[seq(2,length(split1),2)]
    edgelist$Gene <- split1[seq(1,length(split1),2)]
    edgelist <- edgelist[which(edgelist[,2] >cutoff),]
    edgelist <- edgelist[,c(4,5,2)]
    colnames(edgelist) <- c("Peak","Gene","Weight")

    return(edgelist)}

weighted_graph_from_edgelist <- function(edgelist, column=3){
    graph_out <- graph_from_edgelist(as.matrix(edgelist[c(1,2)]))
    E(graph_out)$weight <- 1
    return(graph_out)

}

normalize_peak_gene_matrix <- function(matrix,margin='row'){
    if(margin == "row"){
        sample_count <- apply(matrix,1, max)
        matrix_out <- matrix/sample_count
    }
    else if(margin=="column"){
        sample_count <- apply(matrix,2, max)
        matrix_out <- matrix/sample_count
    }
    return(matrix_out)}
        
        

#object_size <- function(objects){
#    if(is.vector(objects) == TRUE){
        
    

#psetdiff_with_metadata <- function(gr1,gr2){
 #   hits <- findOverlaps}


expand.matrix <- function(A, sparse=TRUE){
    if(sparse ==FALSE){
        m <- nrow(A)
        n <- ncol(A)
        B <- matrix(0,nrow = m, ncol = m)
        C <- matrix(0,nrow = n, ncol = n)
        out <- cbind(rbind(B,t(A)),rbind(A,C))} else{
            m <- nrow(A)
            n <- ncol(A)
            B <- Matrix(0,nrow = m, ncol = m,sparse=TRUE)
            C <- Matrix(0,nrow = n, ncol = n,sparse=TRUE)
            out <- cbind(rbind(B,t(A)),rbind(A,C))}
    return(out)
}

expand.graph.matrix <- function(A,sparse=TRUE){
    if(sparse ==TRUE){
        print("Working 1")
        m <- nrow(A)
        n <- ncol(A)
        B <- Matrix(0,nrow = m, ncol = m,sparse=TRUE)
        C <- Matrix(0,nrow = n, ncol = n,sparse=TRUE)
        int1 <- rbind(B,t(A))
        int2 <- rbind(A,C)
        out <- cbind(int1,int2)} else{
            print("Working 2")
             m <- nrow(A)
        n <- ncol(A)
        B <- matrix(0,nrow = m, ncol = m)
             C <- matrix(0,nrow = n, ncol = n)
             

             out <- cbind(rbind(B,t(A)),rbind(A,C))}


    out<- triu(out)
    return(out)
    
}

get_non_zero <- function(mat,margin="row"){
    out_vec <- c()
    for(i in 1:nrow(mat)){
        out_vec <- c(out_vec,all(mat[i,] == 0))
    }
    return(out_vec)

}

get_most_variable <- function(mat, margin='row',cutoff=0,quantile=NULL){
    if( margin =="row"){
        vec <- apply(mat, 1, var)

        variable <- which(vec >0)

        #out <- mat[variable,]
    }
    return(variable)}

sparse_mat_to_vector <- function(mat){
    mat <- as(mat, "dgTMatrix")
    out <- mat@x
    names(out) <- paste0(rownames(mat)[(mat@i+1)],"~",colnames(mat)[(mat@j+1)])

    return(out)}

cutdown.graph.matrix <- function(graph_mat){
    graph_sub_mat <- graph_mat[grep("[A-Z]_[0-9]", rownames(graph_mat)),]
    graph_sub_mat <- graph_sub_mat[,setdiff(colnames(graph_sub_mat),rownames(graph_sub_mat))]
    return(graph_sub_mat)}

    
make_combined_metric <- function(outdegree_df,expression_df){
    common_genes <- intersect(rownames(outdegree_df),rownames(expression_df))
    colname_index <- intersect(colnames(outdegree_df),colnames(expression_df))

    outdegree_df <- outdegree_df[common_genes,colname_index]
    expression_df <- expression_df[common_genes,colname_index]

    out_df <- (outdegree_df+expression_df)/2

    return(out_df)}
    

    
get_tissue_specific_peaks <- function(abundance_mat,zscore=2.5){
    sub <- apply(abundance_mat,1,function(x) any(((unlist(x)-mean(unlist(x)))/sd(unlist(x))) >=zscore))
    sub <- names(which(sub == TRUE))
    print(paste0(length(sub)," links survived first cut"))
    sub_mat <- abundance_mat[sub,]
    index <- apply(sub_mat,1, function(x) any(x >0.5))
    index <- index[which(index == TRUE)]
    print(paste0(length(index)," links survived second cut"))
    sub_mat <- sub_mat[names(index),]
    return(sub_mat)
}


get_coordinate <- function(peak,gr,peak_column){
    metadata <- unlist(gr@elementMetadata@listData[peak_column])
    index <- grep(peak,metadata)
    chrom <- as.character(gr@seqnames)[index]
    start <- gr@ranges@start[index]
    end <- start+gr@ranges@width[index]
    Peaks <- peak
    coord_df <- cbind(chrom,start,end,Peaks)
    return(coord_df)}
    
PWM_to_Jaspar_norm <- function(PWM){
    PWM_name <- names(PWM@name)

    PWM_mat <- PWM@profileMatrix
    PWM_out <- list(paste0(">",PWM_name))
    for(i in 1:nrow(PWM_mat)){ vec <- paste0(rownames(PWM_mat)[i]," \t"," [ ", paste0(unlist(PWM_mat[i,]),collapse=" ")," ] ")
                               PWM_out <- c(PWM_out,list(vec))}
    return(PWM_out)
}

PIQ_bed_to_gr <- function(bed){
    check <- readLines(bed, n=2)
    search_check <- grep("track", check)
    if(length(search_check) >0){
        new_bed <- gsub(".bed", "_fixed.bed",bed);
        system(paste0("tail -n +2 ",bed, " > ",new_bed))
        bed <- new_bed; print("Bed file was incorrectly formatted")} else{ print("Bed file was correctly formatted")}
     bed <- read.table(bed, stringsAsFactors=F, sep="\t",header=F)
                                 gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),purity=bed$V5,strand=bed$V6)
                                 return(gr)}
 
