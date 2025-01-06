compare_hybrid_network_modules <- function(network1, network2,method=c("jaccard","chisquared"),highlight_gene=NULL,translate_ensembl=NULL){

    require(tempusr)
    require(dplyr)



    network1 <- unique(rbind(dplyr::select(network1, Gene=ensembl_id_1,Module=ensembl_id_1_module),dplyr::select(network1, Gene=ensembl_id_2,Module=ensembl_id_2_module)))
    network2 <- unique(rbind(dplyr::select(network2, Gene=ensembl_id_1,Module=ensembl_id_1_module),dplyr::select(network2, Gene=ensembl_id_2,Module=ensembl_id_2_module)))


    out <- data.frame(stringsAsFactors = F)
    for(i in 1:nrow(network1)){
        gene <- network1[i,1]
        module1 <- network1[i,2]
        module2 <- unique(dplyr::filter(network2,Gene==gene)[,2])
        if(length(module1) >0 && length(module2)>0){

        vec1 <- dplyr::filter(network1,Module==module1)[,1]
        vec2 <- dplyr::filter(network2,Module==module2)[,1]

        if(method=="jaccard"){

            similarity <- length(intersect(vec1,vec2))/length(union(vec1,vec2))
            if(length(vec1)>length(vec2)){ smallest <- vec2; largest <- vec1} else{ smallest <- vec1; largest <- vec2}
            nestedness <- length(intersect(smallest,largest))/length(smallest)

            int <- data.frame(gene, module1,module2,similarity,nestedness)

            colnames(int) <- c("Gene","Module1","Module2","Jaccard_Similarity","Nestedness")

        } else if (method=="chisquared"){
            TP <- length(intersect(vec1,vec2))
            FP <- length(setdiff(vec1,vec2))
            FN <- length(setdiff(vec1,x=vec2))
            TN <- length(union(vec1,vec2))
            mat <- matrix(c(TP,FP,FN,TN),nrow=2,ncol=2)
            mat_result <- chisq.test(mat)
            int <- data.frame(gene,module1,module2,mat_result$statistic,mat_result$p.value)
            colnames(int) <- c("Gene","Module1","Module2","X-squared","P_value")



        }
            out <- rbind(out,int)}}
    if(!is.null(translate_ensembl) &&  translate_ensembl==TRUE){
            gene_df <- tempusr::gene_annotation
            rownames(gene_df) <- gene_df[, 1]
            out$Gene <- gene_df[out$Gene,"hgnc_symbol"]
} else{ out <- out}

    return(out)}



compare_wgcn_modules <- function(gene_info1, gene_info2,method=c("jaccard","chisquared"),highlight_gene=NULL,translate_ensembl=NULL){

    require(tempusr)
    require(dplyr)
    if(length(method)>=2){ method <- "jaccard"} else{ method <- method}





    network1 <- dplyr::select(gene_info1,Gene_ID=ensembl_gene_id,Module_Label=module_label,Gene=hgnc_symbol,Module_Color=module_color)
    network2 <- dplyr::select(gene_info2,Gene_ID=ensembl_gene_id,Module_Label=module_label,Gene=hgnc_symbol,Module_Color=module_color)


    out <- data.frame(stringsAsFactors = F)
    for(i in sort(unique(network1$Module_Label))){

        module1 <- i
        vec1 <- dplyr::filter(network1,Module_Label==module1)$Gene_ID
        int <- data.frame(stringsAsFactors = F)
        for(j in sort(unique(network2$Module_Label))){
            module2 =j
            vec2 <- dplyr::filter(network2,Module_Label==module2)[,1]
            if(method=="jaccard"){
                similarity <- length(intersect(vec1,vec2))/length(union(vec1,vec2))
                if(length(vec1)>length(vec2)){ smallest <- vec2; largest <- vec1} else{ smallest <- vec1; largest <- vec2}
                nestedness <- length(intersect(smallest,largest))/length(smallest)
                sub <- data.frame(module1,module2,similarity,nestedness,length(vec1), length(vec2),length(intersect(vec1,vec2)),length(union(vec1,vec2)))
                int <- rbind(int,sub)

            } else if (method=="chisquared"){

                TP <- length(intersect(vec1,vec2))
                FP <- length(setdiff(vec1,vec2))
                FN <- length(setdiff(vec1,x=vec2))
                TN <- length(union(vec1,vec2))
                mat <- matrix(c(TP,FP,FN,TN),nrow=2,ncol=2)
                mat_result <- chisq.test(mat)
                sub <- data.frame(module1,module2,mat_result$statistic,mat_result$p.value)
                int <- rbind(int,sub)

        }}
            out <- rbind(out,int)}
    if(!is.null(translate_ensembl) &&  translate_ensembl==TRUE){
            gene_df <- tempusr::gene_annotation
            rownames(gene_df) <- gene_df[, 1]
            out$Gene <- gene_df[out$Gene,"hgnc_symbol"]
} else{ out <- out}

    if(method=="jaccard"){colnames(out) <- c("Module1","Module2","Jaccard_Similarity","Nestedness","Module1_Size","Module2_Size","Intersection","Union")
        out <- out %>% group_by(Module1,Module2) %>% mutate(Dist_to_origin=dist(rbind(c(0,0,0),c(Jaccard_Similarity,Nestedness,Intersection)))) %>% ungroup %>% group_by(Module1) %>% mutate(Best_Match=Dist_to_origin==max(Dist_to_origin)) %>% data.frame
    } else if(method=="chisquared"){colnames(out) <- c("Module1","Module2","X-squared","P_value")}
            return(out)}

expand_hybrid_network <- function(hybrid_net){
    require(dplyr)
  sub <- dplyr::filter(hybrid_net, .data$directed == TRUE)
  int <- dplyr::filter(hybrid_net, .data$directed == FALSE)
  sub2 <- int
  genelist1 <- int$ensembl_id_1
  genelist2 <- int$ensembl_id_2
  int$ensembl_id_1 <- genelist2
  int$ensembl_id_2 <- genelist1
  sub2 <- rbind(sub2, int)
  hybrid_net <- rbind(sub, sub2)
  hybrid_net$Edge <- paste0(hybrid_net[,"ensembl_id_1"],"~",hybrid_net[,"ensembl_id_2"])


  return(hybrid_net)
}


get_furthest <- function(coordinate1, coordinate_df){
    distances <- unlist(lapply(1:nrow(coordinate_df), function(x) dist(rbind(coordinate1, coordinate_df[x,]))))

#    str(distances)
    max_dist <- which(distances==max(distances))
    return(max_dist)}


## Get best match for each cluster in network 1. This uses the output from compare wgcn_modules or
## wgcn_comparison <- wgcn_comparison %>% group_by(Module1) %>% mutate(Dist_to_origin=dist(rbind(c(0,0,0),c(Jaccard_Similarity,Nestedness,Intersection))),Best_Match=Dist_to_origin==max(Dist_to_origin)) %>% data.frame


## ##Plotting code ##

## ggplot(test_wgcn_comparison, aes(Jaccard_Similarity,Nestedness,fill=Best_Match))+ geom_point(aes(size=Union),color="gray50",alpha=0.4,shape=21)+theme_anf+geom_text_repel(aes(label=paste0(Module1,"~",Module2)))
