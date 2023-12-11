plot_hybrid_network_modules <- function(hybrid_net,title=NULL,module1="ensembl_id_1_module",module2="ensembl_id_2_module",drop_low_freq=FALSE,drop0=FALSE,highlight_gene=NULL,highlight_module=NULL,net_layout="seham"){
    require(igraph)
    require(ggnetwork)
    require(tidyr)
    require(ggplot2)
    require(tempusr)
    require(scales)



### Transform hybrid network so that undirected edges are converted to bidirectional edges for the sake of visualization
    sub <- dplyr::filter(hybrid_net, directed==TRUE)
    int <- dplyr::filter(hybrid_net, directed==FALSE)
    sub2 <- int
    genelist1 <- int$ensembl_id_1
    genelist2 <- int$ensembl_id_2
    int$ensembl_id_1 <- genelist2
    int$ensembl_id_2 <- genelist1
    sub2 <- rbind(sub2, int)
    hybrid_net <- rbind(sub,sub2)

### Get HGNC symbols
    gene_df <- tempusr::gene_annotation; rownames(gene_df) <- gene_df[,1]

    hybrid_net$Gene1 <- gene_df[hybrid_net[,"ensembl_id_1"],"hgnc_symbol"]
    hybrid_net$Gene2 <- gene_df[hybrid_net[,"ensembl_id_2"],"hgnc_symbol"]


### Remove genes in Module 0/Grey if requested
    if(drop0){ hybrid_net <- dplyr::filter(hybrid_net, ensembl_id_1_module!=0,ensembl_id_2_module!=0)} else{ print("Not removing Module 0/Gray module")}

### highlight gene workflow

        if(!is.null(highlight_gene) & is.null(highlight_module)){
            print("Highlighting gene")
            hybrid_net$Edge_Color <- "black"
            if(all(grepl("ENSG[0-9]{11}",highlight_gene))){ gene_id_col <- c("ensembl_id_1","ensembl_id_2")} else{ gene_id_col <- c("Gene1","Gene2")}
            print(gene_id_col)
            index <- lapply(gene_id_col, function(x) which(hybrid_net[,x] %in% highlight_gene))

            if(length(unlist(index)) ==0 & drop0 ){
                print("Couldn't find your gene of interest, may not be in network and/or assigned to Module 0. Not highlighting any genes")
            } else if(length(unlist(index))==0 & !drop0) {
                print("Couldn't find your gene of interest in the network. Not highlighting any genes")} else{
                                                                                                           hybrid_net$Edge_Color[index[[1]]] <- "red"
                                                                                                           hybrid_net$Edge_Color[index[[2]]] <- "green"
                                                                                                       }} else if (is.null(highlight_gene) & !is.null(highlight_module)){
            print("Highlighting module")
            hybrid_net$Edge_Color <- "black"

            index <- lapply(c(module1,module2), function(x) which(hybrid_net[,x] %in% highlight_module))
            if(length(unlist(index)) ==0 & drop0 ){
                print("Couldn't find your gene of interest, may not be in network and/or assigned to Module 0. Not highlighting any genes")
            } else if(length(unlist(index))==0 & !drop0) {
                print("Couldn't find your gene of interest in the network. Not highlighting any genes")} else{
                                                                                                           hybrid_net$Edge_Color[index[[1]]] <- "red"
                                                                                                           hybrid_net$Edge_Color[index[[2]]] <- "green"}
                                                                                                        } else{ print("Not highlighting anything")
                                                                                                            hybrid_net$Edge_Color <- "black"}


    hybrid_net_summary <- lapply(unique(hybrid_net$Edge_Color), function(x) as.data.frame(table(dplyr::filter(hybrid_net, Edge_Color==x)[c("ensembl_id_1_module","ensembl_id_2_module")])))
    for(i in 1:length(hybrid_net_summary)){ hybrid_net_summary[[i]]$Edge_Color <- unique(hybrid_net$Edge_Color)[i]}
    hybrid_net_summary <- do.call("rbind", hybrid_net_summary)



    hybrid_net_summary[1:2] <- apply(hybrid_net_summary[1:2],2, function(x) as.integer(as.character(x)))
    hybrid_net_summary[,4] <- as.character(hybrid_net_summary[,4])

    colnames(hybrid_net_summary) <- c("Module1","Module2","Edge_Count","Edge_Color")

    hybrid_count <- dplyr::select(hybrid_net, Module1=ensembl_id_1_module,Gene=ensembl_id_1) %>% unique %>% group_by(Module1) %>% dplyr::summarize(Module_Node_Count=n())
    if(drop_low_freq){
        cutoff_val <- sum(hybrid_count)*0.01
        hybrid_net_summary <- dplyr::filter(hybrid_net_summary, Edge_Count>=cutoff_val)   }


    hybrid_net_summary <- left_join(hybrid_net_summary, hybrid_count, by="Module1") %>% group_by(Module1) %>% mutate(Edge_Weight= Edge_Count/Module_Node_Count) %>% ungroup

    hybrid_net_summary$Area <- rescale(hybrid_net_summary$Module_Node_Count,to =c(2,20))

    hybrid_net_summary[1:2] <- apply(hybrid_net_summary[1:2],2, function(x) as.character(x))


    if(!is.null(highlight_gene)){

        print("Working gene")

        net_list <- lapply(unique(hybrid_net_summary$Edge_Color), function(x) dplyr::filter(hybrid_net_summary, Edge_Color==x))
##        for(i in 1:length(net_list)){  sub <- net_list[[i]]; sub$Edge_ID <- 1:nrow(sub); net_list[[i]] <- sub}

        net_list[[1]]$Edge_ID <- 1:nrow(net_list[[1]])
        net_out <- network::as.network(unique(dplyr::filter(net_list[[1]], Module1!=Module2,Edge_Count!=0,!is.na(Edge_Count))),multiple=FALSE)
        net_out <- ggnetwork(net_out, layout = net_layout, arrow.gap=0.03)

        for(i in 2:length(net_list)){  sub <- net_list[[i]]; alt <- net_list[[1]]

            int <- left_join(dplyr::select(sub,Module1,Module2,Edge_Count,Edge_Color, Module_Node_Count, Area),dplyr::select(alt, Module1,Module2, Edge_ID),by=c("Module1","Module2"))

            int2 <- dplyr::filter(net_out, Edge_ID %in% int$Edge_ID)

            int2 <- merge(dplyr::select(int2,-Edge_Count),dplyr::select(int,Edge_ID,Edge_Count), by="Edge_ID")


            int2$Edge_Color <- unique(sub$Edge_Color)
            net_out <- rbind(net_out, int2) }

        net_out <- dplyr::filter(net_out, !is.na(Edge_Count))

    } else{
        print("Working other")
    net_out <- network::as.network(unique(dplyr::filter(hybrid_net_summary, Module1!=Module2,Edge_Count!=0,!is.na(Edge_Count))),multiple=FALSE)
    net_out <- ggnetwork(net_out, layout = net_layout, arrow.gap=0.03)

    net_out <- dplyr::filter(net_out, !is.na(Edge_Count))
##
    }
    saveRDS(net_out,"net_out.rds")

    black <- unique(dplyr::filter(net_out, Edge_Color=="black"))
    red <- unique(dplyr::filter(net_out, Edge_Color=="red"))
    green <- unique(dplyr::filter(net_out, Edge_Color=="green"))

    ggplot(net_out, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(data=black,aes(alpha=Edge_Weight,lwd=Edge_Count),color=black$Edge_Color,curvature=0.1,angle=10,arrow=arrow(length=unit(4, "pt"), type = "closed",angle=25)) +
        geom_edges(data=green,aes(lwd=50*Edge_Count),color="green",curvature=0.4,angle=90,lty=4,arrow=arrow(length=unit(4, "pt"), type = "closed",angle=25),show.legend=F) +
        geom_edges(data=red,aes(lwd=50*Edge_Count),color="red",curvature=0.7,angle=90,lty=2,arrow=arrow(length=unit(4, "pt"), type = "closed",angle=25),show.legend=F) +
        geom_nodes(aes(color=vertex.names,size=Area),alpha=0.8,show.legend=F) +
        geom_nodes(aes(size=Area*2,color=vertex.names),alpha=0.3,show.legend=F)+
        theme_void()+
        geom_nodelabel(aes(label=gsub("Module_","",vertex.names),size=0.1*Area),show.legend=F,fontface=2,color="gray50",label.padding=unit(0.1,"lines"))+
        scale_linewidth_binned(range=c(0.001,3),n.breaks=6)+
        scale_alpha_binned(range=c(0.05,0.7),n.breaks=5)+
        scale_size_area("Module_Node_Count", n.breaks = 10,max_size = 20)+
        labs(title=title)+theme(plot.title=element_text(hjust=0.5))
}

network_colorscale <- function(color_list=NULL,extension=".txt",scales=c("fill","color"),...){
    require(ggplot2)

    if(length(scales)>1){ scales <- "color"}
    if(is.null(color_list)==TRUE){
        cols <- data.frame(c(0:94,paste0("ME",0:94),paste0("Module_",0:94)),c("grey","turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan","midnightblue","lightcyan","grey60","lightgreen","lightyellow","royalblue","darkred","darkgreen","darkturquoise","darkgrey","orange","darkorange","white","skyblue","saddlebrown","steelblue","paleturquoise","violet","darkolivegreen","darkmagenta","sienna3","yellowgreen","skyblue3","plum1","orangered4","mediumpurple3","lightsteelblue1","lightcyan1","ivory","floralwhite","darkorange2","brown4","bisque4","darkslateblue","plum2","thistle2","thistle1","salmon4","palevioletred3","navajowhite2","maroon","lightpink4","lavenderblush3","honeydew1","darkseagreen4","coral1","antiquewhite4","coral2","mediumorchid","skyblue2","yellow4","skyblue1","plum","orangered3","mediumpurple2","lightsteelblue","lightcoral","indianred4","firebrick4","darkolivegreen4","brown2","blue2","darkviolet","plum3","thistle3","thistle","salmon2","palevioletred2","navajowhite1","magenta4","lightpink3","lavenderblush2","honeydew","darkseagreen3","coral","antiquewhite2","coral3","mediumpurple4","skyblue4","yellow3"))
        colnames(cols) <- c("Module","Color")
##        str(cols)




    } else if(is.character(color_list) & !is.vector(color_list)) { if(extension == ".rds" & !is.data.frame(cols)) { print("reading RDS"); cols <- readRDS(color_list)} else if(extension==".txt" & !is.data.frame(cols)){ print("reading txt");cols <- read.table(color_list,sep='\t',stringsAsFactors=F,header=T,comment.char="$")}}

        if(is.data.frame(cols)){ df <- cols; cols <- df[,2]; names(cols) <- df[,1]} else if(is.vector(color_list)==TRUE){ cols <- color_list}


    g <- ggplot2:::manual_scale(scales, values=cols)
    return(g)
}
