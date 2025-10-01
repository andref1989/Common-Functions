#' Plot hybrid networks in a module focused fashion, inspired by figures from Lauren Beck PhD
#'
#' @param hybrid_net Data frame of the hybrid network outputs p0069
#' @param title Plot title (Optional string)
#' @param module1 The name of the column in the data frame containing the module
#' associated with gene1 in the network (Character) Default: ensembl_id_1_module
#' @param module2 The name of the column in the data frame containing the module
#' associated with gene2 in the network (Character) Default: ensembl_id_2_module
#' @param drop_low_freq Whether to remove connections between modules that are
#' infrequently interacting in the network (Boolean) Default: FALSE
#' @param drop0 Whether to remove connections module 0 in the network. Module 0
#' typically contains genes not assigned to other more closely related/interacting
#' and/or pathway-associated modules (Boolean) Default: FALSE
#' @param highlight_gene A gene or vector of genes whose connections you want to
#' highlight. (Character) Default: NULL
#' @param highlight_module A module or vector of modules whose connections you
#' want to highlight (Character) Default: NULL
#' @param highlight_type What kinds of edges should be used for the highlighting of genes/modules
#' @param net_layout The layout to use for the network. Values allowed: circrand,
#' eigen, fruchtermanreingold, geodist ,hall, kamadakawai, mds,princoord, random,
#' rmds, segeo, seham, spring, springrepulse or target (Character) Default: seham
#' @param out_edge_col Color of edges leaving your gene/module(s) of interest (Character) Default: Red
#' @param in_edge_col Color of edges entering your gene/module(s) of interest (Character) Default: Green
#' @param gene_annotation Data frame containing a minimum of ensembl_gene_id and hgnc_symbol
#'
#'
#' @return Plot of your network!
#' @note Ordinary connections between modules are colored black with a linewidth
#' determined by the number of genes interacting between those modules.
#' @note If highlighting a gene or module, red edges are those leaving your
#' highlighted module and green edges are those entering your highlighted module.
#' By default, undirected edges are considered bidirection whereas directed edges
#' are unidirectional.
#' @note Currently edge weight is being calculated as the number of edges leaving
#'  a module A and entering module(s) B-to-Z divided by the number of nodes in
#'  that module. It is essentially a measure of the connectedness level of a
#'  given pair of modules.  This can and likely will change going forward.
#' @note If desired we can provide an option to alter the colors used for
#' gene/module highlighting.
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#'
#' plot_hybrid_network_modules(hybrid_net,
#'   title = "Hybrid network",
#'   drop0 = FALSE,
#'   drop_low_freq = FALSE,
#'   highlight_module = "TP53",
#'   net_layout = "seham")
#'   + network_colorscale()
#' }
#'
plot_hybrid_network_modules <- function(hybrid_net,
                                        title = NULL,
                                        module1 = "ensembl_id_1_module",
                                        module2 = "ensembl_id_2_module",
                                        drop_low_freq = FALSE,
                                        drop0 = FALSE,
                                        highlight_gene = NULL,
                                        highlight_module = NULL,
                                        highlight_type=c("BN","CM","WGCN"),
                                        net_layout = "seham",
                                        out_edge_col = "red",
                                        in_edge_col = "green",
                                        gene_annotation) {


#########################
## Define function to specify where in plot to add colored edges
###########################
insertLayer <- function(plotObj, after=1, ...) {
    if (after < 0)
        after <- after + length(plotObj$layers)

    if (!length(plotObj$layers))
        plotObj$layers <- list(...)
    else
        plotObj$layers <- append(plotObj$layers, list(...), after)

    return(plotObj)
}


##################################################################################################################################
  ### Transform hybrid network so that undirected edges are converted to bidirectional edges for the sake of visualization##########
  ##################################################################################################################################
  if (!any(colnames(hybrid_net) == "directed")) {
    hybrid_net$directed <- ifelse(hybrid_net$edge_type %in% c("BN", "CM"), TRUE, FALSE)
  }

  hybrid_net <- as.data.frame(hybrid_net)





  sub <- dplyr::filter(hybrid_net, .data$directed == TRUE)
  int <- dplyr::filter(hybrid_net, .data$directed == FALSE)
  sub2 <- int
  genelist1 <- int$ensembl_id_1
  genelist2 <- int$ensembl_id_2
  int$ensembl_id_1 <- genelist2
  int$ensembl_id_2 <- genelist1
  sub2 <- rbind(sub2, int)
  hybrid_net <- rbind(sub, sub2)

  ### Get HGNC symbols
  gene_df <- gene_annotation |>
    dplyr::select(.data$ensembl_gene_id, .data$hgnc_symbol) %>% data.frame
  rownames(gene_df) <- gene_df$ensembl_gene_id

  hybrid_net$Gene1 <- gene_df[hybrid_net[, "ensembl_id_1"], "hgnc_symbol"]
  hybrid_net$Gene2 <- gene_df[hybrid_net[, "ensembl_id_2"], "hgnc_symbol"]

    hybrid_net$Mod1 <- hybrid_net[,module1]
    hybrid_net$Mod2 <- hybrid_net[,module2]


  ### Remove genes in Module 0/Grey if requested
  if (drop0) {
      hybrid_net <- dplyr::filter(hybrid_net, .data$Mod1 != 0, .data$Mod2 != 0)
      hybrid_net <- dplyr::filter(hybrid_net, !.data$Mod1 %in% c("Grey","grey","Gray","gray"), !.data$Mod2 %in% c("Grey","grey","Gray","gray"))
  } else {
    message("Not removing Module 0/Gray module")
  }

  ##################################################################################################################################
  ### Find gene(s)/module(s) specified for highlighting in hybrid network. Will automatically check if it's ensembl ID or hgnc format. Does not allow for a mixture of the two. Will create edges colored by whether they are entering or exiting your gene(s)/module(s) of interest
  ##################################################################################################################################
  if (!is.null(highlight_gene) & is.null(highlight_module)) {
    message("Highlighting gene")
    hybrid_net$Edge_Color <- "black"
    if (all(grepl("ENSG[0-9]{11}", highlight_gene))) {
      gene_id_col <- c("ensembl_id_1", "ensembl_id_2")
    } else {
      gene_id_col <- c("Gene1", "Gene2")
    }

    index <- lapply(gene_id_col, function(x) intersect(which(hybrid_net[, x] %in% highlight_gene), which(hybrid_net$edge_type %in% highlight_type)))


    if (length(unlist(index)) == 0 & drop0) {
      message("Couldn't find your gene of interest, may not be in network and/or assigned to Module 0. Not highlighting any genes")
    } else if (length(unlist(index)) == 0 & !drop0) {
      message("Couldn't find your gene of interest in the network. Not highlighting any genes")
    } else {
      hybrid_net$Edge_Color[index[[1]]] <- out_edge_col
      hybrid_net$Edge_Color[index[[2]]] <- in_edge_col
    }
  } else if (is.null(highlight_gene) & !is.null(highlight_module)) {
    message("Highlighting module")
    hybrid_net$Edge_Color <- "black"

    index <- lapply(c(module1, module2), function(x) which(hybrid_net[, x] %in% highlight_module))
    if (length(unlist(index)) == 0 & drop0) {
      message("Couldn't find your gene of interest, may not be in network and/or assigned to Module 0. Not highlighting any genes")
    } else if (length(unlist(index)) == 0 & !drop0) {
      message("Couldn't find your gene of interest in the network. Not highlighting any genes")
    } else {
      hybrid_net$Edge_Color[index[[1]]] <- out_edge_col
      hybrid_net$Edge_Color[index[[2]]] <- in_edge_col
    }
  } else {
    message("Not highlighting anything")
    hybrid_net$Edge_Color <- "black"
  }



  hybrid_net_summary <- lapply(unique(hybrid_net$Edge_Color), function(x) as.data.frame(table(dplyr::filter(hybrid_net, .data$Edge_Color == x)[c("Mod1", "Mod2")])))
  for (i in 1:length(hybrid_net_summary)) {
    hybrid_net_summary[[i]]$Edge_Color <- unique(hybrid_net$Edge_Color)[i]
  }

    hybrid_net_summary <- do.call("rbind", hybrid_net_summary)




  hybrid_net_summary[1:2] <- apply(hybrid_net_summary[1:2], 2, function(x) as.character(x))
  hybrid_net_summary[, 4] <- as.character(hybrid_net_summary[, 4])

  colnames(hybrid_net_summary) <- c("Module1", "Module2", "Edge_Count", "Edge_Color")

  hybrid_count <- dplyr::select(hybrid_net, Module1 = .data$Mod1, Gene = .data$ensembl_id_1) |>
    unique() |>
    dplyr::group_by(.data$Module1) |>

      dplyr::summarize(Module_Node_Count = dplyr::n())
    hybrid_count$Module1 <- as.character(hybrid_count$Module1)





  ##################################################################################################################################
  ### Drops low frequency connections between modules if requested. Should only be used if making a less "busy" figure with the caveat that it is not truly representative
  ##################################################################################################################################

  if (drop_low_freq) {
    cutoff_val <- sum(hybrid_count) * 0.01
    hybrid_net_summary <- dplyr::filter(hybrid_net_summary, .data$Edge_Count >= cutoff_val)
  }

  ##################################################################################################################################
  ### Creates "final" network edgelist with all relevant fields. Calculates an edge weight that is currently based on module hyper-connectivity. WIP
  ##################################################################################################################################
  hybrid_net_summary <- dplyr::left_join(hybrid_net_summary, hybrid_count, by = "Module1") |>
    dplyr::group_by(.data$Module1) |>
    dplyr::mutate(Edge_Weight = .data$Edge_Count / .data$Module_Node_Count) |>
    data.frame()

  hybrid_net_summary$Area <- scales::rescale(hybrid_net_summary$Module_Node_Count, to = c(2, 20))

  hybrid_net_summary[1:2] <- apply(hybrid_net_summary[1:2], 2, function(x) as.character(x))


  if (!is.null(highlight_gene)) {
    message("Working gene")



    ##################################################################################################################################
    ### Work around for issues with ggnetwork supporting bipartite graphs. This allows us to plot multiple set of edges involving the same nodes with different values.
    ##################################################################################################################################
    net_list <- lapply(unique(hybrid_net_summary$Edge_Color), function(x) dplyr::filter(hybrid_net_summary, .data$Edge_Color == x))



    net_list[[1]]$Edge_ID <- 1:nrow(net_list[[1]])

    net_out <- network::as.network(unique(dplyr::filter(net_list[[1]], .data$Module1 != .data$Module2, .data$Edge_Count != 0, !is.na(.data$Edge_Count))), multiple = FALSE)
    net_out <- ggnetwork::ggnetwork(net_out, layout = net_layout, arrow.gap = 0.0125)
    if (length(net_list) > 1) {
      for (i in 2:length(net_list)) {
        sub <- net_list[[i]]
        alt <- net_list[[1]]



        int <- dplyr::left_join(dplyr::select(sub, .data$Module1, .data$Module2, .data$Edge_Count, .data$Edge_Color, .data$Module_Node_Count, .data$Area), dplyr::select(alt, .data$Module1, .data$Module2, .data$Edge_ID), by = c("Module1", "Module2"))


        if (all(sub$Module1 == sub$Module2)) {
          int2 <- int
          message("All edges between highlighted gene(s) are within the module")
        } else {
          int2 <- dplyr::filter(net_out, .data$Edge_ID %in% int$Edge_ID)


          int2 <- merge(dplyr::select(int2, -.data$Edge_Count), dplyr::select(int, .data$Edge_ID, .data$Edge_Count), by = "Edge_ID", all.x = TRUE)



          int2$Edge_Color <- unique(sub$Edge_Color)
          net_out <- rbind(net_out, int2)
          ## str(int)
          ## str(int2)
        }

        net_out <- dplyr::filter(net_out, !is.na(.data$Edge_Count))
      }
    } else {
      net_out <- dplyr::filter(net_out, !is.na(.data$Edge_Count))
    }
  } else {
    message("Working module")
    net_out <- network::as.network(unique(dplyr::filter(hybrid_net_summary, .data$Module1 != .data$Module2, .data$Edge_Count != 0, !is.na(.data$Edge_Count))), multiple = FALSE)
    net_out <- ggnetwork::ggnetwork(net_out, layout = net_layout, arrow.gap = 0.03)

    net_out <- dplyr::filter(net_out, !is.na(.data$Edge_Count))
  }




  message("Finished highlighting genes/modules")
  ##################################################################################################################################
  ### Actually plotting the network ################################################################################################
  ##################################################################################################################################

  black <- unique(dplyr::filter(net_out, .data$Edge_Color == "black"))
  red <- unique(dplyr::filter(net_out, .data$Edge_Color == out_edge_col))
  green <- unique(dplyr::filter(net_out, .data$Edge_Color == in_edge_col))


  p <- ggplot2::ggplot(net_out, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
    ggnetwork::geom_edges(data = black, aes( lwd = .data$Edge_Count), color = "black", curvature = 0.1, angle = 10, arrow = arrow(length = unit(4, "pt"), type = "closed", angle = 35),alpha=0.6) +
    theme_void() +
    scale_linewidth_binned(range = c(0.001, 3), n.breaks = 6) +
    scale_size_area("Module_Node_Count", n.breaks = 8, max_size = 20) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggnetwork::geom_nodes(aes(color = .data$vertex.names, size = .data$Area), alpha = 0.8, show.legend = T) +
    ggnetwork::geom_nodes(aes(size = .data$Area * 2, color = .data$vertex.names), alpha = 0.3, show.legend = F)
  if (nrow(red) >= 1) {
    message("Plotting red edges")
    p <- p |>  insertLayer(after=2,ggnetwork::geom_edges(data = red, aes(lwd = 50 * .data$Edge_Count), color = out_edge_col, curvature = 0.2, angle = 90, arrow = arrow(length = unit(4, "pt"), type = "closed", angle = 25), show.legend = F))
  } else {
    p <- p
    message("No red edges")
  }
  if (nrow(green) >= 1) {
    message("Plotting green edges")
    p <- p |>  insertLayer(after=2, ggnetwork::geom_edges(data = green, aes(lwd = 50 * .data$Edge_Count), color = in_edge_col, curvature = 0.3, angle = 90, arrow = arrow(length = unit(4, "pt"), type = "closed", angle = 25), show.legend = F))
  } else {
    p <- p
    message("No green edges")
  }

    p <- p + ggnetwork::geom_nodelabel(data = net_out, aes(label = gsub("Module_", "", .data$vertex.names), size = 0.1 * .data$Area), show.legend = F, fontface = 2, color = "gray50", label.padding = unit(0.1, "lines"))
    p <- p + guides(fill="none",area="none",color="none",size=guide_legend(override.aes=list(fill="gray50")))
  return(p)
}


network_colorscale <- function(color_list = NULL, extension = ".txt", scales = c("fill", "color"), ...) {
  ## usethis::use_package("ggplot2")

  if (length(scales) > 1) {
    scales <- "color"
  }
  if (is.null(color_list) == TRUE) {
           default_cols <- c("grey", "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen", "skyblue3", "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue", "plum2", "thistle2", "thistle1", "salmon4", "palevioletred3", "navajowhite2", "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4", "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2", "yellow4", "skyblue1", "plum", "orangered3", "mediumpurple2", "lightsteelblue", "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", "brown2", "blue2", "darkviolet", "plum3", "thistle3", "thistle", "salmon2", "palevioletred2", "navajowhite1", "magenta4", "lightpink3", "lavenderblush2", "honeydew", "darkseagreen3", "coral", "antiquewhite2", "coral3", "mediumpurple4", "skyblue4", "yellow3")

    cols <- data.frame(c(0:94, paste0("ME", 0:94), paste0("Module_", 0:94),default_cols), default_cols)
        colnames(cols) <- c("Module", "Color")

  } else if (is.character(color_list) & !is.vector(color_list)) {
    if (extension == ".rds" & !is.data.frame(cols)) {
      message("reading RDS")
      cols <- readRDS(color_list)
    } else if (extension == ".txt" & !is.data.frame(cols)) {
      message("reading txt")
      cols <- utils::read.table(color_list, sep = "\t", stringsAsFactors = F, header = T, comment.char = "$")
    }
  }

  if (is.data.frame(cols)) {
    df <- cols
    cols <- df[, 2]
    names(cols) <- df[, 1]
  } else if (is.vector(color_list) == TRUE) {
    cols <- color_list
  }

 # add manual color scales
  if (scales == "color") {
   g <- scale_color_manual(values = cols)
  } else {
   g <- scale_fill_manual(values = cols)
  }
  return(g)
}
