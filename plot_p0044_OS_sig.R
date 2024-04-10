#' Plot p0044 genewise overall survival associations
#'
#' @param os_signature Data frame or path to rds file of the p0044 os_signatures output for the cohort of interest
#' @param wgcna_gene_info The p0044 module assignments for the cohort of interest
#' @param title Plot title (Optional string)
#' @param highlight_gene Which genes you wish to highlight in the plot
#'
#'
#'
#' @return ggplot object of OS signatures
#' @note
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' plot_os_signatures(os_signature,
#'   wgcna_gene_info,
#'   title = "Test Plot",
#'   highlight_gene = "TP53"
#' )
#' }
#'
plot_os_signatures <- function(os_signature, wgcna_gene_info, title = NULL, highlight_gene = NULL) {
    suppressWarnings()
  if (is.character(os_signature) && file.exists(os_signature)) {
    os_sig <- readRDS(os_signature)
  } else if (is.data.frame(os_signature)) {
    os_sig <- os_signature
  }
  if (is.character(wgcna_gene_info) && file.exists(wgcna_gene_info)) {
    gene_info <- read.table(wgcna_gene_info, header = T, sep = "\t", stringsAsFactors = F)
  } else if (is.data.frame(wgcna_gene_info)) {
    gene_info <- wgcna_gene_info
  }


  os_sig <- left_join(os_sig, network_outputs$wgcna_gene[, 1:3], by = "ensembl_gene_id")
  os_sig <- dplyr::filter(os_sig,!is.na(module_label))

  highlight_df <- rbind(dplyr::filter(os_sig, -log10(adj.p) >= 4, abs(estimate) >= 1.25), dplyr::filter(os_sig, -log10(adj.p) >= 4, abs(estimate) <= 0.75))

  if (!is.null(title)) {
    title <- title
  } else {
    title <- title <- "OS associations by module"
  }

  if (!is.null(highlight_gene)) {
    gene_highlight_df <- dplyr::filter(os_sig, hgnc_symbol %in% highlight_gene)
    p <- ggplot(os_sig, aes(as.factor(module_label), log(estimate), group = as.factor(module_label), label = hgnc_symbol)) +
      geom_boxplot(alpha = 0.6, width = 0.3, outlier.size = 0.01) +
      geom_point(aes(size = -log10(adj.p)), shape = 21, color = "gray50", alpha = 0.5, position = position_jitterdodge(seed = 10), fill = os_sig$module_color) +
      theme_anf +
      xlab("Module") +
      labs(title = title) +
      ylab("log(Hazard Ratio)") +
      geom_text_repel(data = highlight_df, size = 3, fontface = "bold", color = "gray50", position = position_dodge(width=0.2), min.segment.length = 0.1) +
      geom_text_repel(data = gene_highlight_df, size = 4, color = "red", fontface = "bold", position = position_dodge(width=0.2), min.segment.length = 0.1) +
      geom_hline(yintercept = 0, lty = 2, lwd = 2, color = "red", alpha = 0.4)
  } else{
      p <- ggplot(os_sig, aes(as.factor(module_label), log(estimate), group = as.factor(module_label), label = hgnc_symbol)) +
      geom_boxplot(alpha = 0.6, width = 0.3, outlier.size = 0.01) +
      geom_point(aes(size = -log10(adj.p)), shape = 21, color = "gray50", alpha = 0.5, position = position_jitterdodge(seed = 10), fill = os_sig$module_color) +
      theme_anf +
      xlab("Module") +
      labs(title = title) +
      ylab("log(Hazard Ratio)") +
      geom_text_repel(data = highlight_df, size = 3, fontface = "bold", color = "gray50", position = position_dodge(width=0.2), min.segment.length = 0.1) +
         geom_hline(yintercept = 0, lty = 2, lwd = 2, color = "red", alpha = 0.4)
  }

return(p)
}
