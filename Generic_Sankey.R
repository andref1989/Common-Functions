#' Plot generic sequential event data and/or treatment paths in non-Tempus data
#'
#' @param df Data frame of sequential events (data.frame)
#' @param max_order The number of sequential events to include in the plot for the cohort (Integer)
#' @param num_features The number of unique drug names/features/events to include in the plot. The plot will only include top "n" features(Integer)
#' @param drop_unfollowed Whether to remove patients that are  no longer under observation and/or deceased as we progress through sequence of events (Boolean)
#' @param title Plot title (Optional string)
#' @param show_legend Whether to show the color legend for the plot or not.(Boolean)
#' @param order_column The column containing the sequential order of events (per individual data required) (string)
#' @param patient_identifier The column containing the individual identifiers (string)
#' @param feature_name The name of the features being plotting (string)
#' @param sankey_type Type of plot to produce ("alluvial","sankey") default will be alluvial
#' @param show_label Whether to include the drug/class labels in the plot (Boolean)
#' @param na_replacement What to fill in when data is no longer available for a group of entities. (string) Default: "No F/U"
#'
#'
#' @return ggplot object
#' @note Any filtering or other grouping you wish to perform need to be done before passing the dataframe. Facetting has not been tested and likely will not work.
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' plot_generic_sanker(df,
#'   max_order = 5,
#'   num_features = 6,
#'   drop_unfollowed = TRUE,
#'   title = "Test Plot",
#'   show_legend = FALSE,
#'   order_column = "Rank",
#'   patient_identifier = "patient_id",
#'   feature_name = "Drug_Name",
#'   sankey_type="alluvial")
#' }

plot_generic_sankey <- function(df, max_order = 5, num_features = 5, drop_unfollowed = FALSE, title=NULL, show_legend = FALSE, order_column = "Rank",patient_identifier = "patient_id", feature_name= "Drug_Name", sankey_type="alluvial",show_label=F,na_replacement="No F/U"){


  ## require(forcats)
  ## require(dplyr)
  ## require(tidyr)
  ## require(ggsankey)


    df2 <- df[c(patient_identifier, order_column, feature_name)]
  colnames(df2) <- c("patient_id", "Rank", "Feature_Name")
  df_summary <- df2 %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::mutate(Rank2 = paste0(1:length(.data$Rank), "L")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Rank2) %>%
    dplyr::mutate(Feature_Name2 = forcats::fct_lump_n(.data$Feature_Name, n = num_features)) %>% dplyr::ungroup()


  df_summary <- df_summary %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::arrange(.data$patient_id, .data$Rank2) %>%
    dplyr::ungroup() %>%
    data.frame()

  df_summary <- dplyr::select(df_summary, .data$patient_id, Order = .data$Rank2, Feature = .data$Feature_Name2)
  ##        str(df_summary)
  df_sankey <- df_summary %>% tidyr::pivot_wider(id_cols = .data$patient_id, names_from = .data$Order, values_from = .data$Feature)
  ## str(df_sankey)
  all_lines <- unique(df_summary$Order)
  order_cols <- intersect(paste0(1:max_order, "L"), all_lines)

  df_sankey_final <- df_sankey %>% ggsankey::make_long(order_cols)
  if (drop_unfollowed) {
    df_sankey_final <- df_sankey_final %>% tidyr::drop_na(.data$node)
  } else {
    df_sankey_final <- df_sankey_final %>% tidyr::replace_na(list(node = na_replacement))
  }

  title1 <- ifelse(is.null(title), "", title)
  sankey_theme <- ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", size = 10, angle = 15, hjust = 1), axis.text.y = ggplot2::element_text(face = "bold", size = 10), strip.text = ggplot2::element_text(colour = "black", face = "bold", size = 10), plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom")

  guide_name <- ggplot2::guides(fill = ggplot2::guide_legend(title = feature_name))

if(sankey_type=="alluvial"){

    p2 <- ggplot2::ggplot(df_sankey_final, ggplot2::aes(x = .data$x, next_x = .data$next_x, node = .data$node, next_node = .data$next_node, label = .data$node, fill = as.factor(.data$node))) +
    ggsankey::geom_alluvial(show.legend = show_legend, node.color = 1, space = 3,flow.alpha=0.6) +
    sankey_theme +
    ggplot2::labs(title = title1) +
    ggplot2::xlab(order_column) +
    ggplot2::ylab("Num. Patients") +
    guide_name
    if(show_label){ p2 <- p2 + ggsankey::geom_alluvial_label(size = 2, fontface=2, color = "black", hjust = 0, space = 3, show.legend = F)}
} else if(sankey_type=="sankey"){
  p2 <- ggplot2::ggplot(df_sankey_final, ggplot2::aes(x = .data$x, next_x = .data$next_x, node = .data$node, next_node = .data$next_node, label = .data$node, fill = as.factor(.data$node))) +
    ggsankey::geom_sankey(show.legend = show_legend, node.color = 1, space = 3,flow.alpha=0.6) +
    sankey_theme +
    ggplot2::labs(title = title1) +
    ggplot2::xlab(order_column) +
    ggplot2::ylab("Num. Patients") +
    guide_name

  if(show_label){ p2 <- p2 + ggsankey::geom_sankey_label(size = 2, fontface=2,color = "black", hjust = 0, space = 3, show.legend = F)}

    }


  return(p2)
}
