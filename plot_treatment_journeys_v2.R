#' Plot treatment journeys of patients in Tempus data
#'
#' @param cohort Complete path to the Tempus patient cohort or list object created by load_tempus_data (string or list)
#' @param treatment_lines The number of treatment lines to include in the plot for the cohort (Integer)
#' @param num_drug_classes The number of unique drug classes to include in the plot. The plot will only include top "n" classes (Integer)
#' @param num_drug_names The number of unique drug names to include in the plot. The plot will only include top "n" drugs/combinations (Integer)
#' @param drop_untreated Whether to remove patients that are still on the same treatment, no longer under observation and/or deceased as we progress through treatment lines (Boolean)
#' @param plot_option Whether to plot drug classes, individual drug regimens or both for sankey plots ("Drug","Class","Both")
#' @param title Plot title (Optional string)
#' @param show_legend Whether to show the color legend for the plot or not.
#' @param sankey_type Type of plot to produce ("alluvial","sankey") default will be alluvial
#'
#'
#' @return ggplot object if plot_option is set to "Drug" or "Class", Two figures if set to "Both"
#' @note Any filtering or other grouping you wish to perform need to be done before passing the regimen table. Facetting has not been tested and likely will not work. Trying checks again and again and again... seriously will this ever work?? Probably not.
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' plot_treatment_journeys(cohort,
#'   num_drug_names =6
#'   num_drug_classes = 6,
#'   plot_option = "Drug",
#'   title = "Test Plot",
#'   show_legend = FALSE,
#'   sankey_type = "alluvial"
#' )
#' }
#'
plot_treatment_journeys <- function(cohort, treatment_lines = 5, num_drug_classes = 6,
                                    num_drug_names = 5, drop_untreated = FALSE, plot_option = c(
                                      "Drug",
                                      "Class", "Both"
                                      ), title = NULL, show_legend = FALSE, sankey_type ="alluvial",verbose=F) {

#### Defining what we're going to plot 
  if (length(plot_option) > 1) {
    plot_option <- "Both"
  } else {
    plot_option <- plot_option
  }

############  Importing the tempus data and create the plotting dataframe
    if(is.character(cohort)){
        if(verbose){print("Trying to load cohort")}
        input_td <- load_tempus_data(cohort,collection="clinical")
    } else if(is.list(cohort)){
        if(verbose){print("Cohort already loaded, working")}
        input_td <- cohort}

    data_model_check <- any(grepl("regimen", names(input_td)))



    if(data_model_check){
        if(verbose){print("This is a DM1 cohort")}
        regimen_table <- input_td$regimen
        regimen_summary <- regimen_table %>%
            arrange(.data$patient_id,.data$regimen_rank) %>%
            dplyr::select(.data$patient_id,
                          Rank = .data$regimen_rank, Drug_Name = .data$regimen_name, Drug_Class = .data$regimen_class,
      Class_Group = .data$regimen_class_group
    )
    } else{
        if(verbose){print("This is a DM2 cohort")}
        regimen_table <- input_td$onco_line_of_therapy
         regimen_summary <- regimen_table %>%
        arrange(.data$patient_id,.data$line_of_therapy_number) %>%
    dplyr::select(.data$patient_id,
      Rank = .data$line_of_therapy_number, Drug_Name = .data$agents, Drug_Class = .data$therapy_class,
      Class_Group = .data$therapy_class_group
    )
}


 regimen_summary <- regimen_summary %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::mutate(Rank2 = paste0(1:length(.data$Rank), "L")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Rank2) %>%
    dplyr::mutate(Drug_Name2 = forcats::fct_lump(.data$Drug_Name,
      n = num_drug_names
    )) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Rank2) %>%
    dplyr::mutate(Class_Group2 = forcats::fct_lump(.data$Class_Group,
      n = num_drug_classes
    )) %>%
    dplyr::ungroup()
  regimen_summary <- regimen_summary %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::arrange(.data$patient_id, .data$Rank2) %>%
    dplyr::ungroup() %>%
    data.frame()
  regimen_sankey <- regimen_summary %>%
    dplyr::select(.data$patient_id,
      Treatment_Line = .data$Rank2, Drug = .data$Drug_Name2, Class = .data$Class_Group2
    ) %>%
    dplyr::ungroup()
  regimen_sankey_class <- regimen_sankey %>% tidyr::pivot_wider(
    id_cols = .data$patient_id,
    names_from = .data$Treatment_Line, values_from = .data$Class
  )
  regimen_sankey_drug <- regimen_sankey %>% tidyr::pivot_wider(
    id_cols = .data$patient_id,
    names_from = .data$Treatment_Line, values_from = .data$Drug
  )

  all_lines <- unique(regimen_summary$Rank2)
  treatment_cols <- intersect(
    paste0(1:treatment_lines, "L"),
    all_lines
  )
  regimen_sankey_class_final <- regimen_sankey_class %>% ggsankey::make_long(treatment_cols)
  regimen_sankey_drug_final <- regimen_sankey_drug %>% ggsankey::make_long(treatment_cols)
  if (drop_untreated) {
    regimen_sankey_class_final <- regimen_sankey_class_final %>%
      tidyr::drop_na(.data$node)
    regimen_sankey_drug_final <- regimen_sankey_drug_final %>%
      tidyr::drop_na(.data$node)
  } else {
    regimen_sankey_class_final <- regimen_sankey_class_final %>%
      tidyr::replace_na(list(node = "No F/U"))
    regimen_sankey_drug_final <- regimen_sankey_drug_final %>%
      tidyr::replace_na(list(node = "No F/U"))
  }
  title1 <- ifelse(is.null(title), "Patient treatment lines by drug class",
    title
  )
  title2 <- ifelse(is.null(title), "Patient treatment lines by drug name",
    title
  )
  sankey_theme <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 6, angle = 15, hjust = 1
    ), axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 8
    ), strip.text = ggplot2::element_text(
      colour = "black",
      face = "bold", size = 10
    ), plot.title = ggplot2::element_text(hjust = 0.5),
    legend.position = "bottom"
  )
  guide_name <- ggplot2::guides(fill = ggplot2::guide_legend(title = ifelse(plot_option ==
    "Both", "Drug Class/Name", ifelse(plot_option == "Drug",
    "Drug Name", "Drug Class"
  ))))
  if (plot_option == "Drug") {
      if(sankey_type =="alluvial"){
    p2 <- ggplot2::ggplot(regimen_sankey_drug_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_alluvial(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_alluvial_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0
      ) +
      sankey_theme +
      ggplot2::labs(title = title2) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    return(p2)} else if(sankey_type =="sankey"){


                     p2 <- ggplot2::ggplot(regimen_sankey_drug_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_sankey(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_sankey_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0
      ) +
      sankey_theme +
      ggplot2::labs(title = title2) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    return(p2)

                  }
  } else if (plot_option == "Class") {
      if(sankey_type=="alluvial"){
    p <- ggplot2::ggplot(regimen_sankey_class_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_alluvial(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_alluvial_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0
      ) +
      sankey_theme +
      ggplot2::labs(title = title1) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    return(p)} else if(sankey_type =="sankey"){
                 p <- ggplot2::ggplot(regimen_sankey_class_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_sankey(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_sankey_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0
      ) +
      sankey_theme +
      ggplot2::labs(title = title1) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    return(p)

                 }
  } else if (plot_option == "Both") {
      if(sankey_type=="alluvial"){
    p <- ggplot2::ggplot(regimen_sankey_class_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_alluvial(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_alluvial_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0.25, position = ggplot2::position_nudge(x = 0.1)
      ) +
      sankey_theme +
      ggplot2::labs(title = title1) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    print(p)
    p2 <- ggplot2::ggplot(regimen_sankey_drug_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_alluvial(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_alluvial_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0.25, position = ggplot2::position_nudge(x = 0.1)
       ) +
      sankey_theme +
      ggplot2::labs(title = title2) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    print(p2)} else if(sankey_type =="sankey"){
    p <- ggplot2::ggplot(regimen_sankey_class_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_sankey(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_sankey_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0.25, position = ggplot2::position_nudge(x = 0.1)
      ) +
      sankey_theme +
      ggplot2::labs(title = title1) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    print(p)
    p2 <- ggplot2::ggplot(regimen_sankey_drug_final, ggplot2::aes(
      x = .data$x, next_x = .data$next_x,
      node = .data$node, next_node = .data$next_node, label = .data$node,
      fill = as.factor(.data$node)
    )) +
      ggsankey::geom_sankey(
        show.legend = show_legend,
        node.color = 1
      ) +
      ggsankey::geom_sankey_label(
                size = 2, color = "black",show.legend=F,fontface=2,
        hjust = 0.25, position = ggplot2::position_nudge(x = 0.1)
       ) +
      sankey_theme +
      ggplot2::labs(title = title2) +
      ggplot2::xlab("Treatment Line") +
      ggplot2::ylab("Num. Patients") +
      guide_name
    print(p2)

                 }
  } else {
    print("Not a valid plotting option")
    stop()
  }
}
