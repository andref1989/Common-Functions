calc_overlap_cluster_samples <- function(clusters_1,
                                         name_cluster_1,
                                         clusters_2,
                                         name_cluster_2) {

                                        # get list of all clusters

    list_clusters_1 <- unique(clusters_1[[name_cluster_1]])
    list_clusters_2 <- unique(clusters_2[[name_cluster_2]])

                                        # get number of samples in each cluster

    num_samples_cluster_1 <- clusters_1 |>
        dplyr::select(analysis_id_rna,
                      cluster_1 = !!rlang::sym(name_cluster_1)) |>
        dplyr::count(cluster_1,
                     name = "num_samples_cluster_1")

    num_samples_cluster_2 <- clusters_2 |>
        dplyr::select(analysis_id_rna,
                      cluster_2 = !!rlang::sym(name_cluster_2)) |>
        dplyr::count(cluster_2,
                     name = "num_samples_cluster_2")

                                        # join cluster information into 1 data frame

    clusters_all <- dplyr::inner_join(clusters_1 |>
                                      dplyr::select(analysis_id_rna,
                                                    cluster_1 = !!rlang::sym(name_cluster_1)),
                                      clusters_2 |>
                                      dplyr::select(analysis_id_rna,
                                                    cluster_2 = !!rlang::sym(name_cluster_2)),
                                      by = "analysis_id_rna")

                                        # calculate overlap

    overlap <- clusters_all |>
        dplyr::group_by(cluster_1,
                        cluster_2) |>
        dplyr::summarise(num_samples = dplyr::n(),
                         .groups = "drop") |>
        dplyr::full_join(expand.grid(cluster_1 = list_clusters_1,
                                     cluster_2 = list_clusters_2),
                         by = c("cluster_1",
                                "cluster_2")) |>
        dplyr::mutate(num_samples = tidyr::replace_na(num_samples,
                                                      0))

                                        # add column for the number of samples in each cluster

    overlap <- overlap |>
        dplyr::left_join(num_samples_cluster_1,
                         by = "cluster_1") |>
        dplyr::left_join(num_samples_cluster_2,
                         by = "cluster_2") |>
        dplyr::mutate(percent_samples_cluster_1 = num_samples / num_samples_cluster_1 * 100,
                      percent_samples_cluster_2 = num_samples / num_samples_cluster_2 * 100)

                                        # return overlap

    return(overlap)

}

plot_overlap_cluster_samples <- function(overlap,
                                         title_axis_x,
                                         title_axis_y) {

                                        # get list of clusters

    list_clusters_1 <- sort(unique(overlap$cluster_1))
    list_clusters_2 <- sort(unique(overlap$cluster_2))

                                        # define coordinates for bottom triangle

    coordinates_triangle_bottom <- dplyr::bind_rows(
                                        # lower left corner of triangle
                                              expand.grid(cluster_1 = list_clusters_1,
                                                          cluster_2 = list_clusters_2) |>
                                              dplyr::mutate(x = cluster_1,
                                                            y = cluster_2),
                                        # lower right corner of triangle
                                              expand.grid(cluster_1 = list_clusters_1,
                                                          cluster_2 = list_clusters_2) |>
                                              dplyr::mutate(x = cluster_1 + 1,
                                                            y = cluster_2),
                                        # upper right corner of triangle
                                              expand.grid(cluster_1 = list_clusters_1,
                                                          cluster_2 = list_clusters_2) |>
                                              dplyr::mutate(x = cluster_1 + 1,
                                                            y = cluster_2 + 1)
                                          )

                                        # define coordinates for the top triangle

    coordinates_triangle_top <- dplyr::bind_rows(
                                        # lower left corner of triangle
                                           expand.grid(cluster_1 = list_clusters_1,
                                                       cluster_2 = list_clusters_2) |>
                                           dplyr::mutate(x = cluster_1,
                                                         y = cluster_2),
                                        # upper left corner of triangle
                                           expand.grid(cluster_1 = list_clusters_1,
                                                       cluster_2 = list_clusters_2) |>
                                           dplyr::mutate(x = cluster_1,
                                                         y = cluster_2 + 1),
                                        # upper right corner of triangle
                                           expand.grid(cluster_1 = list_clusters_1,
                                                       cluster_2 = list_clusters_2) |>
                                           dplyr::mutate(x = cluster_1 + 1,
                                                         y = cluster_2 + 1)
                                       )

                                        # add percentages to triangle coordinates

    coordinates_triangle_bottom <- overlap |>
        dplyr::select(cluster_1,
                      cluster_2,
                      percent = percent_samples_cluster_1) |>
        dplyr::right_join(coordinates_triangle_bottom,
                          by = c("cluster_1",
                                 "cluster_2"))

    coordinates_triangle_top <- overlap |>
        dplyr::select(cluster_1,
                      cluster_2,
                      percent = percent_samples_cluster_2) |>
        dplyr::right_join(coordinates_triangle_top,
                          by = c("cluster_1",
                                 "cluster_2"))

                                        # define coordinates for text for bottom triangle

    coordinates_text_bottom <- expand.grid(cluster_1 = list_clusters_1,
                                           cluster_2 = list_clusters_2) |>
        dplyr::mutate(x = cluster_1 + 0.9,
                      y = cluster_2 + 0.1)

                                        # define coordinates for text for top triangle

    coordinates_text_top <- expand.grid(cluster_1 = list_clusters_1,
                                        cluster_2 = list_clusters_2) |>
        dplyr::mutate(x = cluster_1 + 0.1,
                      y = cluster_2 + 0.9)

                                        # add percentages to text coordinates

    coordinates_text_bottom <- overlap |>
        dplyr::select(cluster_1,
                      cluster_2,
                      percent = percent_samples_cluster_1) |>
        dplyr::mutate(label = sprintf("%.1f%%",
                                      percent)) |>
        dplyr::right_join(coordinates_text_bottom,
                          by = c("cluster_1",
                                 "cluster_2"))

    coordinates_text_top <- overlap |>
        dplyr::select(cluster_1,
                      cluster_2,
                      percent = percent_samples_cluster_2) |>
        dplyr::mutate(label = sprintf("%.1f%%",
                                      percent)) |>
        dplyr::right_join(coordinates_text_top,
                          by = c("cluster_1",
                                 "cluster_2"))

                                        # create tick labels

    labels_x <- overlap |>
        dplyr::arrange(cluster_1) |>
        dplyr::mutate(label = sprintf("cluster %s\n(%.0f samples)",
                                      cluster_1,
                                      num_samples_cluster_1)) |>
        dplyr::pull(label) |>
        unique()

    labels_y <- overlap |>
        dplyr::arrange(cluster_2) |>
        dplyr::mutate(label = sprintf("cluster %s\n(%.0f samples)",
                                      cluster_2,
                                      num_samples_cluster_2)) |>
        dplyr::pull(label) |>
        unique()

                                        # create the plot

    plot_overlap <- ggplot2::ggplot() +
        ggplot2::geom_polygon(data = coordinates_triangle_bottom,
                              mapping = ggplot2::aes(x = x,
                                                     y = y,
                                                     group = interaction(cluster_1,
                                                                         cluster_2),
                                                     fill = percent),
                              color = "white") +
        ggplot2::geom_polygon(data = coordinates_triangle_top,
                              mapping = ggplot2::aes(x = x,
                                                     y = y,
                                                     group = interaction(cluster_1,
                                                                         cluster_2),
                                                     fill = percent),
                              color = "white") +
        ggplot2::geom_text(data = coordinates_text_bottom,
                           mapping = ggplot2::aes(x = x,
                                                  y = y,
                                                  label = label),
                           color = "white",
                           hjust = 1,
                           vjust = 0) +
        ggplot2::geom_text(data = coordinates_text_top,
                           mapping = ggplot2::aes(x = x,
                                                  y = y,
                                                  label = label),
                           color = "white",
                           hjust = 0,
                           vjust = 1) +
        viridis::scale_fill_viridis(option = "plasma",
                                    direction = -1,
                                    limits = c(0, 100)) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank()) +
        ggplot2::coord_cartesian(expand = 0) +
        ggplot2::labs(x = title_axis_x,
                      y = title_axis_y,
                      title = "Percentage of samples in common",
                      fill = "percentage\nof samples")

                                        # return the plot


    if(!any(is.character(c(list_clusters_1,list_clusters_2)))){
        plot_overlap <- plot_overlap +
                ggplot2::geom_vline(xintercept = seq(min(list_clusters_1),
                                             max(list_clusters_1) + 1,
                                             1),
                            color = "white",
                            linewidth = 2) +
        ggplot2::geom_hline(yintercept = seq(min(list_clusters_2),
                                             max(list_clusters_2) + 1,
                                             1),
                            color = "white",
                            linewidth = 2) +
                    ggplot2::scale_x_continuous(breaks = seq(min(list_clusters_1),
                                                 max(list_clusters_1),
                                                 1) + 0.5,
                                    labels = labels_x) +
        ggplot2::scale_y_continuous(breaks = seq(min(list_clusters_2),
                                                 max(list_clusters_2),
                                                 1) + 0.5,
                                    labels = labels_y)} else {
                                                                  plot_overlap <- plot_overlap +
                ggplot2::geom_vline(xintercept = sort(list_clusters_1),
                            color = "white",
                            linewidth = 2) +
        ggplot2::geom_hline(yintercept = sort(list_clusters_2),
                            color = "white",
                            linewidth = 2) +
                    ggplot2::scale_x_continuous(breaks = sort(list_clusters_1),
                                    labels = labels_x) +
            ggplot2::scale_y_continuous(breaks = sort(list_clusters_2),
                                        labels = labels_y)}

    return(plot_overlap)

    }
