#' Plot generic sequential event data and/or treatment paths in non-Tempus data
#'
#' @param df Data frame of sequential events (data.frame)
#' @param max_order The number of sequential events to include in the plot for the cohort (Integer) Default: 5
#' @param num_features The number of unique drug names/features/events to include in the plot. The plot will only include top "n" features(Integer) Default: 5
#' @param drop_unfollowed Whether to remove patients that are  no longer under observation and/or deceased as we progress through sequence of events (Boolean) default: FALSE
#' @param title Plot title (Optional string) Default:NULL
#' @param show_legend Whether to show the color legend for the plot or not.(Boolean) Default:False
#' @param order_column The column containing the sequential order of events (per individual data required) (string) Default:"Rank"
#' @param patient_identifier The column containing the individual identifiers (string) Default:"patient_id"
#' @param feature_name The name of the features being plotting (string) Default: "Drug_Name"
#' @param show_label Whether to show the text labels in the plot for each element (Boolean) Default: True
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
#'   max_order= 5,
#'   num_features = 6,
#'   drop_unfollowed = TRUE,
#'   title = "Test Plot",
#'   show_legend = FALSE,
#'   order_column = "Rank",
#'   patient_identifier = "patient_id",
#'   feature_name = "Treatment Line"
#' )
#' }
#'


plot_generic_sankey <- function(df, max_order=5,num_features=5,drop_unfollowed=FALSE,title=NULL,show_legend=FALSE,order_column="Rank",patient_identifier="patient_id",feature_name="Drug_Name",show_label=TRUE,na_replacement="No F/U"){


    require(forcats)
    require(dplyr)
    require(tidyr)
    require(ggsankey)


    df2 <- df[c(patient_identifier,order_column,feature_name)]
    colnames(df2) <- c("patient_id","Rank","Feature_Name")
    df_summary <- df2 %>% group_by(patient_id) %>% mutate(Rank2=paste0(1:length(Rank),"L")) %>% ungroup %>%
        group_by(Rank2) %>% mutate(Feature_Name2 = forcats::fct_lump_n(Feature_Name, n = num_features)) %>%
        ungroup

    df_summary <- df_summary %>% group_by(patient_id) %>% arrange(patient_id,Rank2) %>% ungroup %>% data.frame

    df_summary <- dplyr::select(df_summary,patient_id, Order=Rank2,Feature=Feature_Name2)
##        str(df_summary)
    df_sankey <- df_summary %>% tidyr::pivot_wider(id_cols=patient_id,names_from=Order,values_from=Feature)
##str(df_sankey)
    all_lines <- unique(df_summary$Order)
    order_cols <- intersect(paste0(1:max_order,"L"),all_lines)

    df_sankey_final <- df_sankey %>% make_long(order_cols)
    if(drop_unfollowed){
        df_sankey_final <- df_sankey_final %>% drop_na(node)} else {
                                                                drop_index <- intersect(which(!is.na(df_sankey_final$node)),which(is.na(df_sankey_final$next_node)))
                                                                df_sankey_final$next_node[drop_index] <- na_replacement

                                                                for(i in drop_index){ check <- as.character(df_sankey_final$x[i]) > as.character(df_sankey_final$x[i+1])
                                                                    if(!check){ df_sankey_final$node[i+1] <- na_replacement}}


                                                                df_sankey_final <- df_sankey_final %>% tidyr::drop_na(.data$node)
                                                            }

    title1 <- ifelse(is.null(title),"",title)
    sankey_theme <- theme(axis.text.x=element_text(face='bold',size=10,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=10),strip.text = element_text(colour = "black", face = "bold",size=10),plot.title=element_text(hjust=0.5),legend.position="bottom")

    guide_name <- ggplot2::guides(fill = ggplot2::guide_legend(title = feature_name))

if(show_label){
    p2 <- ggplot(df_sankey_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_alluvial(show.legend=show_legend,node.color=1,space=3)+geom_alluvial_label(size = 3, fontface="bold",color = "black", hjust = 0,space=3,show.legend=F)+sankey_theme+labs(title=title1)+xlab(order_column)+ylab("Num. Patients")+guide_name } else {     p2 <- ggplot(df_sankey_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,fill=as.factor(node)))+geom_alluvial(show.legend=show_legend,node.color=1,space=3)+sankey_theme+labs(title=title1)+xlab(order_column)+ylab("Num. Patients")+guide_name}




    return(p2)

    }
