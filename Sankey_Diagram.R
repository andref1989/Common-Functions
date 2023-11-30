plot_sankey <- function(regimen_table,treatment_lines=5, num_drug_classes=6,num_drug_names=5,drop_untreated=FALSE,plot_option=c("Drug","Class","Both")){
    require(forcats)
    require(dplyr)
    require(tidyr)
    require(ggsankey)


    if(length(plot_option) >1){
        plot_option <- "Both"} else{ plot_option <- plot_option}
    regimen_summary <- regimen_table %>%
        dplyr::select(patient_id,Rank=regimen_rank,Drug_Name=regimen_name,
                      Drug_Class=regimen_class,Class_Group=regimen_class_group) %>%
        group_by(patient_id) %>% mutate(Rank2=paste0(1:length(Rank),"L")) %>% ungroup %>%
        group_by(Rank2) %>% mutate(Drug_Name2 = forcats::fct_lump(Drug_Name, n = num_drug_names)) %>%
        ungroup %>% group_by(Rank2) %>%
        mutate(Class_Group2=forcats::fct_lump(Class_Group, n=num_drug_classes)) %>%
        ungroup
    regimen_summary <- regimen_summary %>% group_by(patient_id) %>%
        arrange(patient_id,Rank2) %>% ungroup %>% data.frame


    regimen_sankey <- regimen_summary %>% dplyr::select(patient_id,Treatment_Line=Rank2,Drug=Drug_Name2,Class=Class_Group2) %>% ungroup

    regimen_sankey_class <- regimen_sankey %>% tidyr::pivot_wider(id_cols=patient_id,names_from=Treatment_Line,values_from=Class)
    regimen_sankey_drug <- regimen_sankey %>% tidyr::pivot_wider(id_cols=patient_id,names_from=Treatment_Line,values_from=Drug)


    all_lines <- unique(regimen_summary$Rank2)
    treatment_cols <- intersect(paste0(1:treatment_lines,"L"),all_lines)


    regimen_sankey_class_final <- regimen_sankey_class %>% make_long(treatment_cols)
    regimen_sankey_drug_final <- regimen_sankey_drug %>% make_long(treatment_cols)



    sankey_theme <- theme(axis.text.x=element_text(face='bold',size=6,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=10),plot.title=element_text(hjust=0.5))


    if(drop_untreated){
        regimen_sankey_class_final <- regimen_sankey_class_final %>% drop_na(node)
        regimen_sankey_drug_final <- regimen_sankey_drug_final %>% drop_na(node)} else{
                                                                                    regimen_sankey_class_final <- regimen_sankey_class_final %>% replace_na(list(node="No F/U"))
                                                                                    regimen_sankey_drug_final <- regimen_sankey_drug_final %>% replace_na(list(node="No F/U"))
                                                                                    }

    if(plot_option=="Drug"){
        p2 <- ggplot(regimen_sankey_drug_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0)+sankey_theme+labs(title="Patient treatment lines by drug name")+xlab("Treatment Line")+ylab("Num. Patients")
        print(p2)} else if(plot_option=="Class"){
                     p <- ggplot(regimen_sankey_class_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0)+sankey_theme+labs(title="Patient treatment lines by drug class")+xlab("Treatment Line")+ylab("Num. Patients")
                     print(p)
} else if(plot_option=="Both"){

    p <- ggplot(regimen_sankey_class_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0.25,position=position_nudge(x=0.1))+sankey_theme+labs(title="Patient treatment lines by drug class")+xlab("Treatment Line")+ylab("Num. Patients")
                     print(p)

    p2 <- ggplot(regimen_sankey_drug_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0.25,position=position_nudge(x=0.1))+sankey_theme+labs(title="Patient treatment lines by drug name",)+xlab("Treatment Line")+ylab("Num. Patients")
    print(p2)
} else { print("Not a valid plotting option"); stop()}

    }
