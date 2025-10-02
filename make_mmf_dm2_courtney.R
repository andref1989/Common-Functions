make_mmf_dm2 <- function(td){
   variant_table <- td$onco_result_overview_gene %>%
      # filter(is_gene_significant == 1) %>%
      left_join(
        td$onco_result_snv_indel_passing_filtered %>%
          filter(variant_classification %in% c("LP", "P")) %>%
          mutate(variant_type_code = "SHRTVRNT") %>%
          select(patient_id, tumor_biospecimen_id, gene_symbol,
                 variant_type_code, variant_molecular_consequence,
                 variant_origin, variant_classification),
        by = c("patient_id", "tumor_biospecimen_id", "gene_symbol")
      ) %>%
      left_join(
        td$onco_result_cnv_gene %>%
          # filter(copy_number != 0) %>%
          # filter(is_significant == 1) %>%
          mutate(variant_type_code = "CNALTER") %>%
          select(patient_id, tumor_biospecimen_id, gene_symbol,
                 copy_number, variant_type_code),
        by = c("patient_id", "tumor_biospecimen_id", "gene_symbol")
      ) %>%
      mutate(
        variant_type_code = coalesce(variant_type_code.x, variant_type_code.y),
        result = case_when(
          is.na(copy_number) ~ "Neutral",
          copy_number < 2    ~ "Loss",
          copy_number > 2    ~ "Gain",
          copy_number == 2   ~ "Neutral"
        ),
        biospecimen_id = tumor_biospecimen_id) %>%
      select(-variant_type_code.x, -variant_type_code.y) %>%
      filter(
        grepl("xT|xO|xE", assay),
        gene_symbol %in% xt_genes
      )
   return(variant_table)
}
