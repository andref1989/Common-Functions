#' Select from the filtered MMF (Master Molecular File)
#'
#' Select the genes and the variants of interest.
#'
#' @note The Master Molecular File can be huge, so it would be a good candidate
#' for loading into database so we can query and only pull the needed rows.
#' As is here, this is not much of a function, just a wrapper filter.
#'
#' @param mmf required, the molecular_master_file or filtered table
#' @param analysis_ids optional, dna analysis ids
#' @param somatic_germline_select default c("S","G") to keep both somatic and
#' germline short variants. Allows other variant types (like copy number) to
#' pass.
#' @param functional_impact_select required, The default is c("P","LP"), other
#' options ("B","LB","US","CE"). Allows other variant types (like copy number)
#' to pass.
#' @param variant_types_select required, one or more of ('SHRTVRNT','CNALTER',
#' 'MSISTAT','REARRANG','TMB').
#' The default is SHRTVRNT (short variant)
#' @param genes_of_interest optional, a character vector of gene symbols from
#' the `gene_canonical_name` variable in MMF
#' @param genes_to_exclude optional, a character vector of gene symbols from the
#' `gene_canonical_name` variable in MMF to exclude (the default is c("UGT1A1",
#' "TAP1","TAP2","TNF","DAXX"))
#'
#' @return all columns from `molecular_master_file_filtered` (or whatever you
#' input as the mmf argument)
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Select a gene (EGFR) for ALL analysis_ids
#' # default returns short variants, P/LP, S/G
#' selected_mmf <- select_mmf(mydata$g_molecular_master_file_filtered,
#'                            genes_of_interest = "EGFR")
#'
#' }
select_mmf <- function(
                    cohort,
                    analysis_ids = NULL,
                    somatic_germline_select = c("S","G","somatic","germline"),
                    functional_impact_select = c("P","LP"),
                    variant_types_select = c("SHRTVRNT"),
                    genes_of_interest = NULL,
                    genes_to_exclude = c("UGT1A1","DAXX","TNF","TAP1","TAP2")){

  somatic_germline_select <- match.arg(somatic_germline_select,
                                       several.ok = TRUE)
  functional_impact_select <- match.arg(functional_impact_select,
                                        choices = c("P",
                                                    "LP",
                                                    "US",
                                                    "B",
                                                    "LB",
                                                    "CE"),
                                        several.ok = TRUE)
  if(!is.null(genes_of_interest)) {
    stopifnot(is.character(genes_of_interest))
  }
  if(!is.null(genes_to_exclude)) {
    stopifnot(is.character(genes_to_exclude))
  }


############  Importing the tempus data and create the plotting dataframe
    if(is.character(cohort)){

        input_td <- load_tempus_data(cohort,collection="molecular")
    } else if(is.list(cohort)){

        input_td <- cohort}

    data_model_check <- any(grepl("molecular_master_file", names(input_td)))

    if(data_model_check){


        g_molecular_master_file <- input_td[[grep("molecular_master_file", names(input_td))[1]]]
          variant_types_select <- match.arg(variant_types_select,
                                    choices = c('SHRTVRNT',
                                                'CNALTER',
                                                'MSISTAT',
                                                'REARRANG',
                                                'TMB'),
                                    several.ok = TRUE)


    }
    else {
        g_molecular_master_file <- input_td[[grep("onco_result_snv", names(input_td))[1]]] %>% dplyr::select(patient_id, analysis_id, assay, somatic_germline=variant_origin, gene_canonical_name=gene_symbol, functional_impact=variant_classification,variant_type_code=variant_type_detailed)
          variant_types_select <- match.arg(variant_types_select,
                                    choices = c('snv','SHRTVRNT',
                                                'CNALTER',
                                                'MSISTAT',
                                                'REARRANG',
                                                'TMB'),
                                    several.ok = TRUE)


        }


  filtered <- mmf %>%
    dplyr::filter(
      (.data$somatic_germline %in% somatic_germline_select |
         .data$variant_type_code %in% c('CNALTER','MSISTAT','REARRANG','TMB')),
      .data$variant_type_code %in% variant_types_select,
      (.data$functional_impact %in% functional_impact_select |
        .data$variant_type_code %in% c('CNALTER','MSISTAT','REARRANG','TMB'))
    )
  if (!is.null(genes_of_interest)){
    filtered <- dplyr::filter(filtered,
                              .data$gene_canonical_name %in% genes_of_interest)
  }
  if (!is.null(genes_to_exclude)){
    filtered <- dplyr::filter(filtered,
                              !.data$gene_canonical_name %in% genes_to_exclude)
  }
  if (!is.null(analysis_ids)){
    filtered <- dplyr::filter(filtered, .data$analysis_id %in% analysis_ids)
  }
  filtered
}
