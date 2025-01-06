list_p44_networks <- function(directory="result_010"){

    library(googleCloudStorageR)
    library(gargle)
    output_list <- list()
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)



    p44_result_version <- directory

    networks <- googleCloudStorageR::gcs_list_objects(bucket = "pathos-research-results",
                                                  prefix = sprintf("project_0044/%s",
                                                                   p44_result_version))

    networks <- networks |>
  dplyr::mutate(network = stringr::str_split_i(name,
                                               "\\/",
                                               3)) |>
  dplyr::select(network) |>
  dplyr::distinct()

# exclude CCLE networks

networks <- networks |>
  dplyr::filter(network != "CCLE")



#### determine which networks have a BN ####
    networks <- dplyr::filter(networks, grepl("[0-9]{2}Q[0-9]",network)) %>% arrange(desc(network))
    return(networks)
}
