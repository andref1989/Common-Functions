parse_txt <- function(object) {

  data <- httr::content(object,
                        as = "raw")

  data <- readr::read_tsv(data,
                          col_names = FALSE,
                          show_col_types = FALSE)

  colnames(data) <- c("source",
                      "target")

  return(data)

}


## googleCloudStorageR::gcs_get_object(bucket = bucket_name,
##                                                          object = path,
##                                     parseFunction = parse_txt)
