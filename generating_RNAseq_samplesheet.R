library(googleCloudStorageR)
library(gargle)
library(tidyr)
library(dplyr)
library(glue)

scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)
gs_path <- "gs://pathos-data-internal/"

# read in sample details from experimental design
novogene_isolation <- readxl::read_xlsx("~/Downloads/11-12 poci sample ids plus library qc from novogene_pooled.xlsx", 
                                        sheet = 1, 
                                        skip = 1)
duke_isolation <- readxl::read_xlsx("~/Downloads/11-12 poci sample ids plus library qc from novogene_pooled.xlsx", 
                                    sheet = 1, 
                                    range = "R2:AD95") |> 
  drop_na(`Sample Number`)
  # we want the sample number, cell line, experimental condition/treatment, experimental concentration, RIN

# A_1 indicates that the sample was originally isolated at NovoGene and passed QC
# all other A samples correspond to the A from the second column of samples

# need to parse for the two different columns of data
experimental_design <- novogene_isolation |> 
  dplyr::select(sample = `Sample Number...1`, 
                cell_line = `Cell Line...2`, 
                treatment = `Experimental Conditions...3`, 
                concentration = `Experimental Conditions...4`, 
                rin = `RIN...11`, 
                qc = `Sample QC Results...12`) |> 
  mutate(sample = ifelse(grepl(sample, pattern = "A"), paste0(sample, "_1"), sample)) |> 
  mutate(rna_isolation = "Novogene") |> 
  rbind(duke_isolation |> 
          dplyr::select(sample = `Sample Number`, 
                        cell_line = `Cell Line`, 
                        treatment = `Experimental Conditions...3`, 
                        concentration = `Experimental Conditions...4`, 
                        rin = `RIN`, 
                        qc = `Sample QC Results`) |> 
          mutate(rna_isolation = "Duke"))

fastqs <- gcs_list_objects(bucket = "pathos-data-internal",
                           prefix = "duke-poci-prostate/poci_combination_cell_lines/01.RawData") |>
  dplyr::filter(stringr::str_detect(name, "fq.gz")) |>
  mutate(fastq = case_when(
    grepl("_1.fq",name) ~ 'fastq_1',
    grepl("_2.fq",name) ~ 'fastq_2',
    TRUE ~ NA
  )) |>
  mutate(sample = sub(".*/", "", name),
         sample = sub("_[1,2].fq.gz", "", sample),
         strandedness = 'auto',
         name = paste0(gs_path, name)) |>
  dplyr::select(name, sample, strandedness, fastq) |>
  tidyr::pivot_wider(names_from = fastq, values_from = name) |>
  dplyr::select(sample, fastq_1, fastq_2, strandedness) |>
  # add in details from the experimental design
  left_join(experimental_design, by = "sample")

write.csv(fastqs, file = "~/Desktop/duke_pocenbrodib_cell_line_combinations.csv", row.names = F)

