#' Calculate Most Frequently Used Treatments by Year and Line of Therapy
#'
#' Takes a path or pre-loaded list object, detects whether a treatment regimen
#' table is present (DM1: regimen; DM2: onco_regimen), and returns the top
#' treatments ranked by patient count within each combination of treatment
#' initiation year and line of therapy.
#'
#' @param cohort Complete path to the Tempus patient cohort or list object
#'   created by load_tempus_data (string or list). Only the clinical collection
#'   is required.
#' @param treatment_lines Maximum number of treatment lines to include
#'   (Integer). Lines beyond this are dropped before summarising.
#' @param top_n Number of top treatments to return per year × line combination
#'   (Integer).
#' @param verbose Whether to print progress messages (Boolean).
#'
#' @return A named list with two data frames:
#'   \describe{
#'     \item{Drug}{Top drug regimen names by year and line of therapy. Columns:
#'       Year, Treatment_Line, Treatment, N, Pct.}
#'     \item{Class}{Top drug class groups by year and line of therapy. Columns:
#'       Year, Treatment_Line, Treatment, N, Pct.}
#'   }
#'   \code{N} is the number of patients receiving that treatment, and
#'   \code{Pct} is the percentage of patients at that year × line combination.
#'   Rows are ordered by Year, Treatment_Line, then N (descending).
#'
#' @note
#'   Year is derived from \code{regimen_start_date_year_indexed} (DM1) or
#'   \code{start_date_year_indexed} (DM2). These are indexed/de-identified
#'   years and may not reflect true calendar years depending on the cohort.
#'   Records with a missing Year are excluded from the summary.
#'
#' @note
#'   Any patient-level filtering should be applied to the regimen table before
#'   passing the cohort to this function.
#'
#' @examples
#' \dontrun{
#' result <- calc_top_treatments(cohort, treatment_lines = 4, top_n = 5)
#' result$Drug    # top drug names by year and line
#' result$Class   # top drug classes by year and line
#' }
#'
calc_top_treatments <- function(cohort,
                                treatment_lines = 5,
                                top_n = 5,
                                verbose = FALSE) {

  list_files_dm1 <- c("regimen")
  list_files_dm2 <- c("onco_regimen")

  # ── Load cohort ────────────────────────────────────────────────────────────
  if (is.character(cohort)) {
    if (verbose) message("Trying to load cohort from path")
    input_td <- tryCatch(
      load_tempus_data(cohort, collection = "clinical"),
      error = function(e) {
        tryCatch(
          load_tempus_data(cohort, collection = NULL, list_files = list_files_dm2),
          error = function(f) stop("Could not load a regimen table from the provided path")
        )
      }
    )
  } else if (is.list(cohort)) {
    if (verbose) message("Cohort already loaded, working")
    input_td <- cohort
  } else {
    stop("'cohort' must be a file path (character) or a pre-loaded list object")
  }

  # ── Detect data model ──────────────────────────────────────────────────────
  data_model_check <- tempusr:::calc_data_model(
    input_td,
    list_tables_dm_1 = list_files_dm1,
    list_tables_dm_2 = list_files_dm2
  )
  stopifnot(data_model_check %in% c("1.0", "2.0"))

  if (verbose) message("Detected data model: ", data_model_check)

  # ── DM-specific extraction ─────────────────────────────────────────────────
  prep_dm1 <- function(td) {
    td$regimen %>%
      dplyr::arrange(.data$patient_id, .data$regimen_rank) %>%
      dplyr::select(
        .data$patient_id,
        Rank       = .data$regimen_rank,
        Drug_Name  = .data$regimen_name,
        Class_Group = .data$regimen_class_group,
        Year       = .data$regimen_start_date_year_indexed
      )
  }

  prep_dm2 <- function(td) {
    td$onco_regimen %>%
      dplyr::arrange(.data$patient_id, .data$regimen_sequence) %>%
      dplyr::select(
        .data$patient_id,
        Rank        = .data$regimen_sequence,
        Drug_Name   = .data$agents,
        Class_Group = .data$therapy_class_group,
        Year        = .data$start_date_year_indexed
      )
  }

  regimen_data <- if (data_model_check == "1.0") prep_dm1(input_td) else prep_dm2(input_td)

  # ── Add sequential treatment line label and filter ─────────────────────────
  line_labels <- paste0(seq_len(treatment_lines), "L")

  regimen_data <- regimen_data %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::mutate(Treatment_Line = paste0(seq_along(.data$Rank), "L")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$Treatment_Line %in% line_labels) %>%
    dplyr::mutate(Treatment_Line = factor(.data$Treatment_Line, levels = line_labels))

  # ── Helper: count and rank treatments within Year × Treatment_Line ─────────
  summarise_top_n <- function(data, treatment_col) {
    data %>%
      dplyr::filter(!is.na(.data$Year)) %>%
      dplyr::rename(Treatment = dplyr::all_of(treatment_col)) %>%
      dplyr::filter(!is.na(.data$Treatment)) %>%
      dplyr::count(.data$Year, .data$Treatment_Line, .data$Treatment,
                   name = "N", .drop = FALSE) %>%
      dplyr::group_by(.data$Year, .data$Treatment_Line) %>%
      dplyr::mutate(
        Pct          = round(100 * .data$N / sum(.data$N), 1),
        Rank_in_group = dplyr::row_number(dplyr::desc(.data$N))
      ) %>%
      dplyr::filter(.data$Rank_in_group <= top_n) %>%
      dplyr::arrange(.data$Year, .data$Treatment_Line, .data$Rank_in_group) %>%
      dplyr::select(-.data$Rank_in_group) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }

  drug_summary  <- summarise_top_n(regimen_data, "Drug_Name")
  class_summary <- summarise_top_n(regimen_data, "Class_Group")

  list(Drug = drug_summary, Class = class_summary)
}
