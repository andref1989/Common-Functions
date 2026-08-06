#' Calculate Objective Response Rate (ORR) by Drug Regimen and Line of Therapy
#'
#' Computes ORR (and clinical benefit rate) per regimen × LOT, restricted to
#' regimens used in at least `min_patients` unique patients AND in more than
#' `min_pct` of all patients with treatment data at that line of therapy.
#' Best overall RECIST response (CR > PR > SD > PD) is used per patient per regimen.
#'
#' @param data_cohort File path to a Tempus cohort directory, OR a pre-loaded
#'   list from load_tempus_data() containing onco_regimen. If the list also
#'   contains onco_response_assessment it will be used directly; otherwise the
#'   function will attempt to load it from the cohort path.
#' @param response_data Optional pre-loaded onco_response_assessment data frame.
#'   Required when data_cohort is a list that does not contain it.
#' @param min_patients Minimum unique patients on a regimen × LOT (default 10).
#' @param min_pct Minimum fraction of patients at that LOT who received the
#'   regimen (default 0.10 = 10%).
#' @param max_lot Maximum line of therapy to include (default 6).
#' @param response_col Column in onco_response_assessment with response values
#'   (default "assessment_value_rollup").
#' @param agent_filter Optional regex string applied to the agent/drug name
#'   column (agents in DM2, regimen_name in DM1). Keeps only regimens whose
#'   name matches the pattern. Case-insensitive. E.g. "pembrolizumab|nivolumab".
#' @param therapy_class_filter Optional regex string applied to the therapy
#'   class column (therapy_class in DM2, regimen_class in DM1). E.g.
#'   "Immunotherapy|Anti-PD".
#' @param therapy_class_group_filter Optional regex string applied to the
#'   therapy class group column (therapy_class_group in DM2,
#'   regimen_class_group in DM1). E.g. "Biologic|Targeted".
#' @param assessment_method_filter Optional character vector of allowed
#'   assessment methods (e.g. c("Imaging", "MD Dictated")). When NULL all
#'   methods are included.
#' @param lot_col Optional column name to use as line of therapy. When NULL
#'   the column is chosen automatically: "line_of_therapy_number" for DM2,
#'   "regimen_rank" for DM1.
#' @param patient_filter Optional character vector of patient_ids. When
#'   provided, analysis is restricted to those patients only.
#'
#' @return A data frame sorted by LOT then descending ORR, with columns:
#'   Regimen, LOT, N_patients, N_with_response, N_responders,
#'   ORR (%), CBR (%), N_CR, N_PR, N_SD, N_PD.

calc_orr_by_regimen <- function(data_cohort,
                                 response_data              = NULL,
                                 min_patients               = 10,
                                 min_pct                    = 0.10,
                                 max_lot                    = 6,
                                 response_col               = "assessment_value_rollup",
                                 agent_filter               = NULL,
                                 therapy_class_filter       = NULL,
                                 therapy_class_group_filter = NULL,
                                 assessment_method_filter   = NULL,
                                 lot_col                    = NULL,
                                 patient_filter             = NULL) {

    library(dplyr)

    # ── Helper: best overall RECIST response ───────────────────────────────────
    calc_BOR <- function(data,
                         regimen_data = NULL,
                         response_col = "assessment_value_rollup",
                         date_col     = "assessment_date_indexed") {

        if (is.data.frame(data)) {
            onco_response_assessment <- data
            onco_regimen <- regimen_data
        } else if (is.list(data)) {
            onco_response_assessment <- data$onco_response_assessment
            if (is.null(onco_response_assessment))
                stop("onco_response_assessment not found in data list")
            onco_regimen <- if (!is.null(regimen_data)) regimen_data else data$onco_regimen
        } else {
            stop("data must be a data.frame or a named list containing onco_response_assessment")
        }

        has_date    <- date_col %in% colnames(onco_response_assessment)
        has_regimen <- !is.null(onco_regimen) &&
                       "start_date_indexed" %in% colnames(onco_regimen)

        resp <- onco_response_assessment %>%
            dplyr::rename(response = dplyr::all_of(response_col))

        if (has_date && has_regimen) {
            resp <- resp %>%
                dplyr::left_join(
                    dplyr::select(onco_regimen, regimen_id, start_date_indexed),
                    by = "regimen_id"
                ) %>%
                dplyr::mutate(
                    days_from_tx_start = as.numeric(
                        as.Date(.data[[date_col]]) - as.Date(start_date_indexed)
                    )
                )
        } else {
            resp <- dplyr::mutate(resp, days_from_tx_start = NA_real_)
        }

        recist_priority <- c(
            "Complete Response"   = 1L,
            "Partial Response"    = 2L,
            "Stable Disease"      = 3L,
            "Progressive Disease" = 4L
        )

        resp_valid <- resp %>%
            dplyr::filter(!is.na(response), response != "") %>%
            dplyr::mutate(Priority = dplyr::recode(response,
                                                    !!!recist_priority,
                                                    .default = 5L))

        bor <- resp_valid %>%
            dplyr::group_by(patient_id, regimen_id) %>%
            dplyr::slice_min(Priority, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup() %>%
            dplyr::select(patient_id, regimen_id,
                          response, assessment_method, days_from_tx_start)

        if (has_date) {
            last_resp <- resp_valid %>%
                dplyr::filter(!is.na(.data[[date_col]])) %>%
                dplyr::mutate(.sort_date = as.Date(.data[[date_col]])) %>%
                dplyr::group_by(patient_id, regimen_id) %>%
                dplyr::slice_max(.sort_date, n = 1, with_ties = FALSE) %>%
                dplyr::ungroup() %>%
                dplyr::select(patient_id, regimen_id,
                              response_at_tx_end = response,
                              response_at_tx_end_days_from_tx_start = days_from_tx_start)
        } else {
            last_resp <- dplyr::select(bor, patient_id, regimen_id) %>%
                dplyr::mutate(response_at_tx_end                    = NA_character_,
                              response_at_tx_end_days_from_tx_start = NA_real_)
        }

        dplyr::left_join(bor, last_resp, by = c("patient_id", "regimen_id"))
    }

    # ── 1. Load / extract tables ───────────────────────────────────────────────
    if (is.character(data_cohort)) {

        # tempusr::load_tempus_data handles both CSV and parquet delivery formats
        td <- tempusr::load_tempus_data(
            path_main  = data_cohort,
            collection = NULL,
            list_files = c("onco_regimen", "onco_response_assessment")
        )

        onco_regimen <- td$onco_regimen
        if (is.null(onco_regimen))
            stop("onco_regimen could not be loaded from: ", data_cohort)

        onco_response_assessment <- td$onco_response_assessment
        if (is.null(onco_response_assessment))
            stop("onco_response_assessment could not be loaded from: ", data_cohort)

    } else if (is.list(data_cohort)) {

        onco_regimen <- data_cohort$onco_regimen
        if (is.null(onco_regimen))
            stop("onco_regimen not found in data_cohort list")

        if (!is.null(response_data)) {
            onco_response_assessment <- response_data
        } else if (!is.null(data_cohort$onco_response_assessment)) {
            onco_response_assessment <- data_cohort$onco_response_assessment
        } else {
            stop("Provide response_data or ensure data_cohort list contains onco_response_assessment")
        }

    } else {
        stop("data_cohort must be a file path or a list from load_tempus_data()")
    }

    # ── 2. Apply patient filter ────────────────────────────────────────────────
    if (!is.null(patient_filter)) {
        onco_regimen <- dplyr::filter(onco_regimen, patient_id %in% patient_filter)
        if (nrow(onco_regimen) == 0)
            stop("patient_filter matched no rows in onco_regimen")
        onco_response_assessment <- dplyr::filter(onco_response_assessment, patient_id %in% patient_filter)
    }

    # ── 3. Detect data model and standardise column names ──────────────────────
    if ("agents" %in% colnames(onco_regimen)) {
        # DM2
        dm <- "DM2"
        agent_col  <- "agents"
        class_col  <- if ("therapy_class"       %in% colnames(onco_regimen)) "therapy_class"       else NULL
        group_col  <- if ("therapy_class_group"  %in% colnames(onco_regimen)) "therapy_class_group"  else NULL
        lot_col    <- if (!is.null(lot_col)) lot_col else "line_of_therapy_number"
    } else if ("regimen_name" %in% colnames(onco_regimen)) {
        # DM1
        dm <- "DM1"
        agent_col  <- "regimen_name"
        class_col  <- if ("regimen_class"       %in% colnames(onco_regimen)) "regimen_class"       else NULL
        group_col  <- if ("regimen_class_group"  %in% colnames(onco_regimen)) "regimen_class_group"  else NULL
        lot_col    <- if (!is.null(lot_col)) lot_col else "regimen_rank"
    } else {
        stop("Cannot detect DM1 or DM2 column names in onco_regimen")
    }

    if (!lot_col %in% colnames(onco_regimen))
        stop(sprintf("lot_col '%s' not found in onco_regimen", lot_col))

    # ── 4. Apply modality / treatment type filters to the regimen table ────────
    if (!is.null(agent_filter)) {
        onco_regimen <- dplyr::filter(
            onco_regimen,
            grepl(agent_filter, .data[[agent_col]], ignore.case = TRUE)
        )
        if (nrow(onco_regimen) == 0)
            stop(sprintf("agent_filter '%s' matched no rows in onco_regimen", agent_filter))
    }

    if (!is.null(therapy_class_filter)) {
        if (is.null(class_col))
            warning("therapy_class_filter ignored: no therapy class column found")
        else {
            onco_regimen <- dplyr::filter(
                onco_regimen,
                grepl(therapy_class_filter, .data[[class_col]], ignore.case = TRUE)
            )
            if (nrow(onco_regimen) == 0)
                stop(sprintf("therapy_class_filter '%s' matched no rows", therapy_class_filter))
        }
    }

    if (!is.null(therapy_class_group_filter)) {
        if (is.null(group_col))
            warning("therapy_class_group_filter ignored: no therapy class group column found")
        else {
            onco_regimen <- dplyr::filter(
                onco_regimen,
                grepl(therapy_class_group_filter, .data[[group_col]], ignore.case = TRUE)
            )
            if (nrow(onco_regimen) == 0)
                stop(sprintf("therapy_class_group_filter '%s' matched no rows", therapy_class_group_filter))
        }
    }

    # ── 5. Standardise regimen prep ────────────────────────────────────────────
    regimen_prep <- onco_regimen %>%
        dplyr::select(patient_id, regimen_id,
                      Regimen = dplyr::all_of(agent_col),
                      LOT     = dplyr::all_of(lot_col)) %>%
        dplyr::filter(!is.na(LOT), LOT >= 1, LOT <= max_lot) %>%
        dplyr::mutate(LOT_label = factor(paste0(LOT, "L"),
                                         levels = paste0(seq_len(max_lot), "L")))

    # ── 6. Apply usage thresholds ──────────────────────────────────────────────
    lot_totals <- regimen_prep %>%
        dplyr::group_by(LOT_label) %>%
        dplyr::summarize(N_lot_total = dplyr::n_distinct(patient_id), .groups = "drop")

    regimen_eligible <- regimen_prep %>%
        dplyr::group_by(LOT_label, Regimen) %>%
        dplyr::summarize(N_patients = dplyr::n_distinct(patient_id), .groups = "drop") %>%
        dplyr::left_join(lot_totals, by = "LOT_label") %>%
        dplyr::mutate(Pct_of_LOT = N_patients / N_lot_total) %>%
        dplyr::filter(N_patients >= min_patients, Pct_of_LOT > min_pct) %>%
        dplyr::select(LOT_label, Regimen)

    if (nrow(regimen_eligible) == 0) {
        warning("No regimens met the minimum patient / percentage thresholds.")
        return(data.frame())
    }

    # ── 7. Best RECIST response per patient per regimen ────────────────────────
    eligible_regimen_ids <- regimen_prep %>%
        dplyr::semi_join(regimen_eligible, by = c("LOT_label", "Regimen")) %>%
        dplyr::pull(regimen_id) %>%
        unique()

    resp_filtered <- onco_response_assessment %>%
        dplyr::filter(regimen_id %in% eligible_regimen_ids)

    if (!is.null(assessment_method_filter)) {
        resp_filtered <- dplyr::filter(
            resp_filtered,
            assessment_method %in% assessment_method_filter
        )
        if (nrow(resp_filtered) == 0)
            warning("assessment_method_filter matched no response rows; ORR will be NA for all regimens")
    }

    best_response <- calc_BOR(
        data         = list(onco_response_assessment = resp_filtered,
                            onco_regimen              = onco_regimen),
        response_col = response_col
    ) %>%
        dplyr::select(patient_id, regimen_id, Response = response)

    # ── 8. Join and calculate ORR ──────────────────────────────────────────────
    orr_df <- regimen_prep %>%
        dplyr::semi_join(regimen_eligible, by = c("LOT_label", "Regimen")) %>%
        dplyr::left_join(best_response, by = c("patient_id", "regimen_id")) %>%
        dplyr::group_by(LOT_label, Regimen) %>%
        dplyr::summarize(
            N_patients      = dplyr::n_distinct(patient_id),
            N_with_response = dplyr::n_distinct(patient_id[!is.na(Response)]),
            N_CR = dplyr::n_distinct(patient_id[Response == "Complete Response"]),
            N_PR = dplyr::n_distinct(patient_id[Response == "Partial Response"]),
            N_SD = dplyr::n_distinct(patient_id[Response == "Stable Disease"]),
            N_PD = dplyr::n_distinct(patient_id[Response == "Progressive Disease"]),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            N_responders = N_CR + N_PR,
            ORR = dplyr::if_else(N_with_response > 0,
                                  round(100 * N_responders / N_with_response, 1),
                                  NA_real_),
            CBR = dplyr::if_else(N_with_response > 0,
                                  round(100 * (N_CR + N_PR + N_SD) / N_with_response, 1),
                                  NA_real_)
        ) %>%
        dplyr::arrange(LOT_label, dplyr::desc(ORR)) %>%
        dplyr::select(Regimen, LOT = LOT_label,
                      N_patients, N_with_response, N_responders,
                      ORR, CBR, N_CR, N_PR, N_SD, N_PD)

    return(orr_df)
}
