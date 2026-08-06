#' Calculate Best Overall RECIST Response per Patient per Regimen
#'
#' Computes the best overall response (BOR) using RECIST priority
#' (CR > PR > SD > PD) and the last recorded response relative to treatment
#' start, for every patient × regimen pair in the supplied data.
#'
#' @param data A file path to a Tempus cohort directory (loaded via
#'   \code{tempusr::load_tempus_data}), a data.frame (treated as
#'   onco_response_assessment directly), or a named list containing at minimum
#'   \code{onco_response_assessment} and optionally \code{onco_regimen}
#'   (used to derive days_from_tx_start).
#' @param regimen_data Optional onco_regimen data.frame. Required for
#'   days_from_tx_start when \code{data} is a plain data.frame. Ignored when
#'   \code{data} is a list that already contains onco_regimen.
#' @param response_col Column in onco_response_assessment holding response
#'   values (default "assessment_value_rollup").
#' @param date_col Column in onco_response_assessment holding assessment dates
#'   (default "assessment_date_indexed"). Used to compute days_from_tx_start
#'   and to identify the last recorded response.
#' @param agent_filter Optional filter restricting which regimens are included.
#'   Two forms are accepted:
#'   \itemize{
#'     \item A character vector of \code{regimen_id} values — rows in
#'           onco_response_assessment are kept by exact match.
#'     \item A single grep-compatible search string (e.g.
#'           \code{"tisotumab|enfortumab"}) applied case-insensitively to the
#'           agent name column of onco_regimen (\code{agents} in DM2,
#'           \code{regimen_name} in DM1). The matched regimen_ids are then used
#'           to filter onco_response_assessment. \strong{onco_regimen must be
#'           available} (via the \code{data} list or \code{regimen_data}) when
#'           this form is used.
#'   }
#'
#' @return A data.frame with one row per patient_id × regimen_id containing:
#'   \itemize{
#'     \item \code{patient_id}
#'     \item \code{regimen_id}
#'     \item \code{response} — best overall response abbreviated as CR, PR,
#'           SD, or PD (NA for unrecognised values)
#'     \item \code{assessment_method} — method used at the BOR visit
#'     \item \code{days_from_tx_start} — days from regimen start to the BOR
#'           assessment (NA if date or regimen start unavailable)
#'     \item \code{response_at_tx_end} — response recorded at the latest
#'           assessment with a valid date
#'     \item \code{response_at_tx_end_days_from_tx_start} — days from regimen
#'           start to that last assessment
#'   }

calc_BOR <- function(data,
                     regimen_data = NULL,
                     response_col = "assessment_value_rollup",
                     date_col     = "assessment_date_indexed",
                     agent_filter = NULL) {

    library(dplyr)

    # ── 1. Extract tables ──────────────────────────────────────────────────────
    if (is.character(data)) {
        td <- tempusr::load_tempus_data(
            path_main  = data,
            collection = NULL,
            list_files = c("onco_response_assessment", "onco_regimen")
        )
        onco_response_assessment <- td$onco_response_assessment
        if (is.null(onco_response_assessment))
            stop("onco_response_assessment could not be loaded from: ", data)
        onco_regimen <- if (!is.null(regimen_data)) regimen_data else td$onco_regimen
    } else if (is.data.frame(data)) {
        onco_response_assessment <- data
        onco_regimen <- regimen_data
    } else if (is.list(data)) {
        onco_response_assessment <- data$onco_response_assessment
        if (is.null(onco_response_assessment))
            stop("onco_response_assessment not found in data list")
        onco_regimen <- if (!is.null(regimen_data)) regimen_data else data$onco_regimen
    } else {
        stop("data must be a file path, a data.frame, or a named list containing onco_response_assessment")
    }

    # ── 2. Apply agent_filter ─────────────────────────────────────────────────
    if (!is.null(agent_filter)) {
        # Determine mode: regimen_id vector vs. grep search term.
        # If ANY supplied value matches a regimen_id in the response table,
        # treat the whole vector as regimen_ids; otherwise treat as grep pattern.
        is_id_filter <- any(agent_filter %in% onco_response_assessment$regimen_id)

        if (is_id_filter) {
            onco_response_assessment <- dplyr::filter(
                onco_response_assessment,
                regimen_id %in% agent_filter
            )
        } else {
            # Grep mode — requires onco_regimen
            if (is.null(onco_regimen))
                stop(paste(
                    "agent_filter appears to be a search term but onco_regimen is not available.",
                    "Supply it via the data list or regimen_data."
                ))

            agent_col <- if ("agents" %in% colnames(onco_regimen)) "agents"
                         else if ("regimen_name" %in% colnames(onco_regimen)) "regimen_name"
                         else stop("Cannot find agent name column (agents / regimen_name) in onco_regimen")

            matched_ids <- onco_regimen %>%
                dplyr::filter(grepl(paste(agent_filter, collapse = "|"),
                                    .data[[agent_col]],
                                    ignore.case = TRUE)) %>%
                dplyr::pull(regimen_id) %>%
                unique()

            if (length(matched_ids) == 0)
                stop(sprintf(
                    "agent_filter '%s' matched no regimens in onco_regimen",
                    paste(agent_filter, collapse = "|")
                ))

            onco_response_assessment <- dplyr::filter(
                onco_response_assessment,
                regimen_id %in% matched_ids
            )
        }

        if (nrow(onco_response_assessment) == 0)
            stop("agent_filter matched no rows in onco_response_assessment")
    }

    # ── 3. Attach days_from_tx_start via onco_regimen ─────────────────────────
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
        if (!has_date)
            warning(sprintf(
                "date_col '%s' not found in response data; days_from_tx_start will be NA",
                date_col
            ))
        if (!has_regimen)
            warning(paste(
                "onco_regimen not supplied or missing start_date_indexed;",
                "days_from_tx_start will be NA"
            ))
        resp <- dplyr::mutate(resp, days_from_tx_start = NA_real_)
    }

    # ── 4. RECIST priority map ─────────────────────────────────────────────────
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

    # ── 5. Best overall response (BOR) ────────────────────────────────────────
    bor <- resp_valid %>%
        dplyr::group_by(patient_id, regimen_id) %>%
        dplyr::slice_min(Priority, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(patient_id, regimen_id,
                      response, assessment_method, days_from_tx_start)

    # ── 6. Last response by assessment date ───────────────────────────────────
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

    # ── 7. Combine and abbreviate response labels ─────────────────────────────
    recist_abbrev <- c(
        "Complete Response"   = "CR",
        "Partial Response"    = "PR",
        "Stable Disease"      = "SD",
        "Progressive Disease" = "PD"
    )

    dplyr::left_join(bor, last_resp, by = c("patient_id", "regimen_id")) %>%
        dplyr::mutate(
            response = dplyr::recode(response, !!!recist_abbrev, .default = NA_character_),
            response_at_tx_end = dplyr::recode(response_at_tx_end,
                                                !!!recist_abbrev,
                                                .default = NA_character_)
        )
}
