
#' Internal function to validate the inputs to PR.est.
#' @noRd
PR.est.validate_inputs <- function(panel_data, panel_model, id_name, time_name, outcome_name, endogenous_names, covariate_names, AR1_method, AR1_IV_lags, AR1_IV_outcome, AR1_persistence) {

    # check that the panel_data is a data.table
    if (!is.data.table(panel_data)) {
        stop("panel_data must be a data.table")
    }

    # check that the model is a character
    if (!is.character(panel_model) || length(panel_model) != 1 || !panel_model %in% c("exogenous", "iid", "FE", "MA1", "AR1")) {
        stop("panel_model must be 'exogenous', 'iid', 'FE', 'MA1', or 'AR1'")
    }

    # check that the id_name is a character
    if (!is.character(id_name) || length(id_name) != 1) {
        stop("id_name must be a character")
    }

    # check that the time_name is a character
    if (!is.character(time_name) || length(time_name) != 1) {
        stop("time_name must be a character")
    }

    # check that the outcome_name is a character
    if (!is.character(outcome_name) || length(outcome_name) != 1) {
        stop("outcome_name must be a character")
    }

    # check that the endogenous_names is a character
    if (!is.character(endogenous_names) || length(endogenous_names) == 0) {
        stop("endogenous_names must be a character vector of length > 0")
    }

    # check that the covariate_names is a character or NULL
    if (!is.null(covariate_names) && (!is.character(covariate_names) || length(covariate_names) == 0)) {
        stop("covariate_names must be a character vector or NULL")
    }

    # check that the outcome_name is in the panel_data
    if (!outcome_name %in% names(panel_data)) {
        stop("outcome_name not found in panel_data")
    }

    # check that unit_name and time_name are in the panel_data
    if (!all(c(id_name, time_name) %in% names(panel_data))) {
        stop("id_name and time_name must be in panel_data")
    }

    # check that the endogenous_names are in the panel_data
    if (!all(endogenous_names %in% names(panel_data))) {
        stop("not all endogenous_names found in panel_data")
    }

    # check that the covariate_names are in the panel_data if they are not NULL
    if (!is.null(covariate_names) && !all(covariate_names %in% names(panel_data))) {
        stop("not all covariate_names found in panel_data")
    }

    # check that time_name is integer-like
    if (!is.numeric(panel_data[, get(time_name)]) || any(panel_data[, get(time_name)] %% 1 != 0)) {
        stop("time_name must be an integer")
    }

    # check that time_name is sequential
    if (any(diff(sort(unique(panel_data[, get(time_name)]))) != 1)) {
        stop("time_name must be sequential")
    }

    # check that the panel_data has at least 2 time periods
    if (length(unique(panel_data[, get(time_name)])) < 3) {
        stop("panel_data must have at least 2 time periods")
    }

    # check that the panel_data has at least 2 id_name values
    if (length(unique(panel_data[, get(id_name)])) < 2) {
        stop("panel_data must have at least 2 units")
    }

    # check that the panel_data has no missing values
    if (any(is.na(panel_data))) {
        stop("panel_data must have no missing values")
    }

    # check that the panel_data has no infinite values
    for(vv in c(outcome_name, endogenous_names, covariate_names)) {
        if (any(!is.finite(panel_data[, get(vv)]))) {
            stop("panel_data must have no infinite values")
        }
    }

    # check that the panel_data has no duplicate rows
    if (nrow(panel_data) != nrow(unique(panel_data))) {
        stop("panel_data must have no duplicate rows")
    }

    # check that the panel_data has no duplicate unit_name-time_name combinations
    if (nrow(panel_data) != nrow(unique(panel_data[, list(get(id_name), get(time_name))]))) {
        stop("panel_data must have no duplicate id_name-time_name combinations")
    }

    # check that the AR1_method is a character in "GMM" or "PanelIV"
    if (panel_model == "AR1" && (!is.character(AR1_method) || length(AR1_method) != 1 || !AR1_method %in% c("GMM", "PanelIV"))) {
        stop("AR1_method must be 'GMM' or 'PanelIV'")
    }

    # check that the AR1_IV_lags is an integer in 1 or 2, only if AR1_method = "GMM"
    if (panel_model == "AR1" && AR1_method == "GMM" && (!is.numeric(AR1_IV_lags) || length(AR1_IV_lags) != 1 || !AR1_IV_lags %in% c(1, 2))) {
        stop("AR1_IV_lags must be 1 or 2")
    }

    # check that the AR1_IV_outcome is a logical
    if (panel_model == "AR1" && (!is.logical(AR1_IV_outcome) || length(AR1_IV_outcome) != 1)) {
        stop("AR1_IV_outcome must be TRUE or FALSE")
    }

    # check that the AR1_persistence is a numeric or NULL
    if (panel_model == "AR1" && !is.null(AR1_persistence) && (!is.numeric(AR1_persistence) || length(AR1_persistence) != 1)) {
        stop("AR1_persistence must be a numeric or NULL")
    }

    # if method is "GMM" and AR1_IV_lags = 1, then AR1_IV_outcome must be TRUE
    if (panel_model == "AR1" && AR1_method == "GMM" && AR1_IV_lags == 1 && !AR1_IV_outcome) {
        stop("if AR1_method = 'GMM' and AR1_IV_lags = 1, then AR1_IV_outcome must be TRUE")
    }

    return(TRUE)

}
