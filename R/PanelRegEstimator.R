#' Estimate a panel regression model
#'
#' This function estimates a panel regression model.
#'
#' @import data.table
#' @import fixest
#'
#' @param panel_data (data.table) A data.table containing the panel data. It should be in "long" format, which means that each row is uniquely identified by a (id, time) combination.
#' @param panel_model (character) The panel model to estimate. Must be one of "exogenous", "iid", "FE", "MA1", or "AR1". See the package vignette for details.
#' @param varnames (list) A list containing the names of the variables in the panel_data. It must contain the following elements: id_name, time_name, outcome_name, and endogenous_names. Optionally, it also contains covariate_names and fixedeffect_names. The id_name and time_name are the names of the unit and time variables, respectively. The outcome_name is the name of the outcome variable. The endogenous_names is a character vector of the names of the endogenous variables. The covariate_names is a character vector of the names of the covariates to control linearly, and the fixedeffect_names is a character vector of the names of the covariates to control using discrete indicators. If covariate_names is NULL, then the model is estimated without linear covariates. If fixedeffect_names is NULL, then the model is estimated without fixed effects. 
#' @param AR1_options (list) If panel_model != "AR1", then this list is ignored. If panel_model = "AR1", then this list must be non-empty. It contains the options for estimating the AR(1) model. It must contain the following elements: AR1_method, AR1_IV_lags, AR1_IV_outcome, and AR1_persistence. The AR1_method is the method to use to estimate the AR(1) model, and it must be one of "PanelIV" and "GMM". The AR1_IV_lags is the number of lags to use in the IV regression, which must be 1 or 2, and it is only used if AR1_method = "GMM". The AR1_IV_outcome specifies whether or not to include the lagged outcome variable as an instrument. The particular lags of the outcome used are inferred from AR1_IV_lags. The AR1_persistence is the value to force the AR(1) persistence parameter to equal; if it is NULL, then the AR(1) persistence parameter is estimated. The default is AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE).
#'
#' @return (data.table) The estimates from the panel regression model.
#' @export
PanelReg <- function(panel_data, panel_model, varnames, AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)) {

    ##########################################
    # 1. Input Validation
    ##########################################

    # unpack the variable names
    id_name = varnames$id_name
    time_name = varnames$time_name
    outcome_name = varnames$outcome_name
    endogenous_names = varnames$endogenous_names
    covariate_names = varnames$covariate_names
    fixedeffect_names = varnames$fixedeffect_names

    # unpack the AR1 options
    AR1_method = AR1_options$AR1_method
    AR1_IV_lags = AR1_options$AR1_IV_lags
    AR1_IV_outcome = AR1_options$AR1_IV_outcome
    AR1_persistence = AR1_options$AR1_persistence

    # validate the inputs
    PR.est.validate_inputs(panel_data, panel_model, id_name, time_name, outcome_name, endogenous_names, covariate_names, AR1_method, AR1_IV_lags, AR1_IV_outcome, AR1_persistence)

    # make sure the data is sorted
    panel_data = panel_data[order(get(id_name), get(time_name))]

    ##########################################
    # 2. Build Lagged Variables
    ##########################################

    combined_names = c(outcome_name, endogenous_names, covariate_names)

    # build lagged observable variables
    panel_data_lag = panel_data[, .SD, .SDcols = c(id_name, time_name, combined_names)]
    panel_data_lag[, (time_name) := get(time_name) + 1]
    setnames(panel_data_lag, combined_names, paste0(paste0(combined_names, "_lag")))

    # build twice lagged observable variables
    panel_data_lag2 = panel_data[, .SD, .SDcols = c(id_name, time_name, combined_names)]
    panel_data_lag2[, (time_name) := get(time_name) + 2]
    setnames(panel_data_lag2, combined_names, paste0(paste0(combined_names, "_lag2")))

    # merge the lagged variables into the panel_data
    panel_data = merge(panel_data, panel_data_lag, by = c(id_name, time_name))
    panel_data = merge(panel_data, panel_data_lag2, by = c(id_name, time_name))

    # if the panel model is FE, we need differenced variables
    if (panel_model == "FE") {
        # define differenced outcome and endogenous variables
        panel_data[, (paste0(outcome_name, "_diff")) := get(outcome_name) - get(paste0(outcome_name, "_lag"))]
        for (ii in seq_len(length(endogenous_names))) {
            panel_data[, paste0(endogenous_names[ii], "_diff") := get(endogenous_names[ii]) - get(paste0(endogenous_names[ii], "_lag"))]
        }
        if (length(covariate_names) > 0) {
            for (ii in seq_len(length(covariate_names))) {
                panel_data[, paste0(covariate_names[ii], "_diff") := get(covariate_names[ii]) - get(paste0(covariate_names[ii], "_lag"))]
            }
        }
    }


    ##########################################
    # 3. Estimate the Model
    ##########################################

    # estimate the exogenous model
    if (panel_model == "exogenous") {
        # build the formula
        the_formula = PR.formula.exogenous(outcome_name, endogenous_names, covariate_names, fixedeffect_names)
        # estimate the model
        est = feols(as.formula(the_formula), data = panel_data)
        # return the estimated second stage
        secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
    }

    # estimate the iid model
    if (panel_model == "iid") {
        # build the formula
        the_formula = PR.formula.iid(outcome_name, endogenous_names, covariate_names, fixedeffect_names)
        # estimate the model
        est = feols(as.formula(the_formula), data = panel_data)
        # return the estimated second stage
        secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
        secondstage[, Variable := gsub("fit_", "", Variable)]
    }

    # estimate the fixed effects FE model
    if (panel_model == "FE") {
        # build the formula
        the_formula = PR.formula.FE(outcome_name, endogenous_names, covariate_names, fixedeffect_names)
        # estimate the model
        est = feols(as.formula(the_formula), data = panel_data)
        # return the estimated second stage
        secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
        secondstage[, Variable := gsub("_diff", "", Variable)]
    }

    # estimate the MA1 model
    if (panel_model == "MA1") {
        # build the formula
        the_formula = PR.formula.MA1(outcome_name, endogenous_names, covariate_names, fixedeffect_names)
        # estimate the model
        est = feols(as.formula(the_formula), data = panel_data)
        # return the estimated second stage
        secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
        secondstage[, Variable := gsub("fit_", "", Variable)]
    }


    # estimate the AR1 model
    if (panel_model == "AR1") {
        if (AR1_method == "PanelIV") {
            secondstage = PR.est.AR1.PanelIV(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_outcome = AR1_IV_outcome, AR1_persistence = AR1_persistence)
        }
        if (AR1_method == "GMM") {
            return(PR.est.AR1.GMM(panel_data, outcome_name, endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, AR1_persistence = AR1_persistence))
        }
    }

    # return the second stage estimates
    setnames(secondstage, c("Variable", "Estimate", "SE", "tvalue", "pvalue"))
    secondstage = secondstage[Variable %in% c(endogenous_names, "AR1_persistence")]
    for (vv in c("Estimate", "SE", "tvalue", "pvalue")) {
        secondstage[, (vv) := round(get(vv), 6)]
    }
    return(secondstage)

}
