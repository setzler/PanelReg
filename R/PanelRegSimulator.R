
#' Simulate a panel dataset.
#'
#' This function simulates a panel dataset with a given number of units and time periods. The user can specify the distribution of the error term, the number of units, the number of time periods, and the parameters of the model. The function returns a panel dataset in the form of a data.table.
#'
#' @import data.table
#' @param seed (integer) A seed for the random number generator.
#' @param sample_size (integer) The number of units in the panel dataset.
#' @param min_year (integer) The first year in the panel dataset.
#' @param max_year (integer) The last year in the panel dataset.
#' @param noise_sd (numeric) The standard deviation of the measurement error.
#' @param panel_model (character) The panel model to use. Options are "iid", "FE", "AR1", "MA1".
#' @param true_beta (numeric) The true value of the parameter of interest.
#' @param endog_weight (numeric) The weight of the endogenous variable in the panel model. Must be between 0 and 1.
#' @param MA1_persistence (numeric) The persistence parameter of the MA(1) process, i.e., eps_{it} = shock_{it} + MA1_persistence * shock_{it-1}.
#' @param AR1_persistence (numeric) The persistence parameter of the AR(1) process, i.e., eps_{it} = AR1_persistence * eps_{it-1} + shock_{it}.
#' @param true_gamma (numeric) The true value of the parameter for the exogenous variables. The first entry corresponds to the intercept. By default, true_gamma = c(1).
#' @param return_unobservables (logical) If TRUE, the function returns the full panel dataset, including the unobservable variables. If FALSE, the function returns only the observable variables.
#' @return (data.table) A panel dataset.
#' @export
PanelRegSim <- function(seed = 1234, sample_size = 100, min_year = 2000, max_year = 2010, noise_sd = 1, panel_model = "iid", true_beta = 0.5, endog_weight = 0.5, MA1_persistence = 0.5, AR1_persistence = 0.5, true_gamma = c(1), return_unobservables = FALSE) {


    ###############################################
    # 1. Input Validation
    ###############################################

    # validate the inputs
    PR.sim.validate_inputs(seed, sample_size, min_year, max_year, noise_sd, panel_model, true_beta, endog_weight, MA1_persistence, AR1_persistence, return_unobservables)

    # set the seed
    set.seed(seed)

    ###############################################
    # 2. Create the Panel Data Skeleton
    ###############################################

    # construct a blank dataset with the appropriate identifier combinations
    paneldata = data.table(expand.grid(unit_id = 1:sample_size, time_id = min_year:max_year))

    # draw the noise variable (eta_it)
    paneldata[, measurement_error := rnorm(.N, sd = noise_sd)]
    paneldata[, measurement_error := measurement_error - mean(measurement_error), by = "time_id"]

    # define the shock
    paneldata[, shock := rnorm(.N, sd = 1)]
    paneldata[, shock := shock - mean(shock), by = "time_id"]

    # make sure the data is sorted
    paneldata = paneldata[order(unit_id, time_id)]

    ###############################################
    # 3. Draw the Endogenous Error Process
    ###############################################


    # define panel_models of epsilon process (eps_it)
    if (panel_model %in% c("exogenous", "iid")) {
        paneldata[, eps_it := copy(shock)]
    }
    if (panel_model == "FE") {
        paneldata[, eps_it := rnorm(n = 1, sd = 1), by = "unit_id"]
    }
    if (panel_model == "MA1") {
        # merge in the lagged shock and set initial lag to 0
        paneldata = merge(paneldata, paneldata[, list(unit_id, time_id = time_id + 1, shock_lag = shock)], by = c("unit_id", "time_id"), all.x = TRUE)
        paneldata[is.na(shock_lag), shock_lag := 0.0]
        # define the eps_it as an MA(1) process
        paneldata[, eps_it := shock + MA1_persistence * shock_lag]
    }
    if (panel_model == "AR1") {
        # initialize the eps_it
        paneldata[, eps_it := 0.0]
        paneldata[time_id == min_year, eps_it := shock/ (1 - AR1_persistence^2)]
        for(tt in (min_year + 1):(max_year)) {
            # merge in the lagged eps_it
            paneldata_lag = paneldata[time_id == (tt - 1), list(unit_id, eps_it_lag = eps_it)]
            paneldata = merge(paneldata, paneldata_lag, by = "unit_id", all.x = TRUE)
            # update the autoregression
            paneldata[time_id == tt, eps_it := AR1_persistence * eps_it_lag + shock]
            # clean up the lagged variable
            paneldata[, eps_it_lag := NULL]
            paneldata_lag  = NULL
        }
    }

    ###############################################
    # 4. Construct Endogenous Variables and Outcome
    ###############################################

    num_endog_vars = length(true_beta)

    # draw a random walk to force autocorrelation in each endogenous variable
    for (ii in 1:num_endog_vars) {
        paneldata[, paste0("random_walk", ii) := cumsum(rnorm(.N, sd = 1)), by = "unit_id"]
    }

    # set up endogenity weights
    endog_weights = endog_weight * (num_endog_vars:1 / num_endog_vars)
    if (panel_model == "exogenous") {
        endog_weights = endog_weights * 0.0
    }

    # construct the endogenous variables
    for (ii in 1:num_endog_vars) {
        paneldata[, paste0("endog_var", ii) := 3 + endog_weights[ii] * eps_it + (1 - endog_weights[ii]) * rnorm(n = .N, sd = 1) + get(paste0("random_walk", ii))]
    }

    # define the dependent variable
    paneldata[, outcome := eps_it + measurement_error]
    for (ii in 1:num_endog_vars) {
        paneldata[, outcome := outcome + true_beta[ii] * get(paste0("endog_var", ii))]
    }
    num_covariates = length(true_gamma)
    if (num_covariates > 0) {
        for (ii in 1:num_covariates) {
            if (ii == 1) {
                paneldata[, outcome := outcome + true_gamma[ii]] # add an intercept
            }
            if (ii > 1) {
                paneldata[, (paste0("co_var", ii - 1)) := 3 + rnorm(n = .N, sd = 1) + 0.3 * get(paste0("endog_var", 1))]
                paneldata[, outcome := outcome + true_gamma[ii] * get(paste0("co_var", ii - 1))]
            }
        }
    }

    ###############################################
    # 5. Return the Panel Data
    ###############################################

    if (return_unobservables) {
        return(paneldata)
    }
    keep_vars = c("unit_id", "time_id", "outcome", paste0("endog_var", 1:num_endog_vars))
    if (num_covariates > 1) {
        keep_vars = c(keep_vars, paste0("co_var", 1:(num_covariates - 1)))
    }
    return(paneldata[, .SD, .SDcols = keep_vars])
}


#' Internal function to validate the inputs to PR.sim.
#' @noRd
PR.sim.validate_inputs <- function(seed, sample_size, min_year, max_year, noise_sd, panel_model, true_beta, endog_weight, MA1_persistence, AR1_persistence, return_unobservables) {

    # Check that inputs are valid
    if (!is.numeric(seed) || length(seed) != 1 || seed <= 0 || seed %% 1 != 0) {
        stop("seed must be a positive integer")
    }
    if (!is.numeric(sample_size) || length(sample_size) != 1 || sample_size <= 0 || sample_size %% 1 != 0) {
        stop("sample_size must be a positive integer")
    }
    if (!is.numeric(min_year) || length(min_year) != 1 || min_year %% 1 != 0) {
        stop("min_year must be an integer")
    }
    if (!is.numeric(max_year) || length(max_year) != 1 || max_year %% 1 != 0) {
        stop("max_year must be an integer")
    }
    if (min_year >= max_year) {
        stop("max_year must be greater than min_year")
    }
    if (!is.numeric(noise_sd) || length(noise_sd) != 1 || noise_sd < 0) {
        stop("noise_sd must be a non-negative number")
    }
    if (!is.character(panel_model) || length(panel_model) != 1 || !panel_model %in% c("exogenous", "iid", "FE", "MA1", "AR1")) {
        stop("panel_model must be 'exogenous', 'iid', 'FE', 'MA1', or 'AR1'")
    }
    if (!is.numeric(true_beta)) {
        stop("true_beta must be a number")
    }
    if (!is.numeric(endog_weight) || length(endog_weight) != 1 || endog_weight < 0 || endog_weight > 1) {
        stop("endog_weight must be a number between 0 and 1")
    }

}

#' @importFrom stats arima.sim ecdf lm na.omit nobs rnorm runif var vcov as.formula optim
#' @importFrom data.table .N .SD copy data.table rbindlist setDT setnames `:=`
NULL

# truly ridiculous, variables in data.table objects are forced to be globals
Variable <- eps_it <- eps_it_lag <- measurement_error <- outcome <- shock <- shock_lag <- time_id <- u_it <- unit_id <- NULL
