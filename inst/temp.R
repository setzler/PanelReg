
# simulate example data
simulated_data = PanelRegSim(sample_size = 1e5, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), panel_model = "AR1", noise_sd = 0, seed = 1)

# define varnames
varnames = list(
    id_name = "unit_id",
    time_name = "time_id",
    outcome_name = "outcome",
    endogenous_names = c("endog_var1", "endog_var2"),
    covariate_names = NULL
)

output = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames, AR1_options = list(AR1_method = "GMM", AR1_IV_lags = 2, AR1_IV_outcome = TRUE))


panel_model = "AR1"
AR1_options = list(AR1_method = "GMM", AR1_IV_lags = 2, AR1_IV_outcome = TRUE)


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


beta_guess = c(0.5, -0.2)
delta_guess = c(1)
AR1_persistence_guess = 0.5

GMM_moments = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome)

GMM_means = GMM_moments$moment_mean
GMM_var = diag(GMM_moments$moment_var)
t(GMM_means) %*% solve(GMM_var) %*% GMM_means

V1<-rpois(5,10)
M1<-matrix(diag(V1),ncol=5)
M1

moment_var_diag = matrix(diag(moment_var), ncol = ncol(moment_var), nrow = nrow(moment_var))
