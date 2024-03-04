test_that("Estimate 'exogenous' case", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), panel_model = "exogenous", return_unobservables = TRUE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the exogenous assumption
    est_exogenous = PanelReg(panel_data = copy(simulated_data), panel_model = "exogenous", varnames = varnames)
    expect_equal(est_exogenous$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

    # estimate the model under the iid assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_equal(est_iid$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_equal(est_FE$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_equal(est_MA1$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)  
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 2e-2) # should pass

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.0)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)  
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 2e-2) # should pass

})


test_that("Estimate 'iid' case", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 1e5, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), panel_model = "iid")

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )
    
    # estimate the model under the exogenous assumption
    est_exogenous = PanelReg(panel_data = copy(simulated_data), panel_model = "exogenous", varnames = varnames)
    expect_false(all(abs(est_exogenous$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the iid assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_equal(est_iid$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_false(all(abs(est_FE$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail
    
    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_equal(est_MA1$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 3e-2) # should pass

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.0)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 3e-2) # should pass

})


test_that("Estimate 'FE' case", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 1e5, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), panel_model = "FE")
    
    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the exogenous assumption
    est_exogenous = PanelReg(panel_data = copy(simulated_data), panel_model = "exogenous", varnames = varnames)
    expect_false(all(abs(est_exogenous$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the iid assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_false(all(abs(est_iid$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_equal(est_FE$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_false(all(abs(est_MA1$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail
    
    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.0)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

})

test_that("Estimate 'MA1' case", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 1e5, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), panel_model = "MA1")

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )
    
    # estimate the model under the exogenous assumption
    est_exogenous = PanelReg(panel_data = copy(simulated_data), panel_model = "exogenous", varnames = varnames)
    expect_false(all(abs(est_exogenous$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the iid assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_false(all(abs(est_iid$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_false(all(abs(est_FE$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail
    
    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_equal(est_MA1$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.0)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options) 
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

})




test_that("Estimate 'AR1' case", {
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
    
    # estimate the model under the exogenous assumption
    est_exogenous = PanelReg(panel_data = copy(simulated_data), panel_model = "exogenous", varnames = varnames)
    expect_false(all(abs(est_exogenous$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the iid assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_false(all(abs(est_iid$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_false(all(abs(est_FE$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail
    
    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_false(all(abs(est_MA1$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail
    
    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.5)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass
    
})
