


test_that("Estimation with covariates, iid", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_gamma = c(1, -2), panel_model = "iid", return_unobservables = FALSE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the FE assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_false(all(abs(est_iid$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    varnames$covariate_names = c("co_var1")
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_equal(est_iid$Estimate, c(0.5, -0.2), tolerance = 5e-2) # should pass

})


test_that("Estimation with covariates, FE", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_gamma = c(1, -2), panel_model = "FE", return_unobservables = FALSE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the FE assumption
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_false(all(abs(est_FE$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    varnames$covariate_names = c("co_var1")
    est_FE = PanelReg(panel_data = copy(simulated_data), panel_model = "FE", varnames = varnames)
    expect_equal(est_FE$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})

test_that("Estimation with covariates, MA1", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 20000, min_year = 2005, max_year = 2011, true_beta = c(0.5, -0.2), true_gamma = c(1, -0.5), panel_model = "MA1", return_unobservables = FALSE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_false(all(abs(est_MA1$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    varnames$covariate_names = c("co_var1")
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_equal(est_MA1$Estimate, c(0.5, -0.2), tolerance = 7e-2) # should pass

})


test_that("Estimation with covariates, AR1", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_gamma = c(1, -2), panel_model = "AR1", return_unobservables = FALSE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    varnames$covariate_names = c("co_var1")
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})

test_that("Estimation with covariates, AR1 known rho", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_gamma = c(1, -2), panel_model = "AR1", return_unobservables = FALSE)

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.5)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    varnames$covariate_names = c("co_var1")
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.5)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})




test_that("Estimation with group-year effects, iid", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_FE_var = 1, panel_model = "iid", return_unobservables = FALSE)
    simulated_data[, group_year := .GRP, list(group_id, time_id)]

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the FE assumption
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_false(all(abs(est_iid$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    varnames$fixedeffect_names = c("group_year")
    est_iid = PanelReg(panel_data = copy(simulated_data), panel_model = "iid", varnames = varnames)
    expect_equal(est_iid$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})

test_that("Estimation with group-year effects, MA1", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_FE_var = 1.0, panel_model = "MA1", return_unobservables = FALSE)
    simulated_data[, group_year := .GRP, list(group_id, time_id)]

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the MA1 assumption
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_false(all(abs(est_MA1$Estimate - c(0.5, -0.2)) < 1e-2)) # should fail

    varnames$fixedeffect_names = c("group_year")
    est_MA1 = PanelReg(panel_data = copy(simulated_data), panel_model = "MA1", varnames = varnames)
    expect_equal(est_MA1$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})


test_that("Estimation with group-year effects, AR1", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_FE_var = 1.0, panel_model = "AR1", return_unobservables = FALSE)
    simulated_data[, group_year := .GRP, list(group_id, time_id)]

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    varnames$fixedeffect_names = c("group_year")
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})


test_that("Estimation with group-year effects, AR1 known rho", {
    # simulate example data
    simulated_data = PanelRegSim(sample_size = 100000, min_year = 2005, max_year = 2008, true_beta = c(0.5, -0.2), true_FE_var = 1.0, panel_model = "AR1", return_unobservables = FALSE)
    simulated_data[, group_year := .GRP, list(group_id, time_id)]

    # define varnames
    varnames = list(
        id_name = "unit_id",
        time_name = "time_id",
        outcome_name = "outcome",
        endogenous_names = c("endog_var1", "endog_var2"),
        covariate_names = NULL
    )

    # estimate the model under the AR1 assumption, Panel IV method
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.5)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_false(all(abs(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate - c(0.5, -0.2)) < 5e-2)) # should fail

    varnames$fixedeffect_names = c("group_year")
    AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE, AR1_persistence = 0.5)
    est_AR1_PanelIV = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    expect_equal(est_AR1_PanelIV[Variable %in% c("endog_var1", "endog_var2")]$Estimate, c(0.5, -0.2), tolerance = 1e-2) # should pass

})
