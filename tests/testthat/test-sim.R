test_that("simulated data has the correct dimensions", {
    for (pmod in c("exogenous", "iid", "FE", "MA1", "AR1")) { 
        simulated_data = PanelRegSim(sample_size = 200, min_year = 2005, max_year = 2009, true_beta = c(0.5, 0.2), panel_model = pmod)
        expect_equal(nrow(simulated_data), 200 * (2009 - 2005 + 1))
        expect_equal(ncol(simulated_data), 5)
        expect_equal(names(simulated_data), c("unit_id", "time_id", "outcome", "endog_var1", "endog_var2"))
    }
})

test_that("PR.sim handles bad values correctly", {
  expect_error(PanelRegSim(sample_size = -1), "sample_size must be a positive integer")
  expect_error(PanelRegSim(min_year = 2010, max_year = 2000), "max_year must be greater than min_year")
  expect_error(PanelRegSim(noise_sd = -1), "noise_sd must be a non-negative number")
  expect_error(PanelRegSim(panel_model = "unknown"), "panel_model must be 'exogenous', 'iid', 'FE', 'MA1', or 'AR1'") 
  expect_error(PanelRegSim(endog_weight = -1), "endog_weight must be a number between 0 and 1") 
})
