rm(list = ls())
setwd("~/github/PanelReg/")
library(PanelReg)
library(parallel)
library(data.table)

varnames = list(
  id_name = "unit_id",
  time_name = "time_id",
  outcome_name = "outcome",
  endogenous_names = c("endog_var1", "endog_var2"),
  covariate_names = NULL
)

eval_AR1 <- function(ii) {
    print(ii)
    # read values
    nn = results_holder$N[ii]
    seedval = results_holder$seedval[ii]
    method = results_holder$Method[ii]

    # simulate example data
    simulated_data = PanelRegSim(sample_size = nn, min_year = 2005, max_year = 2007, true_beta = c(0.5, -0.2), panel_model = "AR1", noise_sd = 0, seed = seedval)

    # estimate the model under the AR1 assumption, IV regression
    if (method == "PanelIV") {
        AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = FALSE)
        this_est = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    }
    if (method == "PanelIVy") {
        AR1_options = list(AR1_method = "PanelIV", AR1_IV_outcome = TRUE)
        this_est = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    }
    if (method == "GMM_lag1y") {
        AR1_options = list(AR1_method = "GMM", AR1_IV_outcome = TRUE, AR1_IV_lags = 1)
        this_est = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    }
    if (method == "GMM_lag2") {
        AR1_options = list(AR1_method = "GMM", AR1_IV_outcome = FALSE, AR1_IV_lags = 2)
        this_est = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    }
    if (method == "GMM_lag2y") {
        AR1_options = list(AR1_method = "GMM", AR1_IV_outcome = TRUE, AR1_IV_lags = 2)
        this_est = PanelReg(panel_data = copy(simulated_data), panel_model = "AR1", varnames = varnames, AR1_options = AR1_options)
    }
    # return the estimate
    ests = as.numeric(this_est[Variable  %in% c("endog_var1", "endog_var2")]$Estimate)
    res = data.table(Method = method, N = nn, Seed = seedval, Estimate1 = ests[1], Estimate2 = ests[2])
    return(res)
}

# set up the results holder
results_holder = data.table(expand.grid(seedval = 1:10, N = c(1e3, 2e3, 5e3, 1e4, 2e4), Method = c("PanelIV", "PanelIVy", "GMM_lag1y", "GMM_lag2", "GMM_lag2y")))

# run the simulation
aa = proc.time()[3]
#res = lapply(150:180, eval_AR1)
res = mclapply(seq_len(nrow(results_holder)), eval_AR1, mc.cores = 7)
bb = proc.time()[3] 
print(sprintf("elapsed time: %s minutes", round((bb - aa)/60, 2)))
res = rbindlist(res)
write.csv(res, "inst/AR1_simulation_results.csv", row.names = FALSE)



# plot the results
setwd("~/github/PanelReg/")
library(data.table)
library(ggplot2)
library(latex2exp)
res = setDT(read.csv("inst/AR1_simulation_results.csv"))

# prepare the results
options(scipen = 999)
res[, Size := N / 1e3]
res[, Size := factor(Size, levels = res[, sort(unique(Size))])]
res[, Approach := ""]
res[Method == "PanelIV", Approach := "Panel IV: X only"]
res[Method == "PanelIVy", Approach := "Panel IV: X and Y"]
res[Method == "GMM_lag1", Approach := "GMM, IV: Lag 1 of X"]
res[Method == "GMM_lag1y", Approach := "GMM, IV: Lag 1 of X & Y"]
res[Method == "GMM_lag2", Approach := "GMM, IV: Lag 1-2 of X"]
res[Method == "GMM_lag2y", Approach := "GMM, IV: Lag 1-2 of X & Y"]

gg = ggplot(res, aes(x = Size, y = Estimate1, color = Approach)) +
  geom_boxplot(width = .3) + geom_hline(yintercept = 0.5) +
  labs(title = TeX("AR1 Simulation Results: $beta_1$"), x = "Sample Size (Thousands)", y = "Estimate") +
  theme_bw() +
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = seq(-0.35, 0.9, by = 0.1), limits = c(-0.35, 0.9))
ggsave(filename = "inst/AR1_simulation_results1.png", gg, width = 7, height = 5)

gg = ggplot(res, aes(x = Size, y = Estimate2, color = Approach)) +
  geom_boxplot(width = .3) + geom_hline(yintercept = -0.2) +
  labs(title = TeX("AR1 Simulation Results: $beta_2$"), x = "Sample Size (Thousands)", y = "Estimate") +
  theme_bw() +
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = seq(-0.5, 0, by = 0.05), limits = c(-0.5, 0))
ggsave(filename = "inst/AR1_simulation_results2.png", gg, width = 7, height = 5)
