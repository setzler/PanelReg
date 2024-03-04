#' Internal function to estimate AR(1) model using Panel IV method. There are two cases depending on whether or not the AR(1) persistence parameter is known.
#' @noRd
PR.est.AR1.PanelIV <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_outcome = TRUE, AR1_persistence = NULL) {
    if (is.null(AR1_persistence)) {
        return(PR.est.AR1.PanelIV.unknown_persistence(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_outcome = AR1_IV_outcome))
    } else {
        return(PR.est.AR1.PanelIV.known_persistence(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_persistence = AR1_persistence))
    }
}

#' Internal function to estimate AR(1) model. The formula is y_{i,t} = rho*y_{i,t-1} + X_{i,t}'beta - X_{i,t-1}'(beta*rho) + u_{i,t}, where X_{i,t} depends on u_{i,t}. The instruments are X_{i,t-1}.
#' @noRd
PR.est.AR1.PanelIV.unknown_persistence <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_outcome = TRUE) {
    # build the formula
    the_formula = PR.formula.AR1.PanelIV.unknown_persistence(outcome_name, endogenous_names, covariate_names, fixedeffect_names, AR1_IV_outcome)
    # estimate the model
    est = feols(as.formula(the_formula), data = panel_data)
    # return the estimated second stage
    secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
    setnames(secondstage, c("Variable", "Estimate", "SE", "tvalue", "pvalue"))
    secondstage[, Variable := gsub("fit_", "", Variable)]
    secondstage[Variable == paste0(outcome_name, "_lag"), Variable := "AR1_persistence"]
    return(secondstage)
}

#' Internal function to estimate AR(1) model. It takes the AR(1) persistence parameter as given, constructs the quasi-differenced outcome and endogenous variables, then uses an IV regression to estimate the effect beta. The formula is (y_{i,t} - rho*y_{i,t-1}) = (X_{i,t} - rho*X_{i,t-1})'beta + u_{i,t}, where X_{i,t} depends on u_{i,t}. The instruments are X_{i,t-1}.
#' @noRd
PR.est.AR1.PanelIV.known_persistence <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_persistence, AR1_IV_outcome = TRUE) {
    # define quasi-differenced outcome and endogenous variables
    panel_data[, (paste0(outcome_name, "_quasidiff")) := get(outcome_name) - AR1_persistence  * get(paste0(outcome_name, "_lag"))]
    for (ii in seq_len(length(endogenous_names))) {
        panel_data[, paste0(endogenous_names[ii], "_quasidiff") := get(endogenous_names[ii]) - AR1_persistence * get(paste0(endogenous_names[ii], "_lag"))]
    }
    if (length(covariate_names) > 0) {
        for (ii in seq_len(length(covariate_names))) {
            panel_data[, paste0(covariate_names[ii], "_quasidiff") := get(covariate_names[ii]) - AR1_persistence * get(paste0(covariate_names[ii], "_lag"))]
        }
    }
    # build the formula
    the_formula = PR.formula.AR1.PanelIV.known_persistence(outcome_name, endogenous_names, covariate_names, fixedeffect_names, AR1_IV_outcome)
    # estimate the model
    est = feols(as.formula(the_formula), data = panel_data)
    # return the estimated second stage
    secondstage = cbind(data.table(Variable = row.names(est$coeftable)), as.data.table(est$coeftable))
    secondstage[, Variable := gsub("fit_", "", Variable)]
    secondstage[, Variable := gsub("_quasidiff", "", Variable)]
    return(secondstage)
}

#' Internal function to estimate AR(1) model. It uses GMM to fit the IV moments. For a guess of beta and the AR1_persistence, it computes the unobserved shock using quasi-differences. u_{i,t} = (y_{i,t} - rho*y_{i,t-1}) - (X_{i,t} - rho*X_{i,t-1})'beta. Then it fits mean(u_{i,t} * X_{i,t-1}) = 0 using BFGS.
#' @noRd
PR.est.AR1.GMM <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, AR1_persistence = NULL) {
    # number of endogenous variables
    n_endogenous = length(endogenous_names)
    n_exogenous = length(covariate_names)
    # define the quasi-difference outcome and endogenous variables
    GMM_objective = function(x, AR1_persistence = NULL) {
        # unpack the parameters
        param_index = 1
        beta_guess = x[1:n_endogenous]
        param_index = param_index + n_endogenous
        if (!is.null(AR1_persistence)) {
            AR1_persistence_guess = AR1_persistence
        } else {
            AR1_persistence_guess = x[(n_endogenous + 1)]
            param_index = param_index + 1
        }
        delta_guess = x[param_index:(param_index + n_exogenous)]
        # compute the IV moments
        GMM_moments = PR.est.AR1.GMM_moments(panel_data, outcome_name, endogenous_names, covariate_names, beta_guess, delta_guess, AR1_persistence_guess)
        # extract subset of GMM moments
        moments_to_keep = c("(Intercept)")
        if (AR1_IV_lags >= 1) {
            moments_to_keep = c(moments_to_keep, paste0(endogenous_names, "_lag"))
            if (AR1_IV_outcome) {
                moments_to_keep = c(moments_to_keep, paste0(outcome_name, "_lag"))
            }
        }
        if (AR1_IV_lags >= 2) {
            moments_to_keep = c(moments_to_keep, paste0(endogenous_names, "_lag2"))
            if (AR1_IV_outcome) {
                moments_to_keep = c(moments_to_keep, paste0(outcome_name, "_lag2"))
            }
        }
        # compute the objective value
        fit_vector = GMM_moments[Variable %in% moments_to_keep]$IV_moment
        objval = sum(fit_vector^2)
        return(objval)
    }
    # solve the GMM objective
    startvals = NULL
    best_est = NULL
    sol = NULL
    if (is.null(AR1_persistence)) {
        # set up zeros
        startvals = rep(0, n_endogenous) # endogenous variables
        startvals = c(startvals, 0.0) # intercept
        if (n_exogenous > 0) {
            startvals = c(startvals, rep(0, n_exogenous))
        }
        # solution if AR1_persistence is zero
        sol = optim(startvals, GMM_objective, AR1_persistence = 0.0, method = "Nelder-Mead")
        startvals = sol$par
        sol = optim(startvals, GMM_objective, AR1_persistence = 0.0, method = "BFGS")
        # set up start values with 0 as the persistence parameter
        startvals = sol$par
        startvals = c(startvals[1:n_endogenous], 0.0, startvals[(n_endogenous + 1):length(startvals)])
        # initial solution with free AR1_persistence parameter
        sol = optim(startvals, GMM_objective, AR1_persistence = NULL, method = "Nelder-Mead")
        startvals = sol$par
        best_est = sol$par
        best_val = sol$value
        # try up to 10 times to get convergence, with random starts
        iter = 0
        convergence = 1
        while (convergence != 0 && iter < 10) {
            startvals = sol$par
            startvals = startvals * runif(length(startvals), 0.75, 1.25)
            if (iter < 5) {
                sol = optim(startvals, GMM_objective, AR1_persistence = NULL, method = "BFGS")
            } else {
                sol = optim(startvals, GMM_objective, AR1_persistence = NULL, method = "Nelder-Mead")
            }
            convergence = sol$convergence
            if (sol$value < best_val) {
                best_val = sol$value
                best_est = sol$par
            }
            iter = iter + 1
        }
        # check if the solution is converged
        if (sol$value < best_val) {
            best_val = sol$value
            best_est = sol$par
        }
        if (convergence != 0) {
            stop("GMM did not converge")
        }
    }
    if (!is.null(AR1_persistence)) {
        # set up zeros
        startvals = rep(0, n_endogenous) # endogenous variables
        startvals = c(startvals, 0.0) # intercept
        if (n_exogenous > 0) {
            startvals = c(startvals, rep(0, n_exogenous))
        }
        # solution if AR1_persistence is zero
        sol = optim(startvals, GMM_objective, AR1_persistence = AR1_persistence, method = "Nelder-Mead")
        startvals = sol$par
        iter = 0
        while (sol$convergence != 0 && iter < 10) {
            startvals = sol$par
            startvals = startvals * runif(length(startvals), 0.75, 1.25)
            sol = optim(startvals, GMM_objective, AR1_persistence = AR1_persistence, method = "BFGS")
            iter = iter + 1
        }
        best_est = sol$par
        # check if the solution is converged
        if (sol$convergence != 0) {
            stop("GMM did not converge")
        }
    }

    # unpack the solution
    x = best_est
    param_index = 1
    beta_guess = x[1:n_endogenous]
    param_index = param_index + n_endogenous
    if (!is.null(AR1_persistence)) {
        AR1_persistence_guess = AR1_persistence
    } else {
        AR1_persistence_guess = x[(n_endogenous + 1)]
        param_index = param_index + 1
    }
    delta_guess = x[param_index:(param_index + n_exogenous)]
    # return the best estimates
    best_estimates = data.table(Variable = c(endogenous_names, "AR1_persistence"), Estimate = round(c(beta_guess, AR1_persistence_guess), 6))
    return(best_estimates)
}

#' Internal function to estimate AR(1) model. It computes the IV moments for a guess of beta and the AR1_persistence. The IV moments are mean(u_{i,t} * X_{i,t-1}) for each endogenous variable and its lag.
#' @noRd
PR.est.AR1.GMM_moments <- function(panel_data, outcome_name, endogenous_names, covariate_names, beta_guess, delta_guess, AR1_persistence_guess) {
    #########################
    # quasi-differences
    ########################
    # outcome
    panel_data[, (paste0(outcome_name, "_quasidiff")) := get(outcome_name) - AR1_persistence_guess  * get(paste0(outcome_name, "_lag"))]
    # endogenous variables
    for (ii in seq_len(length(endogenous_names))) {
        panel_data[, paste0(endogenous_names[ii], "_quasidiff") := get(endogenous_names[ii]) - AR1_persistence_guess * get(paste0(endogenous_names[ii], "_lag"))]
    }
    # exogenous variables
    if (length(covariate_names) > 0) {
        for (ii in seq_len(length(covariate_names))) {
            panel_data[, paste0(covariate_names[ii], "_quasidiff") := get(covariate_names[ii]) - AR1_persistence_guess * get(paste0(covariate_names[ii], "_lag"))]
        }
    }
    #########################
    # construct the shock
    ########################
    # quasi-diff y minus quasi-diff X'beta minus quasi-diff intercept
    panel_data$u_it = panel_data[, get(paste0(outcome_name, "_quasidiff"))] - as.matrix(panel_data[, .SD, .SDcols = (paste0(endogenous_names, "_quasidiff"))]) %*% beta_guess  - delta_guess[1] * (1 - AR1_persistence_guess)
    # subtract the quasi-diff exogenous variables
    if (length(covariate_names) > 0) {
        panel_data$u_it = panel_data$u_it - as.matrix(panel_data[, .SD, .SDcols = (paste0(covariate_names, "_quasidiff"))]) %*% delta_guess[2:(length(delta_guess))]
    }
    #########################
    # compute the IV moments
    ########################
    # intercept
    IV_moments = data.table(Variable = "(Intercept)", IV_moment = mean(panel_data[, u_it]))
    # endogenous variables, lag 1
    for (ii in seq_len(length(endogenous_names))) {
        IV_moments = rbind(IV_moments, data.table(Variable = paste0(endogenous_names[ii], "_lag"), IV_moment = mean(panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag"))])))
    }
    # endogenous variables, lag 2
    for (ii in seq_len(length(endogenous_names))) {
        IV_moments = rbind(IV_moments, data.table(Variable = paste0(endogenous_names[ii], "_lag2"), IV_moment = mean(panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag2"))])))
    }
    # outcome, lag 1
    IV_moments = rbind(IV_moments, data.table(Variable = paste0(outcome_name, "_lag"), IV_moment = mean(panel_data[, u_it * get(paste0(outcome_name, "_lag"))])))
    # outcome, lag 2
    IV_moments = rbind(IV_moments, data.table(Variable = paste0(outcome_name, "_lag2"), IV_moment = mean(panel_data[, u_it * get(paste0(outcome_name, "_lag2"))])))
    # exogenous variables, lag 1
    if (length(covariate_names) > 0) {
        for (ii in seq_len(length(covariate_names))) {
            IV_moments = rbind(IV_moments, data.table(Variable = paste0(covariate_names[ii], "_lag"), IV_moment = mean(panel_data[, u_it * get(paste0(covariate_names[ii], "_lag"))])))
        }
    }
    #########################
    # return the IV moments
    ########################
    return(IV_moments)
}
