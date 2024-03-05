
#' Internal function to estimate AR(1) model. It uses GMM to fit the IV moments. For a guess of beta and the AR1_persistence, it computes the unobserved shock using quasi-differences. u_{i,t} = (y_{i,t} - rho*y_{i,t-1}) - (X_{i,t} - rho*X_{i,t-1})'beta. Then it fits mean(u_{i,t} * X_{i,t-1}) = 0 using BFGS.
#' @noRd
PR.est.AR1.GMM <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, AR1_persistence = NULL) {
    # not implemented yet
    if (!is.null(fixedeffect_names)) {
        stop("Stop: Fixed effects not implemented yet for AR(1) GMM estimator.")
    }
    # number of endogenous variables
    n_endogenous = length(endogenous_names)
    n_exogenous = length(covariate_names)
    # define the quasi-difference outcome and endogenous variables
    GMM_objective = function(x, AR1_persistence = NULL, weighting_matrix = NULL) {
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
        GMM_moments = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = weighting_matrix)
        # compute the objective value
        return(GMM_moments$moment_obj)
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
        # unpack the parameters
        beta_guess = startvals[1:n_endogenous]
        delta_guess = startvals[(n_endogenous + 1):(n_endogenous + 1 + n_exogenous)]
        # get an initial weighting matrix
        weighting_matrix = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = 0.0, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = NULL)$moment_var
        weighting_matrix = solve(diag(diag(weighting_matrix))) 
        # solution if AR1_persistence is zero
        sol = optim(startvals, GMM_objective, AR1_persistence = 0.0, weighting_matrix = weighting_matrix, method = "Nelder-Mead")
        startvals = sol$par
        startvals = c(startvals[1:n_endogenous], 0.0, startvals[(n_endogenous + 1):length(startvals)])
        # initial solution with free AR1_persistence parameter
        sol = optim(startvals, GMM_objective, AR1_persistence = NULL, weighting_matrix = weighting_matrix, method = "Nelder-Mead")
        startvals = sol$par
        best_est = sol$par
        best_val = sol$value 
        # unpack the parameters
        x = startvals
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
        # get an updated weighting matrix
        weighting_matrix = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = NULL)$moment_var
        #weighting_matrix = solve(diag(diag(weighting_matrix))) 
        weighting_matrix = solve(weighting_matrix)
        # try up to 10 times to get convergence, with random starts
        iter = 0
        convergence = 1
        while (convergence != 0 && iter < 10) {
            startvals = sol$par
            startvals = startvals * runif(length(startvals), 0.75, 1.25)
            if (iter < 5) {
                sol = optim(startvals, GMM_objective, AR1_persistence = NULL, weighting_matrix = weighting_matrix, method = "BFGS")
            } else {
                sol = optim(startvals, GMM_objective, AR1_persistence = NULL, weighting_matrix = weighting_matrix, method = "Nelder-Mead")
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
PR.est.AR1.GMM_moments <- function(panel_data, outcome_name, endogenous_names, covariate_names, beta_guess, delta_guess, AR1_persistence_guess, AR1_IV_lags = 2, AR1_IV_outcome = TRUE, weighting_matrix = NULL) {
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
    moment_i = panel_data[, u_it]  # intercept
    for (ii in seq_len(length(endogenous_names))) { # endogenous variables, lag 1
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag"))])
    }
    if (AR1_IV_lags >= 1) { # endogenous variables, lag 2
        for (ii in seq_len(length(endogenous_names))) {
            moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag2"))])
        }
    }
    if (AR1_IV_outcome) { # outcome, lag 1
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(outcome_name, "_lag"))])
    }
    if (AR1_IV_lags >= 1 && AR1_IV_outcome) { # outcome, lag 2
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(outcome_name, "_lag2"))])
    }
    if (length(covariate_names) > 0) { # exogenous variables, contemporaneous
        for (ii in seq_len(length(covariate_names))) {
            moment_i = cbind(moment_i, panel_data[, u_it * get(covariate_names[ii])])
        }
    }
    #########################
    # return the IV moments
    ########################
    moment_mean = as.numeric(colMeans(moment_i))
    moment_i = as.matrix(moment_i)
    moment_var = t(moment_i) %*% moment_i / nrow(panel_data)
    if (is.null(weighting_matrix)) {
        weighting_matrix = solve(moment_var)
    } 
    moment_obj = t(moment_mean) %*% weighting_matrix %*% moment_mean
    return(list(moment_mean = moment_mean, moment_var = moment_var, moment_obj = moment_obj))
}
