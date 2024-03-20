
#' Internal function to estimate AR(1) model using GMM method. There are two cases depending on whether or not the AR(1) persistence parameter is known.
#' @noRd
PR.est.AR1.GMM <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_outcome = TRUE, AR1_IV_lags = 1, AR1_persistence = NULL) {
    if (is.null(AR1_persistence)) {
        initial_guess = PR.est.AR1.GMM.unknown_persistence.init_search(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome)
        print(sprintf("Best initial guess of beta parameters: %s. Best initial guess of AR1 persistence: %s.", paste0(round(initial_guess$beta_params, 3), collapse = ", "), round(initial_guess$AR1_persistence, 3))) 
        return(PR.est.AR1.GMM.unknown_persistence(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, beta_guess = initial_guess$beta_params, AR1_persistence_guess = initial_guess$AR1_persistence, delta_guess = c(0, initial_guess$delta_params)))
    } else {
        return(PR.est.AR1.GMM.known_persistence(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, AR1_persistence = AR1_persistence))
    }
}



#' Internal function to estimate AR(1) model. It uses GMM to fit the IV moments. For a guess of beta and the AR1_persistence, it computes the unobserved shock using quasi-differences. u_{i,t} = (y_{i,t} - rho*y_{i,t-1}) - (X_{i,t} - rho*X_{i,t-1})'beta. Then it fits mean(u_{i,t} * X_{i,t-1}) = 0 using BFGS.
#' @noRd
PR.est.AR1.GMM.unknown_persistence <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, AR1_persistence = NULL, beta_guess, AR1_persistence_guess, delta_guess) {
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
        GMM_obj = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = weighting_matrix, get_variance = FALSE)
        # compute the objective value
        return(GMM_obj)
    }
    # get an initial weighting matrix
    var_matrix = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = NULL, get_variance = TRUE)
    weighting_matrix = solve(var_matrix)
    # try up to 10 times to get convergence, with random starts
    startvals = c(beta_guess, AR1_persistence_guess, delta_guess)
    iter = 0
    convergence = 1
    lagged_convergence = 1
    best_val = GMM_objective(startvals, AR1_persistence = NULL, weighting_matrix = weighting_matrix)
    best_est = startvals
    while ((convergence != 0 || lagged_convergence != 0) && iter < 10) {
        # randomize the start values
        if (iter > 0) {
            startvals = startvals * runif(length(startvals), 0.75, 1.25)
        }
        # solve using BFGS, switch to Nelder-Mead if BFGS fails to converge
        if (iter < 5) {
            sol = optim(startvals, GMM_objective, AR1_persistence = NULL, weighting_matrix = weighting_matrix, method = "BFGS")
        } else {
            sol = optim(startvals, GMM_objective, AR1_persistence = NULL, weighting_matrix = weighting_matrix, method = "Nelder-Mead")
        }
        # unpack the parameters
        x = sol$par
        param_index = 1
        beta_guess = x[1:n_endogenous]
        param_index = param_index + n_endogenous
        AR1_persistence_guess = x[(n_endogenous + 1)]
        param_index = param_index + 1
        delta_guess = x[param_index:(param_index + n_exogenous)]

        # update the best solution
        if (sol$value < best_val) {
            best_val = sol$value
            best_est = sol$par
            # get an updated weighting matrix
            var_matrix = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = NULL, get_variance = TRUE)
            weighting_matrix = solve(var_matrix)
        }
        if (iter > 0) {
            lagged_convergence = convergence
        }
        convergence = sol$convergence
        iter = iter + 1
    }
    # check if the solution is converged
    if (convergence != 0) {
        stop("GMM did not converge")
    }
    # unpack the solution
    x = best_est
    param_index = 1
    beta_guess = x[1:n_endogenous]
    param_index = param_index + n_endogenous
    AR1_persistence_guess = x[(n_endogenous + 1)]
    param_index = param_index + 1
    delta_guess = x[param_index:(param_index + n_exogenous)]
    # return the best estimates
    best_estimates = data.table(Variable = c(endogenous_names, "AR1_persistence"), Estimate = round(c(beta_guess, AR1_persistence_guess), 6))
    # update the var matrix
    var_matrix = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = NULL, get_variance = TRUE)
    # prepare standard errors
    SEs = PR.est.AR1.GMM.SE(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, var_matrix = var_matrix)
    best_estimates = merge(best_estimates, SEs, by = "Variable")
    best_estimates[, tvalue := Estimate / SE]
    best_estimates[, pvalue := 2 * (1 - pt(abs(tvalue), N - 1))]
    best_estimates[, N := NULL]
    return(best_estimates)
}


#' Internal function to estimate AR(1) model. It uses GMM to fit the IV moments. For a guess of beta and the AR1_persistence, it computes the unobserved shock using quasi-differences. u_{i,t} = (y_{i,t} - rho*y_{i,t-1}) - (X_{i,t} - rho*X_{i,t-1})'beta. Then it fits mean(u_{i,t} * X_{i,t-1}) = 0 using BFGS.
#' @noRd
PR.est.AR1.GMM.known_persistence <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, AR1_persistence) {
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
        GMM_obj = PR.est.AR1.GMM_moments(panel_data = panel_data, outcome_name = outcome_name, endogenous_names = endogenous_names, covariate_names = covariate_names, fixedeffect_names = fixedeffect_names, beta_guess = beta_guess, delta_guess = delta_guess, AR1_persistence_guess = AR1_persistence_guess, AR1_IV_lags = AR1_IV_lags, AR1_IV_outcome = AR1_IV_outcome, weighting_matrix = weighting_matrix, get_variance = FALSE)
        # compute the objective value
        return(GMM_obj)
    }
    # solve the GMM objective
    startvals = NULL
    best_est = NULL
    sol = NULL
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
PR.est.AR1.GMM_moments <- function(panel_data, outcome_name, endogenous_names, covariate_names, fixedeffect_names, beta_guess, delta_guess, AR1_persistence_guess, AR1_IV_lags = 2, AR1_IV_outcome = TRUE, weighting_matrix = NULL, get_variance = FALSE) {
    #########################
    # quasi-differences
    #########################
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
    #########################
    # quasi-diff y minus quasi-diff X'beta minus quasi-diff intercept
    panel_data$u_it = panel_data[, get(paste0(outcome_name, "_quasidiff"))] - as.matrix(panel_data[, .SD, .SDcols = (paste0(endogenous_names, "_quasidiff"))]) %*% beta_guess  - delta_guess[1] * (1 - AR1_persistence_guess)
    # subtract the quasi-diff exogenous variables
    if (length(covariate_names) > 0) {
        panel_data$u_it = panel_data$u_it - as.matrix(panel_data[, .SD, .SDcols = (paste0(covariate_names, "_quasidiff"))]) %*% delta_guess[2:(length(delta_guess))]
    }
    if (length(fixedeffect_names) > 0) {
        panel_data$u_it = resid(feols(as.formula(sprintf("u_it ~ 1 | %s", fixedeffect_names)), data = panel_data))
    }
    #########################
    # compute the IV moments
    #########################
    moment_i = panel_data[, u_it]  # intercept
    for (ii in seq_len(length(endogenous_names))) { # endogenous variables, lag 1
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag"))])
    }
    if (AR1_IV_lags > 1) { # endogenous variables, lag 2
        for (ii in seq_len(length(endogenous_names))) {
            moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(endogenous_names[ii], "_lag2"))])
        }
    }
    if (AR1_IV_outcome) { # outcome, lag 1
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(outcome_name, "_lag"))])
    }
    if (AR1_IV_lags > 1 && AR1_IV_outcome) { # outcome, lag 2
        moment_i = cbind(moment_i, panel_data[, u_it * get(paste0(outcome_name, "_lag2"))])
    }
    if (length(covariate_names) > 0) { # exogenous variables, contemporaneous
        for (ii in seq_len(length(covariate_names))) {
            moment_i = cbind(moment_i, panel_data[, u_it * get(covariate_names[ii])])
        }
    }
    #########################
    # return the IV moments
    #########################
    moment_i = as.matrix(moment_i)
    # mean
    moment_mean = as.numeric(colMeans(moment_i))
    # variance
    if (get_variance) {
        moment_var = 0
        for (ii in seq_len(nrow(moment_i))) {
            moment_var = moment_var + (moment_i[ii, ] - moment_mean) %*% t(moment_i[ii, ] - moment_mean)
        }
        moment_var = moment_var / nrow(moment_i)
        return(moment_var)
    }
    # objective
    moment_obj = t(moment_mean) %*% weighting_matrix %*% moment_mean
    return(moment_obj)
}

PR.est.AR1.GMM.SE <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, beta_guess, delta_guess, AR1_persistence_guess, var_matrix) {
    ###########################
    # get the instrument matrix
    ###########################
    instrument_names = c(paste0(endogenous_names, "_lag"), covariate_names)
    if (AR1_IV_lags > 1) {
        instrument_names = c(instrument_names, paste0(endogenous_names, "_lag2"))
    }
    if (AR1_IV_outcome) { # outcome, lag 1
        instrument_names = c(instrument_names, paste0(outcome_name, "_lag"))
    }
    if (AR1_IV_lags > 1 && AR1_IV_outcome) { # outcome, lag 2
        instrument_names = c(instrument_names, paste0(outcome_name, "_lag2"))
    }
    if (length(covariate_names) > 0) { # exogenous variables, contemporaneous
        instrument_names = c(instrument_names, covariate_names)
    }
    instrument_matrix = panel_data[, .SD, .SDcols = instrument_names]
    instrument_matrix$ones = 1
    instrument_matrix = as.matrix(instrument_matrix)
    ##############################################
    # quasi-differences are needed for the moments
    ##############################################
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
    # rho-related difference: -y_{i,t-1} + x_{i,t-1}'beta + delta'x_{i,t}
    rho_diff = panel_data[, -get(paste0(outcome_name, "_lag"))] + as.matrix(panel_data[, mget(paste0(endogenous_names, "_lag"))]) %*% beta_guess + delta_guess[1]
    if (length(covariate_names) > 0) {
        rho_diff = rho_diff + as.matrix(panel_data[, mget(covariate_names)]) %*% delta_guess[2:length(delta_guess)]
    }
    #########################
    # compute the IV moments
    #########################
    IV_moments = NULL
    # endogenous quasi-differences
    for (ii in seq_len(length(endogenous_names))) {
        this_partial_eta = -panel_data[, get(paste0(endogenous_names[ii], "_quasidiff"))]
        moments_ii = c()
        for (zz in seq_len(ncol(instrument_matrix))) {
            moments_ii = c(moments_ii, mean(this_partial_eta * instrument_matrix[, zz]))
        }
        IV_moments = cbind(IV_moments, moments_ii)
    }
    # rho-related difference
    this_partial_eta = rho_diff
    moments_ii = c()
    for (zz in seq_len(ncol(instrument_matrix))) {
        moments_ii = c(moments_ii, mean(this_partial_eta * instrument_matrix[, zz]))
    }
    IV_moments = cbind(IV_moments, moments_ii)
    # intercept quasi-difference
    moments_ii = colMeans(-(1 - AR1_persistence_guess) * instrument_matrix)
    IV_moments = cbind(IV_moments, moments_ii)
    # exogenous quasi-differences
    for (ii in seq_len(length(covariate_names))) {
        this_partial_eta = -panel_data[, get(paste0(covariate_names[ii], "_quasidiff"))]
        moments_ii = c()
        for (zz in seq_len(ncol(instrument_matrix))) {
            moments_ii = c(moments_ii, mean(this_partial_eta * instrument_matrix[, zz]))
        }
        IV_moments = cbind(IV_moments, moments_ii)
    }
    # variance of the estimates, G'(Omega^-1)G / sample size
    G = as.matrix(IV_moments)
    var_coefs = (solve(t(G) %*% solve(var_matrix) %*% as.matrix(G))) / nrow(instrument_matrix)
    se_coefs = sqrt(diag(var_coefs))
    # return the standard errors
    the_names = c(endogenous_names, "AR1_persistence", "intercept")
    if (length(covariate_names) > 0) {
        the_names = c(the_names, covariate_names)
    }
    se_coefs = data.table(Variable = the_names, SE = se_coefs, N = nrow(instrument_matrix))
    return(se_coefs)
}


#' Internal function to estimate AR(1) model. It uses GMM to fit the IV moments. For a guess of beta and the AR1_persistence, it computes the unobserved shock using quasi-differences. u_{i,t} = (y_{i,t} - rho*y_{i,t-1}) - (X_{i,t} - rho*X_{i,t-1})'beta. Then it fits mean(u_{i,t} * X_{i,t-1}) = 0 using BFGS.
#' @noRd
PR.est.AR1.GMM.unknown_persistence.init_search <- function(panel_data, outcome_name, endogenous_names, covariate_names = NULL, fixedeffect_names = NULL, AR1_IV_lags = 1, AR1_IV_outcome = TRUE, AR1_persistence = NULL, num_restarts = 10) {
    # number of endogenous variables
    n_endogenous = length(endogenous_names)
    n_exogenous = length(covariate_names)
    instrument_names = c(paste0(endogenous_names, "_lag"))
    if (AR1_IV_outcome) { # outcome, lag 1
        instrument_names = c(instrument_names, paste0(outcome_name, "_lag"))
    }
    if (AR1_IV_lags > 1) {
        instrument_names = c(instrument_names, paste0(endogenous_names, "_lag2"))
    }
    if (AR1_IV_outcome && AR1_IV_lags > 1) { # outcome, lag 2
        instrument_names = c(instrument_names, paste0(outcome_name, "_lag2"))
    }
    if (length(covariate_names) > 0) { # exogenous variables
        instrument_names = c(instrument_names, covariate_names)
        instrument_names = c(instrument_names, paste0(covariate_names, "_lag"))
        if (AR1_IV_lags > 1) { # exogenous variables, lag 2
            instrument_names = c(instrument_names, paste0(covariate_names, "_lag2"))
        }
    }
    # verify that number of instruments is sufficient
    if (length(instrument_names) <= (n_endogenous + n_exogenous)) {
        stop("Number of instruments must be greater than the number of endogenous variables and covariates")
    }
    # calculate the moment conditions
    order_of_covariances = c(outcome_name, paste0(outcome_name, "_lag"), endogenous_names, paste0(endogenous_names, "_lag"))
    if (length(covariate_names) > 0) {
        order_of_covariances = c(order_of_covariances, covariate_names)
    }
    # we only need the number of instruments = number of parameters
    #instrument_names = instrument_names[1:(n_endogenous + 1 + n_exogenous)]
    # calculate the covariance moments
    covariance_moments = NULL
    for (zz in instrument_names){
        row_of_covariance_moments = c()
        for (var in order_of_covariances){
            row_of_covariance_moments = c(row_of_covariance_moments, panel_data[, cov(get(zz), get(var))])
        }
        covariance_moments = rbind(covariance_moments, row_of_covariance_moments)
    }
    rownames(covariance_moments) = instrument_names
    colnames(covariance_moments) = order_of_covariances
    # calculate the fit of the moments
    evaluate_fit <- function(x) {
        # beta and rho
        param_index = 1
        beta_guess = x[1:n_endogenous]
        param_index = param_index + n_endogenous
        AR1_persistence_guess = x[(n_endogenous + 1)]
        pvec = c(1, -AR1_persistence_guess, -beta_guess, beta_guess * AR1_persistence_guess)
        # delta
        if (n_exogenous > 0) {
            param_index = param_index + 1
            delta_guess = x[param_index:(param_index + n_exogenous - 1)]
            pvec = c(pvec, -delta_guess, delta_guess * AR1_persistence_guess)
        }
        # evaluation the moment fit
        Amat = t(covariance_moments) %*% (covariance_moments)
        pvec = matrix(pvec)
        g_eval = t(pvec) %*% Amat %*% pvec
        return(as.vector(g_eval))
    }
    # calculate the gradient
    evaluate_gradient <- function(x) {
      # beta and rho
      param_index = 1
      beta_guess = x[1:n_endogenous]
      param_index = param_index + n_endogenous
      AR1_persistence_guess = x[(n_endogenous + 1)]
      pvec = c(1, -AR1_persistence_guess, -beta_guess, beta_guess * AR1_persistence_guess)
      # delta
      if (n_exogenous > 0) {
        param_index = param_index + 1
        delta_guess = x[param_index:(param_index + n_exogenous - 1)]
        pvec = c(pvec, -delta_guess, delta_guess * AR1_persistence_guess)
      }
      # evaluation the moment fit
      Amat = t(covariance_moments) %*% (covariance_moments)
      pvec = matrix(pvec)
      # construct derivative matrix
      pvec_deriv = NULL
      for (ii in 1:length(endogenous_names)){
        beta_term = rep(0, length(endogenous_names))
        beta_term[ii] = 1
        this_pvec = c(0, 0, -beta_term, beta_term * AR1_persistence_guess)
        pvec_deriv = cbind(pvec_deriv, this_pvec)
      }
      beta_term = rep(0, length(endogenous_names))
      rho_deriv = c(0, -1, rep(0, length(endogenous_names)), beta_guess)
      pvec_deriv = cbind(pvec_deriv, rho_deriv)
      g_eval = t(pvec) %*% (2*Amat) %*% pvec_deriv
      return(as.vector(g_eval))
    }
    # minimize the quadratic form
    results = NULL
    fits = c()
    for (ii in 1:num_restarts) {
        set.seed(ii)
        startvals = rep(NA, n_endogenous + 1 + n_exogenous)
        soln = pso::psoptim(par = startvals * NA, fn = evaluate_fit, gr = evaluate_gradient)
        results = rbind(results, soln$par)
        fits = c(fits, sum(abs(soln$value)))
    }
    results_round = unique(round(cbind(results, fits), 3))
    results_round = results_round[results_round[, "fits"] == min(results_round[, "fits"]), ]
    if (!("numeric" %in% class(results_round))) {
      print("Warning: best solution is not unique.")
      print(results_round)
      results_round = results[results[, "fits"] == min(results[, "fits"]), ]
      results_round = results_round[1,]
    }
    if ("numeric" %in% class(results_round)) {
        best_params = results_round[1:(n_endogenous + 1 + n_exogenous)]
        beta_guess = best_params[1:n_endogenous]
        AR1_persistence_guess = best_params[(n_endogenous + 1)]
        delta_guess = NULL
        if (n_exogenous > 0) {
            delta_guess = best_params[(n_endogenous + 2):(length(best_params))]
        }
        return(list(beta_params = beta_guess, AR1_persistence = AR1_persistence_guess, delta_params = delta_guess))
    }
}

