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
