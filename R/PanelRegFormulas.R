#' Internal function to create the formula for the "exogenous" panel model
#' @noRd
PR.formula.exogenous <- function(outcome_name, endogenous_names, covariate_names, fixedeffect_names) {
    the_formula = sprintf("%s ~ %s", outcome_name, paste(endogenous_names, collapse = " + "))
    if (length(covariate_names) > 0) {
        the_formula = sprintf("%s + %s", the_formula, paste(covariate_names, collapse = " + "))
    }
    if (length(fixedeffect_names) > 0) {
        the_formula = sprintf("%s | %s", the_formula, paste(fixedeffect_names, collapse = " + "))
    }
    return(the_formula)
}

#' Internal function to create the formula for the "iid" panel model
#' @noRd
PR.formula.iid <- function(outcome_name, endogenous_names, covariate_names, fixedeffect_names) {
    # covariate part of the formula
    covariate_formula = "1"
    if (length(covariate_names) > 0) {
        covariate_formula = paste0(covariate_formula, " + ", paste(covariate_names, collapse = " + "))
    }
    if (length(fixedeffect_names) > 0) {
        covariate_formula = sprintf("%s | %s", covariate_formula, paste(fixedeffect_names, collapse = " + "))
    }
    # main formula
    the_formula = sprintf("%s ~ %s | %s ~ %s", 
                            outcome_name, 
                            covariate_formula, 
                            paste(endogenous_names, collapse = " + "), 
                            paste(paste0(endogenous_names, "_lag"), collapse = " + ")
                        )
    return(the_formula)
}

#' Internal function to create the formula for the "FE" panel model
#' @noRd
PR.formula.FE <- function(outcome_name, endogenous_names, covariate_names, fixedeffect_names) {
    # covariate part of the formula
    the_formula = sprintf("%s ~ %s",
                              paste0(outcome_name, "_diff"),
                              paste0(paste0(endogenous_names, "_diff"), collapse = " + "))
    if (length(covariate_names) > 0) {
        the_formula = sprintf("%s + %s", the_formula, paste(paste0(covariate_names, "_diff"), collapse = " + "))
    }
    if (length(fixedeffect_names) > 0) {
        the_formula = sprintf("%s | %s", the_formula, paste(fixedeffect_names, collapse = " + "))
    }
    return(the_formula)
}

#' Internal function to create the formula for the "MA1" panel model
#' @noRd
PR.formula.MA1 <- function(outcome_name, endogenous_names, covariate_names, fixedeffect_names) {
    # covariate part of the formula
    covariate_formula = "1"
    if (length(covariate_names) > 0) {
        covariate_formula = paste0(covariate_formula, " + ", paste(covariate_names, collapse = " + "))
    }
    if (length(fixedeffect_names) > 0) {
        covariate_formula = sprintf("%s | %s", covariate_formula, paste(fixedeffect_names, collapse = " + "))
    }
    # main formula
    the_formula = sprintf("%s ~ %s | %s ~ %s", 
                            outcome_name, 
                            covariate_formula,
                            paste(endogenous_names, collapse = " + "),
                            paste(paste0(endogenous_names, "_lag2"), collapse = " + ")
                        )
    return(the_formula)
}