Panel Regressions in R
================





## Overview

This package estimates panel regression models, which are models of the
form,
$y_{it} = x_{it}' \beta + w_{it}' \gamma + \epsilon_{it} + \nu_{it}$.
The parameter of interest is $\beta$, $x_{it}$ is the vector of
endogenous regressors, $w_{it}$ is the vector of exogenous covariates,
$\epsilon_{it}$ is the endogenous unobservable, and $\eta_{it}$ is the
exogenous unobservable (i.e. measurement error).

The purpose of this package is to provide a convenient estimator of
$\beta$ given the most common time-series process assumptions on
$\epsilon_{it}$, which are:

1.  Fixed Effects (“FE”): $\epsilon_{it} = \alpha_i$.
2.  Serially Independent (“iid”): $\epsilon_{it} = \eta_{it}$, where
    $\eta_{it}$ is exogenous.
3.  Moving Average of Order 1 (“MA1”):
    $\epsilon_{it} = \eta_{it} + \theta \eta_{it-1}$, where $\eta_{it}$
    is exogenous.
4.  First Order Autoregressive (“AR1”):
    $\epsilon_{it} = \rho \epsilon_{it-1} + \eta_{it}$, where
    $\eta_{it}$ is exogenous.

There are only 2 functions in this package:

1.  `PanelReg()`: This function estimates the panel regression.
    Documentation is available
    [here](https://setzler.github.io/PanelReg/reference/PanelReg.html)
    or by running `?PanelReg` in R.
2.  `PanelRegSim()`: This function simulates data from a panel
    regression model. Documentation is available
    [here](https://setzler.github.io/PanelReg/reference/PanelRegSim.html)
    or by running `?PanelRegSim` in R.

## Basic Usage

Before estimation, set up a variable list with the names of your
variables:

``` r
varnames = list(
  id_name = "id",
  time_name = "year",
  outcome_name = "Y",
  endogenous_names = c("X1", "X2"),
  covariate_names = c("W1")
)
```

To estimate the panel regression, use syntax such as the following:

``` r
reg = PanelReg(panel_data, panel_model = "MA1", varnames)
```

The argument `panel_model` must be one of `"exogenous"`, `"FE"`,
`"iid"`, `"MA1"`, or `"AR1"`.

For complete examples and details on the estimators, see the [Get
Started](https://setzler.github.io/PanelReg/articles/PanelReg.html)
page.

## Installation

To install the package from Github:

``` r
devtools::install_github("setzler/PanelReg")
```

To use the package after it is installed:

``` r
library(PanelReg)
```

Make sure these packages have been installed as well:

``` r
library(data.table)
library(fixest)
```
