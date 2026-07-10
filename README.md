# pNimble

`pNimble` is an `R` package to run NIMBLE models in parallel.

This package is based on the routines available at https://github.com/MigueBeneito/pNimble

## Installation

`pNimble` can be installed from GitHub as follows:

```r
remotes::install_github("bsmiguelangel/pNimble")
```

## Current version

The current version of `pNimble` includes tools to run NIMBLE models in parallel and to use a Leroux CAR distribution in NIMBLE model code.

Main features currently included:

* `pNimble()`, a function to run independent NIMBLE MCMC chains in parallel.
* Optional use of HMC sampling through `nimbleHMC`.
* Support for user-defined monitors.
* Support for replacing default NIMBLE samplers.
* Optional posterior summaries using `MCMCvis`, with safe error handling so that posterior samples are still returned if summaries cannot be calculated.
* Optional calculation of WAIC from posterior samples, with safe error handling so that posterior samples are still returned if WAIC cannot be calculated.
* Notification system using `ntfy`.
* `dcar_leroux()`, a Leroux CAR density function for use in NIMBLE models.
* `rcar_leroux()`, the corresponding random generation function required by NIMBLE. Random generation is currently not implemented.
* Support for a zero-mean constraint in the Leroux CAR distribution through the `zero_mean` argument.
* Adaptation of the Leroux CAR distribution to support HMC methods, including the use of `ADbreak()` to avoid unnecessary derivatives.

## History

### Version 0.4.0

* Confirmed compatibility with the CRAN version of `nimble`. The development version of `nimble` is no longer required for installing or using `pNimble`.
* Renamed the `sd.theta` argument in the Leroux CAR distribution to `sd`.
* Added safe error handling for WAIC and posterior summaries, so that posterior samples are still returned when these calculations fail.

## To do

* Make the strength of the zero-mean constraint adaptive according to the number of small areas.
* Write a helper function to construct the Leroux CAR objects from a neighbourhood matrix `W`, including `from.to` and the eigenvalues of `D - W`.
* Ask the `nimble-users` Google Group about the `R CMD check` NOTE caused by assigning `dcar_leroux()` and `rcar_leroux()` to the global environment.
* Add a diagnostic traceplot function that only plots parameters with problematic MCMC behaviour, for example `Rhat > 1.02` or `n.eff < 400`.
* Modify `rcar_leroux()` so that it can generate random values from the Leroux CAR distribution.
* Submit the package to CRAN.
