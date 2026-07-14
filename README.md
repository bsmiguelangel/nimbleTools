# pNimble

`pNimble` is an `R` package to run NIMBLE models in parallel.

This package is based on the routines available at https://github.com/MigueBeneito/pNimble. Earlier versions of these routines were numbered `0.1`, `0.2` and `0.3`, and the version history below starts at `0.4.0` with their development as an `R` package.

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
* Support for user-specified monitored variables.
* Support for replacing selected default NIMBLE samplers.
* Optional posterior summaries using `MCMCvis`, with safe error handling so that posterior samples are still returned if summaries cannot be calculated.
* Optional calculation of WAIC from posterior samples, with safe error handling so that posterior samples are still returned if WAIC cannot be calculated.
* `MCMCproblems()`, a function to identify parameters with problematic MCMC behaviour and optionally produce traceplots only for those parameters.
* Notification system using `ntfy`.
* `dcar_leroux()`, a Leroux CAR density function for use in NIMBLE models.
* `rcar_leroux()`, the corresponding random generation function required by NIMBLE. It currently generates independent normal values for compatibility with NIMBLE, not exact simulations from the Leroux CAR distribution.
* `lerouxObjects()`, a helper function to construct the objects required by `dcar_leroux()` from a binary neighbourhood matrix.
* Support for a zero-mean constraint in the Leroux CAR distribution through the `zero_mean` argument.
* Adaptation of the Leroux CAR distribution to support HMC methods, including the use of `ADbreak()` to avoid unnecessary derivatives.

## History

### Version `0.4.0`

* Confirmed compatibility with the CRAN version of `nimble`. The development version of `nimble` is no longer required for installing or using `pNimble`.
* Renamed the `sd.theta` argument in the Leroux CAR distribution to `sd`.
* Made the zero-mean constraint in the Leroux CAR distribution adaptive to the number of small areas.
* Added safe error handling for WAIC and posterior summaries, so that posterior samples are still returned when these calculations fail.
* Added `lerouxObjects()` to construct the objects required by `dcar_leroux()` from a binary neighbourhood matrix.
* Added `MCMCproblems()` to identify parameters with problematic MCMC behaviour and optionally produce traceplots only for those parameters.
* Added a simple random generation function for `rcar_leroux()` to improve compatibility with NIMBLE user-defined distributions.

## To do

* Consider a new package name, such as `nimbleTools`, as the package now includes tools beyond parallelisation.
* Extend `lerouxObjects()` so that Leroux CAR objects can also be constructed from the usual `WinBUGS` neighbourhood objects `adj` and `num`, in addition to a binary neighbourhood matrix.
* Implement exact random generation from the Leroux CAR distribution, following GMRF simulation methods such as those described by Rue and Held.
* Release version `1.0.0` when the package is ready for submission to CRAN.
