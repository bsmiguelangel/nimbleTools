# nimbleTools

`nimbleTools` is an `R` package providing tools for running and diagnosing `NIMBLE` models.

This package is based on the routines available at https://github.com/MigueBeneito/pNimble. Earlier versions of these routines were numbered `0.1`, `0.2` and `0.3`, and the version history below starts at `0.4.0` with their development as an `R` package.

## Installation

`nimbleTools` can be installed from `GitHub` as follows:

```r
remotes::install_github("bsmiguelangel/nimbleTools")
```

## Current version

The current version of `nimbleTools` includes tools to run `NIMBLE` models in parallel, diagnose MCMC output, and use a Leroux CAR distribution in `NIMBLE` model code.

Main features currently included:

* `pNimble()`, a function to run independent `NIMBLE` MCMC chains in parallel.
* Optional use of HMC sampling through `nimbleHMC`.
* Support for user-specified monitored variables.
* Support for replacing selected default `NIMBLE` samplers.
* Optional posterior summaries using `MCMCvis`, with safe error handling so that posterior samples are still returned if summaries cannot be calculated.
* Optional calculation of WAIC from posterior samples, with safe error handling so that posterior samples are still returned if WAIC cannot be calculated.
* `MCMCproblems()`, a function to identify parameters with problematic MCMC behaviour and optionally produce traceplots only for those parameters.
* Notification system using `ntfy`.
* `dcar_leroux()`, a Leroux CAR density function for use in `NIMBLE` models.
* `rcar_leroux()`, the corresponding random generation function required by `NIMBLE`. It currently generates independent normal values for compatibility with `NIMBLE`, not exact simulations from the Leroux CAR distribution.
* `lerouxObjects()`, a helper function to construct the objects required by `dcar_leroux()` from different neighbourhood representations, including a binary neighbourhood matrix, `WinBUGS` objects `adj` and `num`, and graph objects.
* Support for a zero-mean constraint in the Leroux CAR distribution through the `zero_mean` argument.
* Adaptation of the Leroux CAR distribution to support HMC methods, including the use of `ADbreak()` to avoid unnecessary derivatives.

## History

### Version `0.4.0`

* Confirmed compatibility with the CRAN version of `nimble`. The development version of `nimble` is no longer required for installing or using `nimbleTools`.
* Renamed the `sd.theta` argument in the Leroux CAR distribution to `sd`.
* Made the zero-mean constraint in the Leroux CAR distribution adaptive to the number of small areas.
* Added safe error handling for WAIC and posterior summaries, so that posterior samples are still returned when these calculations fail.
* Added `lerouxObjects()` to construct the objects required by `dcar_leroux()` from different neighbourhood representations, including a binary neighbourhood matrix, `WinBUGS` objects `adj` and `num`, and graph objects.
* Added `MCMCproblems()` to identify parameters with problematic MCMC behaviour and optionally produce traceplots only for those parameters.
* Added a simple random generation function for `rcar_leroux()` to improve compatibility with `NIMBLE` user-defined distributions.

## To do

* Improve the speed of `MCMCproblems()` by using an existing posterior summary when available, instead of recalculating diagnostics unnecessarily.
* Extend `MCMCproblems()` so that selected parameters can be plotted directly through the `params` argument, regardless of their `Rhat` or `n.eff` values.
* Consider whether `MCMCproblems()` should be renamed to `MCMCmonitor()` in a future version.
* Check that the `monitors` argument in `pNimble()` only returns the variables requested by the user.
* Check the handling of indexed parameter names in `MCMCproblems()`, such as `beta[1]`, especially when using `ISB` and `exact` arguments.
* Move the registration of the Leroux CAR distribution to package loading, instead of registering and deregistering it every time `pNimble()` is called.
* Explore the use of `nimbleFunction` setup code to give the Leroux CAR distribution a more standard input format, closer to the ICAR distribution used by `NIMBLE`.
* Explore whether `nimble::as.carAdjacency()` can be used to construct the neighbourhood objects required by the Leroux CAR distribution.
* Add checks for the symmetry and validity of neighbourhood structures used by the Leroux CAR distribution.
* Explore the use of `getTimes` to measure computation times and identify possible bottlenecks in `pNimble()`.
* Consider including total computation time in the object returned by `pNimble()`.
* Check that, when `WAIC = TRUE`, `pNimble()` automatically includes the stochastic parent nodes of the data nodes in the monitored variables.
* Reduce unnecessary console messages during model loading, compilation and distribution registration when possible.
* Assess the sensitivity of the adaptation parameter used in the zero-mean constraint, currently fixed at `10`.
* Implement exact random generation from the Leroux CAR distribution, following GMRF simulation methods such as those described by Rue and Held.
* Release version `1.0.0` when the package is ready for submission to `CRAN`.
