# Dearing and Blevins (2024): Monte Carlo Experiments

Dearing and Blevins (2024).
Efficient and Convergent Sequential Pseudo-Likelihood Estimation of Dynamic Discrete Games.
_Review of Economic Studies_.

## Overview

This subdirectory contains the replication code for the Monte Carlo
experiments in Appendix B.1 of the paper.  The simulations are based
on the model of Pesendorfer and Schmidt-Dengler (2008), a dynamic game
with multiple equilibria.

## Data Availability

The Monte Carlo experiments in this paper do not involve analysis of external
data (i.e., the only data used are generated via simulation in the code).

## Software Requirements

Matlab 2018a was used to produce the results in this section.  The
code has also been tested with Matlab 2022b.

## Hardware Requirements

The complete series of Monte Carlo experiments for this section can be
completed on a modern laptop in approximately one to two hours.

## Description of Code

Below is a complete list of program and function files used in the
Monte Carlo simulations in Appendix B.1.

Main program:

- `main.m` - Main control program that executes all Monte Carlo
  experiments for multiple sample sizes (`nobs`), multiple equilibria
  (`eqm_dgp`), and multiple levels of noise in the estimates (`c`).
  Produces log files named according to the pattern
  `psd2008_estimate_eqmN.log` where `N` is the equilibrium index
  (`eqm_dgp` which can be 1, 2, or 3).  Results are also stored
  in `.mat` files for further processing if desired.

Main support functions:

- `psd2008_estimate.m` - Estimates the model using the generated data
  and stores results.  Requires setting `nsims` (= 1000) before running.
  Produces

- `generate_data.m` - Solves for equilibria and generates data for the
  Pesendorfer and Schmidt-Dengler (2008) model.  Requires setting
  `nsims` (number of simulations; `1000` in the paper) `nobs` (number
  of observations; `250` and `1000` in the paper), and `eqm_dgp`
  (equilibrium number; 1, 2, or 3) in the Matlab session before running
  directly.  Produces `psd2008_data.mat`.

Other support functions:

- `EqmDiff.m` - Equilibrium conditions for the model.
- `Gderiv.m` - Derivative of G function with respect to v.
- `LL_EPL.m` - Log likelihood for EPL estimation.
- `LL_NPL.m` - Log likelihood for NPL estimation.
- `normsurplus.m` - Social surplus function under the normal
  distribution.

## List of tables, figures, and programs

Tables 12-16 in the paper are produced by `psd2008_estimate_all.m`.
This script loops over all equilibria and sample sizes and estimates
the model under several scenarios:

- Table 12 corresponds to the case where `eqm_dgp = 1`,
  `c = 0.0`, `nstart = 1`, and `nobs = [ 250, 1000 ]`.

- Table 13 corresponds to the case where `eqm_dgp = 2`,
  `c = 0.0`, `nstart = 1`, and `nobs = [ 250, 1000 ]`.

- Table 14 corresponds to the case where `eqm_dgp = 3`,
  `c = 0.0`, `nstart = 1`, and `nobs = 1000`.

- Table 15 corresponds to the case where `eqm_dgp = [ 1, 2 ]`,
  `c = 0.5`, `nstart = 1`, and `nobs = 250`.

- Table 16 corresponds to the case where `eqm_dgp = [ 1, 2 ]`,
  `c = 1.0`, `nstart = 5`, and `nobs = [ 250, 1000 ]`.

### Results

The results in the paper were produced by running `main.m` described
above using Matlab 2018a.  The log files are included in the
replication package in the `results` subdirectory:

- `psd2008_estimate_eqm1.log`
- `psd2008_estimate_eqm2.log`
- `psd2008_estimate_eqm3.log`

## References

-   Pesendorfer, M. and P. Schmidt-Dengler (2008).
    [Asymptotic least squares estimators for dynamic games](https://doi.org/10.1111/j.1467-937X.2008.00496.x).
    _Review of Economic Studies_ 75, 901-928.
