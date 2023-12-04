# Dearing and Blevins (2024): Monte Carlo Experiments

Dearing and Blevins (2024).
Efficient and Convergent Sequential Pseudo-Likelihood Estimation of Dynamic Discrete Games.
_Review of Economic Studies_.

## Overview

This subdirectory contains the replication code for the Monte Carlo
experiments in Appendix B.2 of the paper.  The simulations are based
on the model of Pesendorfer and Schmidt-Dengler (2010), a static game
with an unstable equilibrium leading to non-convergence of NPL.

## Data Availability

The Monte Carlo experiments in this paper do not involve analysis of external
data (i.e., the only data used are generated via simulation in the code).

## Software Requirements

The results of this section were produced using Matlab.  The code has
been verified to run in Matlab 2022b.

## Hardware Requirements

The Monte Carlo experiments for this section can be completed on a
modern laptop in approximately five minutes.

## Description of Code

Below is a complete list of program and function files used in the Monte Carlo
simulations.

-   `main.m` - Main Monte Carlo simulation and estimation program.
-   `LL_MLE.m` - Log likelihood function for MLE.
-   `LL_EPL.m` - Log likelihood function for EPL.
-   `LL_NPL.m` - Log likelihood function for NPL.
-   `myCDF.m` - Error CDF specified in Pesendorfer and Schmidt-Dengler (2010).

To replicate the Monte Carlo experiments, run `main` in Matlab.

## List of tables, figures, and programs

-   Table 17 in the paper is replicated by `main.m`.

### Results

The results in the paper were produced by running `main.m` described
above.  A log file produced by Matlab 2018b is included in the
replication package: `mc-psd2010.log`.

## References

-   Pesendorfer, M. and P. Schmidt-Dengler (2010).
    [Sequential estimation of dynamic discrete games: A comment](https://doi.org/10.3982/ECTA7633).
    _Econometrica_ 78, 833-842.
