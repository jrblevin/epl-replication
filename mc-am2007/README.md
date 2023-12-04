# Dearing and Blevins (2024): Monte Carlo Simulations

Dearing and Blevins (2024).
Efficient and Convergent Sequential Pseudo-Likelihood Estimation of Dynamic Discrete Games.
_Review of Economic Studies_.

## Overview

This subdirectory contains the replication code for the Monte Carlo experiments in
Section 4 of the paper.  The simulations are based on the dynamic game of entry
and exit in Example 1 of the paper, parameterized to match a model with five
heterogeneous firms from Aguirregabiria and Mira (2007).

## Data Availability

The Monte Carlo experiments in this paper do not involve analysis of external
data (i.e., the only data used are generated via simulation in the code).

## Software Requirements

Matlab is required to replicate the results in this section.  The code
has been verified to be compatible with Matlab 2022b.

## Hardware Requirements

The complete series of Monte Carlo experiments for this section can be
completed on a modern laptop in approximately three to five hours.

## Description of Code

Below is a complete list of program and function files used in the Monte Carlo
simulations.

Monte Carlo experiments:

- `mc_main.m`: Control program that runs all Monte Carlo simulations (across all
  experiments and sample sizes) and saves the results.
- `mc_epl.m`: The primary function which carries out each NPL and EPL Monte
  Carlo experiment for a fixed sample size and experiment number.
- `Gfunc.m`: Evaluates the equilibrium conditions of the dynamic entry/exit
  oligopoly model.
- `epldygam.m`: Estimates the structural parameters using EPL.
- `npldygam.m`: Estimates the structural parameters using NPL.
- `freqprob.m`: Frequency probabilty estimation.
- `milogit.m`: Maximum Likelihood Estimation of Logit Model.
- `clogit.m`: Maximum Likelihood estimation of McFadden's Conditional Logit
  model.
- `loglogit.m`: Log likelihood and gradient functions for multinomial logit.
- `simdygam.m`: Simulates data from the dynamic game.
- `mpestat.m`: Sample summary statistics.

Post-estimation programs to generate figures:

- `iter_histograms.m`: Produces histograms of the number of iterations taken by
  the EPL and NPL estimators across Monte Carlo replications.
- `parameter_histograms.m`: Produces histograms of the paramter estimates across
  Monte Carlo replications.
- `time_histograms.m`: Produces histograms of the runtimes of the EPL and NPL
  estimators across Monte Carlo replications.

## List of tables, figures, and programs

- Tables 1-6 in the paper are produced by `mc_main.m`.
- Figure 1 in the paper is produced by `parameter_histograms.m`.
- Figure 2 in the paper is produced by `iter_histograms.m` and
  `time_histograms.m`.

## Results

The `results` subdirectory contains the output of all Monte Carlo experiments
reported in the paper. Equilibria for each experiment 1-3 are computed by the
`mc_epl.m` program and stored in the files `mc_epl_eq1.mat`, `mc_epl_eq2.mat`,
and `mc_epl_eq3.mat` respectively, if they do not already exist.  If they exist,
the equilibrium choice probabilities are loaded reused for computational
efficiency.

Each of the log files named `mc_epl_exper_J_N_obs.log`, contain the results of
1000 replications single Monte Carlo replication for experiment `J` (1, 2, or 3)
and sample size `N` (1600 or 6400). These log files are large, so they have
been gzip compressed to save space. The raw Matlab data files containing the
parameter estimates, number of iterations, and runtimes from each replication
are saved in .mat files named `mc_epl_exper_J_N_obs.mat`. The parameter, time,
and iteration histograms are saved in PDF files with similar names.

## References

-   Aguirregabiria, V., and P. Mira (2007).
    [Sequential Estimation of Dynamic Discrete Games](https://doi.org/10.1111/j.1468-0262.2007.00731.x).
    _Econometrica_ 75, 1-53.
