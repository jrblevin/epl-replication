# Dearing and Blevins (2024) Replication Files

Dearing and Blevins (2024).
Efficient and Convergent Sequential Pseudo-Likelihood Estimation of
Dynamic Discrete Games.
_Review of Economic Studies_.

## Overview

This repository contains the replication code for all Monte Carlo
experiments and the Application of Dearing and Blevins (2024).  The
replication files for each section are contained in separate
subdirectories.  Here we give a brief overview of the contents of the
package.

-   The `mc-am2007` subdirectory contains the replication code for the
    Monte Carlo experiments in Section 4 of the paper.  The
    simulations are based on the dynamic game of entry and exit in
    Example 1 of the paper, parameterized to match a model with five
    heterogeneous firms from Aguirregabiria and Mira (2007).

-   The `application` subdirectory contains the replication code for
    the application to U.S.  wholesale club store competition in
    Section 5 of the paper.  The simulations are based on the dynamic
    game of entry and exit in Example 1 of the paper.

-   The `mc-psd2008` contains the replication code for the Monte Carlo
    experiments in Appendix B.1 of the paper.  The simulations are
    based on the model of Pesendorfer and Schmidt-Dengler (2008), a
    dynamic game with multiple equilibria.

-   The `mc-psd2010` subdirectory contains the replication code for
    the Monte Carlo experiments in Appendix B.2 of the paper.  The
    simulations are based on the model of Pesendorfer and
    Schmidt-Dengler (2010), a static game with an unstable equilibrium
    leading to non-convergence of NPL.

## Data Availability

The Monte Carlo experiments in this paper do not involve analysis of external
data (i.e., the only data used are generated via simulation in the code).

The application in the paper uses data from several sources:

-   [IPUMS NHGIS Table CL8: Total Population (1990, 2000, 2010, 2020)][nhgis]
-   [U.S. Department of Housing and Urban Development Zip Code to County Crosswalk Files][crosswalk]
-   [Data Axle Historical Business Database][infogroup]

We have included in this replication package the NHGIS population data
extract from Table CL8: Total Population (`application/nhgis0012_ts_geog2010_zcta.csv`)
as well as the HUD crosswalk file (`application/ZIP_COUNTY_122021.xlsx`).

The source files from the Data Axle Historical Business Database
cannot be made publicly available due to the terms of use.  While we
cannot redistribute the original dataset, we have provided our
anonymized analysis sample which users can use to replicate the tables
in the paper.

The full raw dataset can be obtained freely through certain university
libraries or directly from Data Axle directly for a fee.  Researchers
interested in access to the data may contact Data Axle at
<https://www.data-axle.com/contact-us/>.

The data preparation programs are also provided and can be used to
generate the analysis sample once certain data files have been
obtained from Data Axle.  These files are named following the pattern
`BUSINESS_HISTORICAL_YYYY.zip` for the years 1997 through 2020 and
`YYYY_Business_Academic_QCQ.txt.gz` for the year 2021.

   [nhgis]: https://www.nhgis.org/tabular-data-sources
   [crosswalk]: https://www.huduser.gov/portal/datasets/usps_crosswalk.html
   [infogroup]: https://www.data-axle.com/our-data/business-data/

### Details on each Data Source

| Data Name                    | Data Files                       | Location                     | Provided | Citation              |
|------------------------------|----------------------------------|------------------------------|----------|-----------------------|
| Zip Code Population          | `nhgis0012_ts_geog2010_zcta.csv` | `application/nhgis0012_csv/` | True     | Manson, et al. (2022) |
| Zip Code to County Crosswalk | `ZIP_COUNTY_122021.xlsx`         | `application/`               | True     | HUD (2021)            |
| Historical Business Database | Not available                    | `application/`               | False    | Data Axle (2022)      |

## Software Requirements

Matlab is required to replicate all Monte Carlo experiments and the
application.  All code has been verified to be compatible with Matlab
2022b.

Python and Stata were used to prepare the analysis sample for the
application.

The Python script `application/extract.py` imports the following
Python modules: `zipfile`, `gzip`, and `pandas`.  Although `zipfile`
and `gzip` are part of the standard library, Pandas is an external
package.  Depending on your Python installation, it can likely be
installed using the following command:

```python
pip install pandas
```

The Stata scripts used for the application require two additional
packages to be installed: `xttrans2` and `mat2txt`.  These can be
installed using the following commands:

```stata
net install xttrans2
ssc install mat2txt
```

The analysis sample used for the application in the paper was produced
using Python 3.11.6, Pandas 2.0.3, Stata MP 16.1, xttrans 1.2.0, and
mat2txt 1.1.2.

## Hardware Requirements

All of the programs in this replication package can be completed on a
modern laptop (e.g., a circa 2020 MacBook Pro) without issue.

-   The complete series of Monte Carlo experiments reported in Section
    4 (`mc-am2007` subdirectory) can be completed in approximately
    three to five hours.

-   The estimation and counterfactual simulations for the application
    (`application` subdirectory), including 250 bootstrap replications
    can be completed in approximately 30 minutes to one hour.

-   The complete series of Monte Carlo experiments reported in
    Appendix B.1 (`mc-psd2008` subdirectory) can be completed in
    approximately one to two hours.

-   The complete series of Monte Carlo experiments reported in
    Appendix B.2 (`mc-psd2010` subdirectory) can be completed on a
    modern laptop in approximately five minutes.

## Description of Code

### Monte Carlo Experiments (Section 4)

Below is a complete list of program and function files used in the
Monte Carlo simulations reported in Section 4.  These files are
contained in the `mc-am2007` subdirectory.

Monte Carlo experiments:

-   `mc_main.m`: Control program that runs all Monte Carlo simulations
    (across all experiments and sample sizes) and saves the results.
-   `mc_epl.m`: The primary function which carries out each NPL and
    EPL Monte Carlo experiment for a fixed sample size and experiment
    number.
-   `Gfunc.m`: Evaluates the equilibrium conditions of the dynamic
    entry/exit oligopoly model.
-   `epldygam.m`: Estimates the structural parameters using EPL.
-   `npldygam.m`: Estimates the structural parameters using NPL.
-   `freqprob.m`: Frequency probability estimation.
-   `milogit.m`: Maximum Likelihood Estimation of Logit Model.
-   `clogit.m`: Maximum Likelihood estimation of McFadden's
    Conditional Logit model.
-   `loglogit.m`: Log likelihood and gradient functions for
    multinomial logit.
-   `simdygam.m`: Simulates data from the dynamic game.
-   `mpestat.m`: Sample summary statistics.

Post-estimation programs to generate figures:

-   `iter_histograms.m`: Produces histograms of the number of
    iterations taken by the EPL and NPL estimators across Monte Carlo
    replications.
-   `parameter_histograms.m`: Produces histograms of the parameter
    estimates across Monte Carlo replications.
-   `time_histograms.m`: Produces histograms of the runtimes of the
    EPL and NPL estimators across Monte Carlo replications.

### Application (Section 5)

Below is a complete list of data, program, and function files used for
the application described in Section 5 of the paper.

#### Raw Data

-   `ZIP_COUNTY_122021.xlsx` - U.S. Department of Housing and Urban
    Development Zip Code to County Crosswalk File.

-   `nhgis0012_csv/` - This subdirectory contains the Zip code
    population data and codebook from IPUMS NHGIS.

-   Data from the Data Axle Historical Business Database
    (`BUSINESS_HISTORICAL_YYYY.zip` and `YYYY_Business_Academic_QCQ.txt.gz`)
    could not be included in the replication package due to restrictions
    on redistribution.

#### Data Processing

The following Python and Stata scripts were used to process the raw
data and produce the included analysis sample `clubstore_county.csv`:

1.   `extract.py` - Extracts records related to wholesale club stores
     from the Data Axle Historical Business Database files.  Produces
     `extract.dta`.

2.   `population_zip.do` - Generate zip-code population data in Stata.
     Produces `population_zip.dta`.

3.   `population_county.do` - Generate county population data in
     Stata.  Produces `zip_county_crosswalk.dta` and
     `population_county.dta`.

4.   `clubstore_county.do` - Generate final analysis dataset using Stata.
     Produces the analysis sample `clubstore_county.csv` and
     `ptrans.txt` (both included).

#### Analysis Sample

-   `clubstore_county.csv` - This file contains the county-level
    analysis sample used in the application.  This file was produced
    following the steps in the previous Data Processing section.

-   `ptrans.txt` - This file contains the state-to-state transition
    matrix for market size used in the application.

#### Estimation

-   `main.m` - This is the main Matlab control file for the application.
    This file produces Tables 7-11 in the paper.
-   `Gfunc.m`: Evaluates the equilibrium conditions of the dynamic entry/exit
    oligopoly model.
-   `epldygam.m`: Estimates the structural parameters using EPL.
-   `npldygam.m`: Estimates the structural parameters using NPL.
-   `freqprob.m`: Frequency probability estimation.
-   `clogit.m`: Maximum Likelihood estimation of McFadden's Conditional
    Logit model.
-   `milogit.m`: Maximum Likelihood Estimation of multinomial logit model.
-   `loglogit.m`: Log likelihood and gradient functions for multinomial logit model.
-   `mpestat.m`: Sample summary statistics.
-   `forwardsim.m`: Simulates data on state transition and decisions, used
    for counterfactual simulations.

### Monte Carlo Experiments (Appendix B.1)

Below is a complete list of program and function files used in the
Monte Carlo simulations reported in Appendix B.1 of the paper.  These
files are contained in the `mc-psd2008` subdirectory.

Main program:

-   `main.m` - Main control program that executes all Monte Carlo
    experiments for multiple sample sizes (`nobs`), multiple
    equilibria (`eqm_dgp`), and multiple levels of noise in the
    estimates (`c`).  Produces log files named according to the
    pattern `psd2008_estimate_eqmN.log` where `N` is the equilibrium
    index (`eqm_dgp` which can be 1, 2, or 3).  Results are also
    stored in `.mat` files for further processing if desired.

Main support functions:

-   `psd2008_estimate.m` - Estimates the model using the generated
    data and stores results.  Requires setting `nsims` (= 1000) before
    running.  Produces

-   `generate_data.m` - Solves for equilibria and generates data for
    the Pesendorfer and Schmidt-Dengler (2008) model.  Requires
    setting `nsims` (number of simulations; `1000` in the paper)
    `nobs` (number of observations; `250` and `1000` in the paper),
    and `eqm_dgp` (equilibrium number; 1, 2, or 3) in the Matlab
    session before running directly.  Produces `psd2008_data.mat`.

Other support functions:

-   `EqmDiff.m` - Equilibrium conditions for the model.
-   `Gderiv.m` - Derivative of G function with respect to v.
-   `LL_EPL.m` - Log likelihood for EPL estimation.
-   `LL_NPL.m` - Log likelihood for NPL estimation.
-   `normsurplus.m` - Social surplus function under the normal
    distribution.

### Monte Carlo Experiments (Appendix B.2)

Below is a complete list of program and function files used in the
Monte Carlo simulations reported in Appendix B.2 of the paper.  These
files are contained in the `mc-psd2010` subdirectory.

-   `main.m` - Main Monte Carlo simulation and estimation program.
-   `LL_MLE.m` - Log likelihood function for MLE.
-   `LL_EPL.m` - Log likelihood function for EPL.
-   `LL_NPL.m` - Log likelihood function for NPL.
-   `myCDF.m` - Error CDF specified in Pesendorfer and Schmidt-Dengler (2010).

To replicate the Monte Carlo experiments, run `main` in Matlab.

## List of Tables and Figures

-   Tables 1-6 were produced by `mc-am2007/mc_main.m`.

-   Tables 7-11 in the paper were produced by `application/main.m`.

-   Tables 12-16 were produced by `mc-psd2008/psd2008_estimate_all.m`.
    This script loops over all equilibria and sample sizes and
    estimates the model under several scenarios:

    -   Table 12 corresponds to the case where `eqm_dgp = 1`,
        `c = 0.0`, `nstart = 1`, and `nobs = [ 250, 1000 ]`.

    -   Table 13 corresponds to the case where `eqm_dgp = 2`,
        `c = 0.0`, `nstart = 1`, and `nobs = [ 250, 1000 ]`.

    -   Table 14 corresponds to the case where `eqm_dgp = 3`,
        `c = 0.0`, `nstart = 1`, and `nobs = 1000`.

    -   Table 15 corresponds to the case where `eqm_dgp = [ 1, 2 ]`,
        `c = 0.5`, `nstart = 1`, and `nobs = 250`.

    -   Table 16 corresponds to the case where `eqm_dgp = [ 1, 2 ]`,
        `c = 1.0`, `nstart = 5`, and `nobs = [ 250, 1000 ]`.

-   Table 17 was produced by `mc-psd2010/main.m`.

-   Figure 1 was produced by `mc-am2007/parameter_histograms.m`.
    The subfigures shown correspond to the following PDF files:

    -   `mc_epl_exper_1_6400_obs_param_7_hist.pdf`
    -   `mc_epl_exper_2_6400_obs_param_7_hist.pdf`
    -   `mc_epl_exper_3_6400_obs_param_7_hist.pdf`
    -   `mc_epl_exper_1_6400_obs_param_8_hist.pdf`
    -   `mc_epl_exper_2_6400_obs_param_8_hist.pdf`
    -   `mc_epl_exper_3_6400_obs_param_8_hist.pdf`

-   Figure 2 was produced by `mc-am2007/iter_histograms.m` and
    `mc-am2007/time_histograms.m`.  The subfigures shown correspond to
    the following PDF files:

    -   `mc_epl_exper_1_6400_obs_iter_hist.pdf`
    -   `mc_epl_exper_2_6400_obs_iter_hist.pdf`
    -   `mc_epl_exper_3_6400_obs_iter_hist.pdf`
    -   `mc_epl_exper_1_6400_obs_time_hist.pdf`
    -   `mc_epl_exper_2_6400_obs_time_hist.pdf`
    -   `mc_epl_exper_3_6400_obs_time_hist.pdf`

## Results

### Monte Carlo Experiments (Section 4)

The `mc-am2007/results` subdirectory contains the output of all Monte Carlo
experiments reported in Section 4 of the paper. Equilibria for each
experiment 1-3 are computed by the `mc_epl.m` program and stored in
the files `mc_epl_eq1.mat`, `mc_epl_eq2.mat`, and `mc_epl_eq3.mat`
respectively, if they do not already exist.  If they exist, the
equilibrium choice probabilities are loaded reused for computational
efficiency.

Each of the log files named `mc_epl_exper_J_N_obs.log`, contain the results of
1000 replications single Monte Carlo replication for experiment `J` (1, 2, or 3)
and sample size `N` (1600 or 6400). These log files are large, so they have
been gzip compressed to save space. The raw Matlab data files containing the
parameter estimates, number of iterations, and runtimes from each replication
are saved in .mat files named `mc_epl_exper_J_N_obs.mat`. The parameter, time,
and iteration histograms are saved in PDF files with similar names.

### Application (Section 5)

The results in Section 5 of the paper were produced by running
`application/main.m` described above using Matlab R2022b.  The log
file and results are included in the replication package:

-   `application/clubstore.log` - Contains the estimation log file for
    the application, including estimation with the observed sample
    (replication 1) as well as 250 bootstrap replications
    (replications 2-251).

-   `application/clubstore.mat` - The raw Matlab data file for the
    application, including estimation with the observed sample and all
    bootstrap replications

### Monte Carlo Experiments (Appendix B.1)

The results in this section were produced by running `mc-psd2008/main.m`
described above using Matlab 2018a.  The log files are included in the
replication package in the `mc-psd2008/results` subdirectory:

-   `psd2008_estimate_eqm1.log`
-   `psd2008_estimate_eqm2.log`
-   `psd2008_estimate_eqm3.log`

### Monte Carlo Experiments (Appendix B.2)

The results in the paper were produced by running `mc-psd2010/main.m`
described above.  A log file produced by Matlab 2018b is included in
the replication package: `mc-psd2010/mc-psd2010.log`.

## References

-   Aguirregabiria, V., and P. Mira (2007).
    [Sequential Estimation of Dynamic Discrete Games](https://doi.org/10.1111/j.1468-0262.2007.00731.x).
    _Econometrica_ 75, 1-53.

-   Data Axle. U.S. Historical Businesses, 1997-2021 [Data set].
    Retrieved November 15, 2022 from The Ohio State University Libraries
    Research Commons.
    <https://library.ohio-state.edu/record=e1002559~S7>

-   Manson, S., J. Schroeder, D. Van Riper, T. Kugler, and S. Ruggles.
    IPUMS National Historical Geographic Information System: Version 17.0
    [Data set]. Minneapolis, MN: IPUMS. 2022.
    Retrieved March 9, 2023 from IPUMS NHGIS.
    <http://doi.org/10.18128/D050.V17.0>

-   Pesendorfer, M. and P. Schmidt-Dengler (2008).
    [Asymptotic least squares estimators for dynamic games](https://doi.org/10.1111/j.1467-937X.2008.00496.x).
    _Review of Economic Studies_ 75, 901-928.

-   Pesendorfer, M. and P. Schmidt-Dengler (2010).
    [Sequential estimation of dynamic discrete games: A comment](https://doi.org/10.3982/ECTA7633).
    _Econometrica_ 78, 833-842.

-   U.S. Department of Housing and Urban Development (HUD).
    HUD-USPS ZIP Code Crosswalk data 2021-Q4 [Data set].
    Retrieved May 16, 2023 from the HUD Office of Policy Development and Research (PD&R).
    <https://www.huduser.gov/portal/datasets/usps_crosswalk.html>
