# Dearing and Blevins (2024): Application

Dearing and Blevins (2024).
Efficient and Convergent Sequential Pseudo-Likelihood Estimation of Dynamic
Discrete Games.
_Review of Economic Studies_.

## Overview

This subdirectory contains the replication code for the application to U.S.
wholesale club store competition in Section 5 of the paper.  The simulations are
based on the dynamic game of entry and exit in Example 1 of the paper.

## Data Availability and Provenance Statements

The application in the paper uses data from several sources:

-   [IPUMS NHGIS Table CL8: Total Population (1990, 2000, 2010, 2020)][nhgis]
-   [U.S. Department of Housing and Urban Development Zip Code to County Crosswalk Files][crosswalk]
-   [Data Axle Historical Business Database][infogroup]

We have included in this replication package the NHGIS population data
extract from Table CL8: Total Population (`nhgis0012_ts_geog2010_zcta.csv`)
as well as the HUD crosswalk file (`ZIP_COUNTY_122021.xlsx`).

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

## Files

Below is a complete list of program and function files used for the
application.

### Raw Data

-   `ZIP_COUNTY_122021.xlsx` - U.S. Department of Housing and Urban
    Development Zip Code to County Crosswalk File.

-   `nhgis0012_csv/` - This subdirectory contains the Zip code
    population data and codebook from IPUMS NHGIS.

-   Data from the Data Axle Historical Business Database
    (`BUSINESS_HISTORICAL_YYYY.zip` and `YYYY_Business_Academic_QCQ.txt.gz`)
    could not be included in the replication package due to restrictions
    on redistribution.

### Data Processing

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

The Python script `extract.py` imports the following Python modules:
`zipfile`, `gzip`, and `pandas`.  Although `zipfile` and `gzip` are
part of the standard library, Pandas is an external package.
Depending on your Python installation, it can likely be installed
using the following command:

```python
pip install pandas
```

The Stata scripts also require some additional packages to be
installed: `xttrans2` and `mat2txt`.  These can be installed using the
following commands:

```stata
net install xttrans2
ssc install mat2txt
```

### Analysis Sample

The analysis sample used for the application in the paper was produced
by following the steps above using Python 3.11.6, Pandas 2.0.3, Stata
MP 16.1, xttrans 1.2.0, and mat2txt 1.1.2.

-   `clubstore_county.csv` - This file contains the county-level
    analysis sample used in the application.  This file was produced
    following the steps in the previous Data Processing section.

-   `ptrans.txt` - This file contains the state-to-state transition
    matrix for market size used in the application.

### Estimation

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

### Results

The results in the paper were produced by running `main.m` described
above using Matlab R2022b.  The log file and results are included in
the replication package:

-   `clubstore.log` - Contains the estimation log file for the
    application, including estimation with the observed sample
    (replication 1) as well as 250 bootstrap replications
    (replications 2-251).

-   `clubstore.mat` - The raw Matlab data file for the application,
    including estimation with the observed sample and all bootstrap
    replications

## References

-   Steven Manson, Jonathan Schroeder, David Van Riper, Tracy Kugler, and Steven Ruggles.
    IPUMS National Historical Geographic Information System: Version 17.0
    [Data set]. Minneapolis, MN: IPUMS. 2022.
    Retrieved March 9, 2023 from IPUMS NHGIS.
    <http://doi.org/10.18128/D050.V17.0>

-   U.S. Department of Housing and Urban Development (HUD).
    HUD-USPS ZIP Code Crosswalk data 2021-Q4 [Data set].
    Retrieved May 16, 2023 from the HUD Office of Policy Development and Research (PD&R).
    <https://www.huduser.gov/portal/datasets/usps_crosswalk.html>

-   Data Axle. U.S. Historical Businesses, 1997-2021 [Data set].
    Retrieved November 15, 2022 from The Ohio State University Libraries
    Research Commons.
    <https://library.ohio-state.edu/record=e1002559~S7>
