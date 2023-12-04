// Construct county-level population by aggregating zip code data
// Run population_zip.do first.

clear
log using population_county.log, replace

// Prepare crosswalk file
import excel "ZIP_COUNTY_122021.xlsx", sheet("SQLT0004") firstrow
rename zip zip_code
destring zip_code, replace
destring county, replace

// Construct an indictor for primary county
by zip (tot_ratio), sort: gen primary_county = (_n == _N)

// Save for later
save zip_county_crosswalk.dta, replace

// Merge county population with crosswalk
clear
use population_zip
joinby zip_code using zip_county_crosswalk

// Sum weighted population of zip codes to create county population
generate wpop = tot_ratio * pop
drop pop
collapse (sum) pop=wpop, by(county year)

// Set panel structure
xtset county year, yearly
summarize

save population_county.dta, replace

log close
