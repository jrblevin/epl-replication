// Stata post-processing commands

// 1. Extract from raw data using Python: extract.py
// 2. Generate zip-code population data in Stata: population_zip.do
// 3. Generate county population data in Stata: population_county.do
// 4. Generate club store dataset using Stata: clubstore_county.do (this file)

// Requires xttrans2 and mat2txt
// net install xttrans2
// ssc install mat2txt

local min_pop = 20000
local max_pop = 600000
local min_year = 2009
local max_year = 2021
local popvar "lpop"
local bintype "equal" // "percentiles" or "equal"

clear
log using clubstore_county.log, replace

use extract.dta

destring abi, replace
destring archive_version_year, generate(year)
destring census_tract, replace
destring fips, replace
destring census_block, replace
destring cbsa_code, replace
destring csa_code, replace
destring zip_code, replace
destring zip4, replace
destring area_code, replace
destring primary_naics_code, replace
destring primary_sic_code, replace
destring employee_size, replace
destring sales_volume, replace
destring cbsa_level, replace
destring longitude latitude, replace
destring subsidiary_number parent_number, replace
destring business_status_code, replace

// Fix .0 in population_code
gen pop_code = subinstr(population_code,".0","",.)
destring pop_code, replace
drop population_code

// Fix quoting in company names
gen company_new = subinstr(company,"''","'",.)
drop company
rename company_new company

// PHARMACY - 591205
drop if strpos(company, "PHARMACY") > 0
drop if primary_sic_code >= 591000 & primary_sic_code < 592000

// BAKERY - 546102
drop if strpos(company, "BAKERY") > 0
drop if primary_sic_code >= 546000 & primary_sic_code < 547000

// TIRE - 553123
drop if strpos(company, "TIRE") > 0
drop if primary_sic_code >= 553000 & primary_sic_code < 554000

// OPTICAL - 599502
drop if strpos(company, "OPTICAL") > 0
drop if primary_sic_code >= 599000 & primary_sic_code < 600000

// JEWELRY - 594409
drop if strpos(company, "JEWELRY") > 0
drop if primary_sic_code >= 594000 & primary_sic_code < 595000

drop if strpos(company, "DISTRIBUTION") > 0
drop if strpos(company, "DISTR") > 0
drop if strpos(company, "DSTR") > 0

drop if strpos(company, "VISION") > 0

drop if strpos(company, "PHOTO") > 0
drop if strpos(company, "PHOT") > 0

drop if strpos(company, "GAS") > 0
drop if strpos(company, "FUEL") > 0

drop if strpos(company, "FLORAL") > 0

drop if strpos(company, "LIQUOR") > 0

drop if strpos(company, "EXPRESS") > 0

drop if strpos(company, "CONNECTION CTR") > 0

drop if strpos(company, "SATELLITE") > 0

drop if strpos(company, "LIBRARY") > 0

// Normalize company names
gen company_new = subinstr(company," MEMBERS ONLY","",.)

// FIRM ID
gen firm_id = 5
label define firm_id_label 1 "Sam's Club" 2 "Costco" 3 "BJ's" 4 "Direct Buy" 5 "Other"
label values firm_id firm_id_label

// Sam's Club = 1
replace company_new = "Sam's Club" if strpos(company,"SAM") > 0
replace firm_id = 1 if company_new == "Sam's Club"

// Costco = 2
replace company_new = "Costco" if strpos(company,"COSTCO") > 0
replace firm_id = 2 if company_new == "Costco"

// BJ's = 3
replace company_new = "BJ's" if strpos(company,"BJ") > 0
replace company_new = "BJ's" if strpos(company,"B J") > 0
replace firm_id = 3 if company_new == "BJ's"

// Direct Buy = 4
replace company_new = "Direct Buy" if strpos(company,"DIRECTBUY") > 0
replace company_new = "Direct Buy" if strpos(company,"DIRECT BUY") > 0
replace firm_id = 4 if company_new == "Direct Buy"

// Other = 5
replace company_new = "Other" if firm_id == 5

// Drop Puerto Rico
drop if state == "PR"

// Drop if zip code is missing
drop if zip_code == .

// Merge zip-to-county crosswalk
joinby zip_code using zip_county_crosswalk.dta

// Only count store for the primary county
keep if primary_county == 1

// Merge in county population data
merge m:1 county year using population_county
drop if _merge < 3
drop _merge

// Generate pseudo-id
sort county year firm_id
gen pseudo_id = 100000 * firm_id + county
by county year firm_id: gen nstore = _N
gen active = nstore > 0
drop nstore

drop zip_code
duplicates drop pseudo_id year, force
xtset pseudo_id year, yearly
tsfill, full
gen firm_id_new = int(pseudo_id/100000)
label values firm_id_new firm_id_label
drop firm_id
rename firm_id_new firm_id
gen county_new = pseudo_id - 100000 * firm_id
drop county
rename county_new county

replace active = 0 if active == .

duplicates drop county year firm_id, force
keep county year firm_id active
order county year firm_id active

reshape wide active, i(year county) j(firm_id)
sort county year

// Merge in population data
merge 1:1 county year using population_county

// Generate market ID
sort county year
gen market = ceil(_n/25)

// Panel structure
xtset market year, yearly
summarize market year

// Activity indicators
replace active1 = 0 if active1 == .
replace active2 = 0 if active2 == .
replace active3 = 0 if active3 == .
replace active4 = 0 if active4 == .
replace active5 = 0 if active5 == .
gen nactive = active1 + active2 + active3 + active4 + active5

// Lagged activity indicators
gen lactive1 = l.active1
gen lactive2 = l.active2
gen lactive3 = l.active3
gen lactive4 = l.active4
gen lactive5 = l.active5

// Statistics about popluation
summarize pop, detail
sort county year
by county: egen all_time_max_pop = max(pop)
by county: egen all_time_min_pop = min(pop)

// Select markets by population
summarize all_time_min_pop all_time_max_pop if year == `min_year', detail
drop if all_time_min_pop < `min_pop'
drop if all_time_max_pop > `max_pop'
summarize all_time_max_pop if year == `min_year', detail

// Discretize population
gen lpop = log(pop)
summarize pop lpop
gen dpop = 1

if "`bintype'" == "equal" {
   local nbins = 5
   sum `popvar', meanonly
   local max = r(max)
   local min = r(min)
   local binwidth = (`max' - `min') / `nbins'
   di as txt "Discretizing population (`popvar') into equal width bins of width `binwidth'"
   replace dpop = 2 if `popvar' >= (`min' + `binwidth')
   replace dpop = 3 if `popvar' >= (`min' + 2*`binwidth')
   replace dpop = 4 if `popvar' >= (`min' + 3*`binwidth')
   replace dpop = 5 if `popvar' >= (`min' + 4*`binwidth')
}
else {
   di as txt "Discretizing population (`popvar') into percentile bins"
   _pctile `popvar', p(20 40 60 80)
   local p1 = r(r1)
   local p2 = r(r2)
   local p3 = r(r3)
   local p4 = r(r4)

   replace dpop = 2 if `popvar' >= `p1'
   replace dpop = 3 if `popvar' >= `p2'
   replace dpop = 4 if `popvar' >= `p3'
   replace dpop = 5 if `popvar' >= `p4'
}

drop pop
rename dpop pop
xttrans2 pop, freq matcell(ptrans)
mat2txt, matrix(ptrans) saving(ptrans.txt) replace

summarize

// Regenerate market ID
sort market year
replace market = ceil(_n/25)

// Statistics on number of active firms (top 5)
tab nactive if year == `min_year'
tab nactive if year == `max_year'

// Statistics on number of active firms (top 3)
drop nactive
gen nactive = active1 + active2 + active3
tab nactive if year == `min_year'
tab nactive if year == `max_year'

// Final dataset
keep market year active1 active2 active3 lactive1 lactive2 lactive3 pop
order market year active1 active2 active3 lactive1 lactive2 lactive3 pop

drop if year <= `min_year'
drop if year > `max_year'

di as txt "Minimum population: `min_pop' (level)"
di as txt "Maximum population: `max_pop' (level)"
di as txt "Year range: `min_year'-`max_year'"
di as txt "Population variable: `popvar' (log or level)"
di as txt "Bin type: `bintype' (percentiles or equal)"

xtdescribe
summarize

save clubstore_county.dta, replace
outsheet using clubstore_county.csv, comma replace

log close
