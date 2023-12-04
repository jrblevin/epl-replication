// Process NHGIS zip-code-level population data

clear
log using population_zip.log, replace

insheet using nhgis0012_csv/nhgis0012_ts_geog2010_zcta.csv
keep zctaa cl8aa1990 cl8aa2000 cl8aa2010 cl8aa2020
reshape long cl8aa, i(zctaa) j(year)
rename cl8aa pop
rename zctaa zip_code
xtset zip_code year, yearly delta(10)

// convert to annual frequency with interpolated values
expand 10, generate(added)
sort zip_code year added
by zip_code year: gen new_year = year + sum(added)

drop year
rename new_year year
tsset zip_code year, yearly
replace pop = . if added == 1
drop added
by zip_code: ipolate pop year, generate(new_pop) epolate
drop pop
rename new_pop pop
drop if year > 2021
drop if year < 1997

summarize

save population_zip.dta, replace

log close
