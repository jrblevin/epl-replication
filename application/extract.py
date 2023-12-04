#!/usr/bin/env python3

# "COMPANY","ADDRESS LINE 1","CITY","STATE","ZIPCODE","ZIP4","COUNTY CODE","AREA CODE","IDCODE","LOCATION EMPLOYEE SIZE CODE","LOCATION SALES VOLUME CODE","PRIMARY SIC CODE","SIC6_DESCRIPTIONS","PRIMARY NAICS CODE","NAICS8 DESCRIPTIONS","SIC CODE","SIC6_DESCRIPTIONS (SIC)","SIC CODE 1","SIC6_DESCRIPTIONS (SIC1)","SIC CODE 2","SIC6_DESCRIPTIONS(SIC2)","SIC CODE 3","SIC6_DESCRIPTIONS(SIC3)","SIC CODE 4","SIC6_DESCRIPTIONS(SIC4)","ARCHIVE VERSION YEAR","YELLOW PAGE CODE","EMPLOYEE SIZE (5) - LOCATION","SALES VOLUME (9) - LOCATION","BUSINESS STATUS CODE","INDUSTRY SPECIFIC FIRST BYTE","YEAR ESTABLISHED","OFFICE SIZE CODE","COMPANY HOLDING STATUS","ABI","SUBSIDIARY NUMBER","PARENT NUMBER","PARENT ACTUAL EMPLOYEE SIZE","PARENT ACTUAL SALES VOLUME","PARENT EMPLOYEE SIZE CODE","PARENT SALES VOLUME CODE","SITE NUMBER","ADDRESS TYPE INDICATOR","POPULATION CODE","CENSUS TRACT","CENSUS BLOCK","LATITUDE","LONGITUDE","MATCH CODE","CBSA CODE","CBSA LEVEL","CSA CODE","FIPS CODE"

# Company, Address Line 1, City, State, ZipCode, Zip4, County Code, Area Code, IDCode, Location Employee Size Code, Location Sales Volume Code, Primary SIC Code, SIC6_Descriptions, Primary NAICS Code, NAICS8 Descriptions, SIC Code, SIC6_Descriptions (SIC), SIC Code 1, SIC6_Descriptions (SIC1), SIC Code 2, SIC6_Descriptions(SIC2), SIC Code 3, SIC6_Descriptions(SIC3), SIC Code 4, SIC6_Descriptions(SIC4), Archive Version Year, Yellow Page Code, Employee Size (5) - Location, Sales Volume (9) - Location, Business Status Code, Industry Specific First Byte, Year Established, Office Size Code, Company Holding Status, ABI, Subsidiary Number, Parent Number, Parent Actual Employee Size, Parent Actual Sales Volume, Parent Employee Size Code, Parent Sales Volume Code, Site Number, Address Type Indicator, Population Code, Census Tract, Census Block, Latitude, Longitude, Match Code, CBSA Code, CBSA Level, CSA Code, FIPS Code

# Company, Address Line 1, City, State, ZipCode, Zip4, County Code, Area Code, IDCode, Location Employee Size Code, Location Sales Volume Code, Primary SIC Code, SIC6_Descriptions, Primary NAICS Code, NAICS8 Descriptions, SIC Code, SIC6_Descriptions (SIC), SIC Code 1, SIC6_Descriptions (SIC1), SIC Code 2, SIC6_Descriptions(SIC2), SIC Code 3, SIC6_Descriptions(SIC3), SIC Code 4, SIC6_Descriptions(SIC4), Archive Version Year, Yellow Page Code, Employee Size (5) - Location, Sales Volume (9) - Location, Business Status Code, Industry Specific First Byte, Year Established, Office Size Code, Company Holding Status, ABI, Subsidiary Number, Parent Number, Parent Actual Employee Size, Parent Actual Sales Volume, Parent Employee Size Code, Parent Sales Volume Code, Site Number, Address Type Indicator, Population Code, Census Tract, Census Block, Latitude, Longitude, Match Code, CBSA Code, CBSA Level, CSA Code, FIPS Code

import zipfile
import gzip
import pandas as pd

FILES = [
    'BUSINESS_HISTORICAL_1997.zip',
    'BUSINESS_HISTORICAL_1998.zip',
    'BUSINESS_HISTORICAL_1999.zip',
    'BUSINESS_HISTORICAL_2000.zip',
    'BUSINESS_HISTORICAL_2001.zip',
    'BUSINESS_HISTORICAL_2002.zip',
    'BUSINESS_HISTORICAL_2003.zip',
    'BUSINESS_HISTORICAL_2004.zip',
    'BUSINESS_HISTORICAL_2005.zip',
    'BUSINESS_HISTORICAL_2006.zip',
    'BUSINESS_HISTORICAL_2007.zip',
    'BUSINESS_HISTORICAL_2008.zip',
    'BUSINESS_HISTORICAL_2009.zip',
    'BUSINESS_HISTORICAL_2010.zip',
    'BUSINESS_HISTORICAL_2011.zip',
    'BUSINESS_HISTORICAL_2012.zip',
    'BUSINESS_HISTORICAL_2013.zip',
    'BUSINESS_HISTORICAL_2014.zip',
    'BUSINESS_HISTORICAL_2015.zip',
    'BUSINESS_HISTORICAL_2016.zip',
    'BUSINESS_HISTORICAL_2017.zip',
    'BUSINESS_HISTORICAL_2018.zip',
    'BUSINESS_HISTORICAL_2019.zip',
    'BUSINESS_HISTORICAL_2020.zip',
    '2021_Business_Academic_QCQ.txt.gz'
]

FIELDS = [ "company", "address1", "city", "state", "zip_code", "zip4", "county_code", "area_code", "idcode", "employee_size_code", "sales_volume_code", "primary_sic_code", "sic6_descriptions", "primary_naics_code", "naics8_descriptions", "sic_code", "sic6_descriptions_sic", "sic_code_1", "sic6_descriptions_sic1", "sic_code_2", "sic6_descriptions_sic2", "sic_code_3", "sic6_descriptions_sic3", "sic_code_4", "sic6_descriptions_sic4", "archive_version_year", "yellow_page_code", "employee_size", "sales_volume", "business_status_code", "industry_specific_first_byte", "year_established", "office_size_code", "company_holding_status", "abi", "subsidiary_number", "parent_number", "parent_actual_employee_size", "parent_actual_sales_volume", "parent_employee_size_code", "parent_sales_volume_code", "site_number", "address_type_indicator", "population_code", "census_tract", "census_block", "latitude", "longitude", "match_code", "cbsa_code", "cbsa_level", "csa_code", "fips_code" ]

INT_FIELDS = [ "zip_code", "zip4", "area_code", "primary_sic_code", "primary_naics_code", "employee_size", "sales_volume", "business_status_code", "census_tract", "census_block" ]

KEEP_FIELDS = [
    "company",
    "address1",
    "city",
    "state",
    "zip_code",
    "zip4",
    "area_code",
    "employee_size_code",
    "sales_volume_code",
    "primary_sic_code",
    "sic6_descriptions",
    "primary_naics_code",
    "naics8_descriptions",
    "archive_version_year",
    "employee_size",
    "sales_volume",
    "business_status_code",
    "year_established",
    "abi",
    "subsidiary_number",
    "parent_number",
    "population_code",
    "census_tract",
    "census_block",
    "latitude",
    "longitude",
    "match_code",
    "cbsa_code",
    "cbsa_level",
    "csa_code",
    "fips_code"]

def filter_csv(f, nrows):
    iter_csv = pd.read_csv(f, iterator=True, chunksize=10000,
                           header=0, names=FIELDS, nrows=nrows,
                           engine='c', encoding='ISO-8859-1', on_bad_lines='warn')

    df = pd.concat([chunk[
        chunk['company'].astype(str).str.startswith("SAM'S CLUB") |
        chunk['company'].astype(str).str.startswith("B J'S WHOLESALE CLUB") |
        chunk['company'].astype(str).str.startswith("COSTCO WHOLESALE INC") |
        (chunk['primary_sic_code'].fillna(0.0).astype(int) == 531110) # Warehouse Clubs
    ] for chunk in iter_csv])

    # Convert to integers
    for field in INT_FIELDS:
        df[field] = df[field].fillna(0.0).astype(int)

    # Find out which columns are of type 'object'
    #print(df.select_dtypes(include=['object']).columns)

    # Print summary
    print(df[['abi', 'company', 'primary_sic_code', 'primary_naics_code']])

    df = df[KEEP_FIELDS]

    # Convert all fields to strings
    for field in KEEP_FIELDS:
        df[field] = df[field].fillna('.').astype(str)

    return df

def filter_zip(f):
    dfs = []
    if f.endswith('zip'):
        with zipfile.ZipFile(f) as zf:
            for name in zf.namelist():
                if name.endswith('csv'):
                    print("\n%s / %s\n" % (f, name))
                    # Handle truncated file Business-1998-NH.csv
                    nrows=None
                    if name.endswith('Business-1998-NH.csv'):
                        nrows=61634
                    df = filter_csv(zf.open(name), nrows=nrows)
                    dfs.append(df)

    elif f.endswith('txt.gz'):
        print("\n%s\n" % f)
        df = filter_csv(gzip.open(f, mode='rt'), nrows=None)
        dfs.append(df)

    return pd.concat(dfs)

dfs = []
for f in FILES:
    df = filter_zip(f)
    dfs.append(df)

df = pd.concat(dfs)

df.to_csv("extract.csv")
df.to_excel("extract.xlsx")
df.to_stata("extract.dta")

