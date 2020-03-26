#!/usr/bin/python

# Download WHO geographic distribution of COVID-19 cases worldwide
# Source : European Centre for Disease Prevention and Control
# Plot cases and deaths for selected countries
# The downloaded spreadsheet is stored locally in Covid-19.xlsx
# To use cached local spreadsheet, use "-l" option
# Intermediate data for Cases/Deaths and also for each country are stored in relevant .csv files
# All plots can be aligned to :
#   First date of detection or death, in that country (default)
#   First date of detection in China, 2019-12-31 (-n)
# Data can be plotted as daily values (default) cumulative values (-c)
# Countries to plot and line colours are specified in the appropriate tables at the top of this file

# Dependencies : pandas, matplotlib, numpy, google-auth-httplib2, beautifulsoup4, xlrd


import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import urllib.request

from bs4 import BeautifulSoup


topLevelPage  = "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
localFileName = "Covid-19.xlsx"

countries = ["China","Germany","Italy","United_Kingdom","United_States_of_America", "Australia"]
colours   = ["red",  "black",  "green","blue",          "orange",                   "pink"]

# Extract cases and deaths and align day 0 to first date of detection or death
def extractAligned(covidData, country, noAlignFlag):
    print("Country: " + country)

                    # Extract the data for the required country
    countryData = covidData[covidData["Countries and territories"].str.match(country)]
    countryData = countryData.iloc[::-1]

                    # Create cumulative Cases column and cumulative Deaths column - Rename column titles
    countryDataCS = countryData.cumsum(axis = 0)
    countryDataCS = countryDataCS.rename(columns={"Cases": "Cases_Cumulative", "Deaths": "Deaths_Cumulative"})

                    # Copy cumulative columns to countryData
    countryData['Cases_Cumulative']  = countryDataCS['Cases_Cumulative']
    countryData['Deaths_Cumulative'] = countryDataCS['Deaths_Cumulative']


    outputFileName = country + ".csv"
    countryData.to_csv(outputFileName, index=False)

                    # Print first data of Cases
    dc = countryData.index[countryData['Cases'] != 0].tolist()
    print("First Case            : " + str(dc[0]).replace(' 00:00:00',''))

                    # Print first data of Deaths
    dd = countryData.index[countryData['Deaths'] != 0].tolist()
    print("First Death           : " + str(dd[0]).replace(' 00:00:00',''))

                    # Remove leading zeros from Cumulative_Cases
                    # Get names of indexes for which column Cases_Cumulative has value 0
    indexNames = countryData[ countryData['Cases_Cumulative'] == 0 ].index
                    # Delete these row indexes from dataFrame
    ecnlz = countryData.drop(indexNames)
    if noAlignFlag == True:
        ecnlz.index = pd.to_datetime(ecnlz.index)
    else:
        ecnlz = ecnlz.reset_index()

                    # Remove leading zeros from Cumulative_Deaths
                    # Get names of indexes for which column Deaths_Cumulative has value 0
    indexNames = countryData[ countryData['Deaths_Cumulative'] == 0 ].index
                    # Delete these row indexes from dataFrame
    ednlz = countryData.drop(indexNames)
    if noAlignFlag == True:
        ednlz.index = pd.to_datetime(ednlz.index)
    else:
        ednlz = ednlz.reset_index()

    totalCases=countryData['Cases_Cumulative'].iloc[-1]
    totalDeaths=countryData['Deaths_Cumulative'].iloc[-1]
    fatalityRate=totalDeaths*100./totalCases


    print('Total number of Cases : ' + str(totalCases))
    print('Total number of Deaths: ' + str(totalDeaths))
    print("Fatality rate         : %.2f %%" % (fatalityRate))
    print('')
    return country, ecnlz, ednlz;


# main
def main(useCachedFileFlag, cumulativeResultsFlag, noAlignFlag):

    cachedFilePresentFlag = False

    try:
        f = open(localFileName)
        cachedFilePresentFlag = True
        f.close()
    except IOError:
        cachedFilePresentFlag = False


                    # If cached file not present or we have requested to refresh then get the file
    if (cachedFilePresentFlag == False) or (useCachedFileFlag == False):
        resp = urllib.request.urlopen(topLevelPage)
        soup = BeautifulSoup(resp, "html.parser", from_encoding=resp.info().get_param('charset'))

        for link in soup.find_all('a', href=True):
            # print(link['href'])
            if (".xlsx" in link['href']):
                # print(link['href'])
                xlsxfileurl = link['href']

        if (xlsxfileurl):
            urllib.request.urlretrieve(xlsxfileurl, localFileName)
            cachedFilePresentFlag = True
            print("Cached file updated")
        else:
            print("Spreadsheet file not found on website")
            exit()

    numberOfCountries = len(countries)

    ecountry = {}   # Create empty dictionaries to store result data frames for each country
    ecnlz = {}
    ednlz = {}

    if (cachedFilePresentFlag == True):
        covidData = pd.read_excel(localFileName, index_col=0)
                    # Spreadsheet columns :
                    # DateRep	Day	Month	Year	Cases	Deaths	Countries and territories	GeoId

        clen = 0    # For longest sequency
        dlen = 0

        countryIndex = 0
        for country in countries:
                    # For each country - extract the aligned data
                    # Data can be aligned on 2019-12-29 or first instance
            ecountry[countryIndex], ecnlz[countryIndex], ednlz[countryIndex] = extractAligned(covidData, country, noAlignFlag)

            clen = np.maximum(clen, ecnlz[countryIndex].shape[0])
            dlen = np.maximum(dlen, ednlz[countryIndex].shape[0])

            countryIndex = countryIndex+1


                    # Create DataFrame of countries cases and deaths
        if noAlignFlag == True:
            dl=ecnlz[0].filter(['DateRep'], axis=1)         # Extract dates to create index
            dl=dl.reset_index()
            dates=list(dl['DateRep'])

            combinedCases  = pd.DataFrame(index = dates)    # Create dataframes
            combinedDeaths = pd.DataFrame(index = dates)

        else:
            c_idx = np.arange(0, clen, 1)
            d_idx = np.arange(0, dlen, 1)

            combinedCases  = pd.DataFrame(index = c_idx)    # Create dataframes
            combinedDeaths = pd.DataFrame(index = d_idx)

        if (cumulativeResultsFlag == True):         # Select daily or cumulative results
            casesType  = 'Cases_Cumulative'
            deathsType = 'Deaths_Cumulative'
        else:
            casesType  = 'Cases'
            deathsType = 'Deaths'


        countryIndex = 0
        for country in countries:
                    # Copy Cases and Deaths columns to summary DataFrame
            combinedCases [ecountry[countryIndex]] = ecnlz[countryIndex][casesType]
            combinedDeaths[ecountry[countryIndex]] = ednlz[countryIndex][deathsType]
            countryIndex = countryIndex+1


                    # Write data to .csv files
        if (cumulativeResultsFlag == True):
            combinedCases.to_csv("cumulativeCases.csv",   index=False)
            combinedDeaths.to_csv("cumulativeDeaths.csv", index=False)
        else:
            combinedCases.to_csv("dailyCases.csv",   index=False)
            combinedDeaths.to_csv("dailyDeaths.csv", index=False)


                    # Get last date in combinedCases
        lastDate = str(covidData.first_valid_index())
        lastDate = lastDate.replace(' 00:00:00','')


                    # Print columns of countries cases and deaths
        if (cumulativeResultsFlag == True):
            titleStr='Covid-19 Cumulative Cases: ' + str(lastDate)
        else:
            titleStr='Covid-19 Daily Cases: ' + str(lastDate)

        ax = plt.gca()          # Create plot - get current axis
        countryIndex = 0
        for country in countries:
            ecnlz[countryIndex].plot(kind='line',y=casesType,title=titleStr,label=ecountry[countryIndex],color=colours[countryIndex], ax=ax)
            countryIndex = countryIndex+1
        plt.show()

        if (cumulativeResultsFlag == True):
            titleStr='Covid-19 Cumulative Deaths: ' + str(lastDate)
        else:
            titleStr='Covid-19 Daily Deaths: ' + str(lastDate)

        ax = plt.gca()          # Create plot - get current axis
        countryIndex = 0
        for country in countries:
            ednlz[countryIndex].plot(kind='line',y=deathsType,title=titleStr,label=ecountry[countryIndex],color=colours[countryIndex], ax=ax)
            countryIndex = countryIndex+1
        plt.show()

    else:
        print("Cached spreadsheet file not found on computer")
        exit()


if __name__ == '__main__':
    useCachedFileFlag       = False
    cumulativeResultsFlag   = False
    noAlignFlag             = False

    if len(countries) != len(colours):
        print("The number of colours must equal the number of countries")
        exit()


    parser = argparse.ArgumentParser(description='Covid-19 Visualizer')
    parser.add_argument("-c", "--cumulative", action="store_true", help="Display cumulative results")
    parser.add_argument("-l", "--local",      action="store_true", help="Use local cached Covid-19.xlsx")
    parser.add_argument("-n", "--noalign",    action="store_true", help="Do not align first instance dates - all graphs start 2019-12-31")
    args = parser.parse_args()

    if (args.cumulative):
        cumulativeResultsFlag = True
        print("Cumulative Results = True")

    if (args.local):
        useCachedFileFlag = True
        print("Use cached file = True")

    if (args.noalign):
        noAlignFlag = True
        print("Do not align first instance date = True")

    main(useCachedFileFlag, cumulativeResultsFlag, noAlignFlag)
