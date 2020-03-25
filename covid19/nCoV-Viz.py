#!/usr/bin/python

# Download WHO geographic distribution of COVID-19 cases worldwide
# Plot cases and deaths for selected countries
# All plots align day 0 to first date of detection or death
# The downloaded spreadsheet is stored locally in Covid-19.xlsx
# To use cached local spreadsheet, execute : nCoV-Viz.py -l

# Dependencies : pandas, matplotlib, numpy, google-auth-httplib2, beautifulsoup4, xlrd

topLevelPage = "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
localFileName = "Covid-19.xlsx"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import urllib.request

from bs4 import BeautifulSoup

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
    countryData['Cases_Cumulative'] = countryDataCS['Cases_Cumulative']
    countryData['Deaths_Cumulative'] = countryDataCS['Deaths_Cumulative']


    outputFileName = country + ".csv"
    countryData.to_csv(outputFileName, index=False)

                    # Print first data of Cases
    dc = countryData.index[countryData['Cases'] != 0].tolist()
    print("First Case : " + str(dc[0]))

                    # Print first data of Deaths
    dd = countryData.index[countryData['Deaths'] != 0].tolist()
    print("First Death : " + str(dd[0]))

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
        soup = BeautifulSoup(resp, from_encoding=resp.info().get_param('charset'))

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

    if (cachedFilePresentFlag == True):
        covidData = pd.read_excel(localFileName, index_col=0)

                    # DateRep	Day	Month	Year	Cases	Deaths	Countries and territories	GeoId
        cn_country, cn_ecnlz, cn_ednlz = extractAligned(covidData, "China", noAlignFlag)
        de_country, de_ecnlz, de_ednlz = extractAligned(covidData, "Germany", noAlignFlag)
        it_country, it_ecnlz, it_ednlz = extractAligned(covidData, "Italy", noAlignFlag)
        uk_country, uk_ecnlz, uk_ednlz = extractAligned(covidData, "United_Kingdom", noAlignFlag)
        us_country, us_ecnlz, us_ednlz = extractAligned(covidData, "United_States_of_America", noAlignFlag)

                    # Find out which sequence is longest
        clen = cn_ecnlz.shape[0]
        dlen = cn_ednlz.shape[0]

        clen = np.maximum(clen, de_ecnlz.shape[0])
        dlen = np.maximum(dlen, de_ednlz.shape[0])

        clen = np.maximum(clen, it_ecnlz.shape[0])
        dlen = np.maximum(dlen, it_ednlz.shape[0])

        clen = np.maximum(clen, uk_ecnlz.shape[0])
        dlen = np.maximum(dlen, uk_ednlz.shape[0])

        clen = np.maximum(clen, us_ecnlz.shape[0])
        dlen = np.maximum(dlen, us_ednlz.shape[0])

                    # Create DataFrame of countries cases and deaths
        c_idx = np.arange(0, clen, 1)
        d_idx = np.arange(0, dlen, 1)

        combinedCases  = pd.DataFrame(index = c_idx, columns = [cn_country, de_country, it_country, uk_country, us_country])
        combinedDeaths = pd.DataFrame(index = d_idx, columns = [cn_country, de_country, it_country, uk_country, us_country])

        if (cumulativeResultsFlag == True):         # Select daily or cumulative results
            casesType  = 'Cases_Cumulative'
            deathsType = 'Deaths_Cumulative'
        else:
            casesType  = 'Cases'
            deathsType = 'Deaths'


                    # Copy Cases columns to summary DataFrame
        combinedCases[cn_country] = cn_ecnlz[casesType]
        combinedCases[de_country] = de_ecnlz[casesType]
        combinedCases[it_country] = it_ecnlz[casesType]
        combinedCases[uk_country] = uk_ecnlz[casesType]
        combinedCases[us_country] = us_ecnlz[casesType]
                    # Copy Deaths columns to summary DataFrame
        combinedDeaths[cn_country] = cn_ednlz[deathsType]
        combinedDeaths[de_country] = de_ednlz[deathsType]
        combinedDeaths[it_country] = it_ednlz[deathsType]
        combinedDeaths[uk_country] = uk_ednlz[deathsType]
        combinedDeaths[us_country] = us_ednlz[deathsType]

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
        cn_ecnlz.plot(kind='line',y=casesType,title=titleStr,label=cn_country,color='red', ax=ax)
        # de_ecnlz.plot(kind='line',y=casesType,title=titleStr,label=de_country,color='black', ax=ax)
        it_ecnlz.plot(kind='line',y=casesType,title=titleStr,label=it_country,color='green', ax=ax)
        uk_ecnlz.plot(kind='line',y=casesType,title=titleStr,label=uk_country,color='blue',ax=ax)
        us_ecnlz.plot(kind='line',y=casesType,title=titleStr,label=us_country,color='orange', ax=ax)
        plt.show()

        if (cumulativeResultsFlag == True):
            titleStr='Covid-19 Cumulative Deaths: ' + str(lastDate)
        else:
            titleStr='Covid-19 Daily Deaths: ' + str(lastDate)

        ax = plt.gca()          # Create plot - get current axis
        cn_ednlz.plot(kind='line',y=deathsType,title=titleStr,label=cn_country,color='red', ax=ax)
        # de_ednlz.plot(kind='line',y=deathsType,title=titleStr,label=de_country,color='black', ax=ax)
        it_ednlz.plot(kind='line',y=deathsType,title=titleStr,label=it_country,color='green', ax=ax)
        uk_ednlz.plot(kind='line',y=deathsType,title=titleStr,label=uk_country,color='blue',ax=ax)
        us_ednlz.plot(kind='line',y=deathsType,title=titleStr,label=us_country,color='orange', ax=ax)
        plt.show()

    else:
        print("Cached spreadsheet file not found on computer")
        exit()


if __name__ == '__main__':
    useCachedFileFlag       = False
    cumulativeResultsFlag   = False
    noAlignFlag             = False

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
