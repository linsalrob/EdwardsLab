#!/usr/bin/python

# Download WHO geographic distribution of COVID-19 cases worldwide
# Source : European Centre for Disease Prevention and Control
# Plot cases and deaths for selected countries
# The downloaded spreadsheet is stored locally in Covid-19.csv
# To use cached local spreadsheet, use "-l" option
# Intermediate data for cases/deaths and also for each country are stored in relevant .csv files
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
import os
import urllib.request

from bs4 import BeautifulSoup


topLevelPage  = "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
localFileName = "Covid-19.csv"

# countries = ["China","Germany","Italy","United_Kingdom","United_States_of_America", "Australia", "Spain",   "Sweden"]
# colours   = ["red",  "black",  "green","blue",          "orange",                   "pink",      "Violet",  "grey"]
countries = ["China","Germany","Italy","United_Kingdom","United_States_of_America", "Australia"]
colours   = ["red",  "black",  "green","blue",          "orange",                   "pink"]
country_single = ["United_Kingdom"]             # Default value, can be overwritten


# Extract cases and deaths and align day 0 to first date of detection or death
def extractCountries(covidData, country, dates, noAlignFlag):
    print("Country: " + country)

    countryData  = pd.DataFrame(index = dates)      # Create dataframe for country data

                    # Extract the data for the required country
                    # We need to copy it to the new countryData array so that dates are pre-pended back to 2019-12-31
    countryData_tmp = covidData[covidData["countriesAndTerritories"].str.match(country)]
    countryData_tmp = countryData_tmp.iloc[::-1]    # Invert table - top to bottom
    countryData[list(countryData_tmp.columns.values)] = countryData_tmp[list(countryData_tmp.columns.values)]
    countryData=countryData.fillna(0)               # Replace NaN with 0

#    countryFileName = country + '.csv'
#    # countryData.to_csv (countryFileName, index = False, header=True)
#    countryData.to_csv (countryFileName, index = True, header=True)

                                    # Fill columns : countriesAndTerritories geoId countryterritoryCode  popData2019
    countryData['countriesAndTerritories'] = countryData['countriesAndTerritories'].iloc[-1]
    countryData['geoId']                   = countryData['geoId'].iloc[-1]
    countryData['countryterritoryCode']    = countryData['countryterritoryCode'].iloc[-1]
    countryData['popData2019']             = countryData['popData2019'].iloc[-1]

                                    # Replace NaN with 0
    countryData=countryData.fillna(0)

                                    # Create cumulative cases column and cumulative deaths column - Rename column titles
    countryDataCS = countryData.cumsum(axis = 0)
    countryDataCS = countryDataCS.rename(columns={"cases": "casesCumulative", "deaths": "deathsCumulative"})

                                    # Copy cumulative columns to countryData
    countryData['casesCumulative']  = countryDataCS['casesCumulative']
    countryData['deathsCumulative'] = countryDataCS['deathsCumulative']

                                    # Calculate moving averages
    countryData['casesMA'] = countryData['cases'].rolling(7).mean()
    countryData['deathsMA'] = countryData['deaths'].rolling(7).mean()
    countryData['casesCumulativeMA'] = countryData['casesCumulative'].rolling(7).mean()
    countryData['deathsCumulativeMA'] = countryData['deathsCumulative'].rolling(7).mean()

                                    # Calculate fatality rates and clip to 100%
    countryData['fatalityPercentage'] = countryData['deaths'] * 100./countryData['cases']
    countryData['fatalityPercentage'] = countryData['fatalityPercentage'].where(countryData['fatalityPercentage'] <= 100., 100.)
    countryData.loc[countryData.cases == 0, "fatalityPercentage"] = 0                           # When cases 0= 0 set percentage to 0

    countryData['fatalityPercentageCumulative'] = countryData['deathsCumulative'] * 100./countryData['casesCumulative']
    countryData['fatalityPercentageCumulative'] = countryData['fatalityPercentageCumulative'].where(countryData['fatalityPercentageCumulative'] <= 100., 100.)
    countryData.loc[countryData.casesCumulative == 0, "fatalityPercentageCumulative"] = 0       # When cases 0= 0 set percentage to 0

    countryData['fatalityPercentageMA'] = countryData['fatalityPercentage'].rolling(7).mean()
    countryData['fatalityPercentageCumulativeMA'] = countryData['fatalityPercentageCumulative'].rolling(7).mean()

    outputFileName = country + ".csv"
    countryData.to_csv(outputFileName, index=True)

                                    # Print latest cases and deaths count
    latestCases=countryData['cases'].iloc[-1]
    latestDeaths=countryData['deaths'].iloc[-1]
    print('Latest cases          : ' + str(latestCases))
    print('Latest deaths         : ' + str(latestDeaths))
    print('Latest fatality rate  : ' + str((latestDeaths*100.)/latestCases))

                                    # Print first data of cases
    dc = countryData.index[countryData['cases'] != 0].tolist()
    print("First Case            : " + str(dc[0]).replace(' 00:00:00',''))

                                    # Print first data of deaths
    dd = countryData.index[countryData['deaths'] != 0].tolist()
    print("First Death           : " + str(dd[0]).replace(' 00:00:00',''))

    totalCases=countryData['casesCumulative'].iloc[-1]
    totalDeaths=countryData['deathsCumulative'].iloc[-1]
    fatalityRate=totalDeaths*100./totalCases
    population=countryData['popData2019'].iloc[0]
    print('Total number of cases : ' + str(totalCases))
    print('Total number of deaths: ' + str(totalDeaths))
    print("Fatality rate         : %.2f %%" % (fatalityRate))
    print("Population (2018)     : %.2f (Million)" % (population / 1000000.))
    print('')

                    # If we are not aligning first case or death then just return the data
    if noAlignFlag == True:
        return country, population, countryData, countryData;

                    # Align to first case or death by removing leading zeros and resetting index
                    # Get names of indexes for which column casesCumulative has value 0
    else:
                                    # Remove leading zeros from cases
        indexNames = countryData[ countryData['casesCumulative'] == 0 ].index
        extractedCases = countryData.drop(indexNames)
        extractedCases = extractedCases.reset_index()

                                    # Remove leading zeros from deaths
        indexNames = countryData[ countryData['deathsCumulative'] == 0 ].index
        extractedDeaths = countryData.drop(indexNames)
        extractedDeaths = extractedDeaths.reset_index()

        return country, population, extractedCases, extractedDeaths;


# main
def main(useCachedFileFlag, cumulativeResultsFlag, noAlignFlag, noPlotFlag, fileSavePlotFlag, popNormalizeFlag):
    global countries

    cachedFilePresentFlag = False

    try:
        f = open(localFileName)
        cachedFilePresentFlag = True
        f.close()
    except IOError:
        cachedFilePresentFlag = False


                                    # If cached file not present or we have requested to refresh then get the file
    if (cachedFilePresentFlag == False) or (useCachedFileFlag == False):
        cachedFilePresentFlag = False

        resp = urllib.request.urlopen(topLevelPage)
        soup = BeautifulSoup(resp, "html.parser", from_encoding=resp.info().get_param('charset'))

        for link in soup.find_all('a', href=True):
            # print(link['href'])
            if ("xlsx" in link['href']):                # If data in .xlsx format then retrieve and store as local .csv format
                xlsxfileurl = link['href']
                try:
                    xlsx_tmp = pd.read_excel(xlsxfileurl, index_col=0)
                    xlsx_tmp.to_csv(localFileName, index=True)
                    cachedFilePresentFlag = True
                    print("Cached spreadheet updated (xlsx)")
                except:
                    cachedFilePresentFlag = False
                    print("Spreadheet failed to download (xlsx)")
                break

        if (cachedFilePresentFlag == False):            # If data NOT in .xlsx format then retrieve and store .csv format
            for link in soup.find_all(class_="btn btn-primary", href=True):
                if ("csv" in link['href']):
                    csvfileurl = link['href']
                    try:
                        urllib.request.urlretrieve(csvfileurl, localFileName)
                        cachedFilePresentFlag = True
                        print("Cached spreadheet updated (csv)")
                    except:
                        cachedFilePresentFlag = False
                        print("Spreadheet failed to download (csv)")
                    break

        if (cachedFilePresentFlag == False):
            print("No spreadsheet found at the URL, use \"-l\" to use local cached file")
            exit(0)


    numberOfCountries = len(countries)

    extractedCountry = {}           # Create empty dictionaries to store result data frames for each country
    extractedPopulation = {}
    extractedCases = {}
    extractedDeaths = {}

    if (cachedFilePresentFlag == True):
        covidData = pd.read_csv(localFileName, index_col=0, encoding="iso8859_1")
                    # Spreadsheet columns :
                    # dateRep	day	month	year	cases	deaths	countriesAndTerritories	geoId	countryterritoryCode	popData2019

        covidData=covidData.fillna(0)               # Replace NaN with 0

        clen = 0                                    # For longest sequency
        dlen = 0

                                                    # Extract Chinese dates to create index - this allows for countries that do not have full data supplied
        dates_tmp = covidData[covidData["countriesAndTerritories"].str.match("China")]
        dates_tmp = dates_tmp.iloc[::-1]            # Invert table - top to bottom
        dates_tmp=dates_tmp.reset_index()
        dates=list(dates_tmp['dateRep'])

        countryIndex = 0
        for country in countries:
                                                    # Extract the data for each country
                                                    # Data can be aligned on 2019-12-29 or first instance
            extractedCountry[countryIndex], extractedPopulation[countryIndex], extractedCases[countryIndex], extractedDeaths[countryIndex] = \
                extractCountries(covidData, country, dates, noAlignFlag)

            # print(extractedCases)
            # print(extractedDeaths)

            if (popNormalizeFlag == True):
                extractedCases[countryIndex]['cases']               = extractedCases[countryIndex]['cases']                 * (10000000.0 / extractedPopulation[countryIndex])
                extractedCases[countryIndex]['deaths']              = extractedCases[countryIndex]['deaths']                * (10000000.0 / extractedPopulation[countryIndex])
                extractedCases[countryIndex]['casesCumulative']     = extractedCases[countryIndex]['casesCumulative']       * (10000000.0 / extractedPopulation[countryIndex])
                extractedCases[countryIndex]['deathsCumulative']    = extractedCases[countryIndex]['deathsCumulative']      * (10000000.0 / extractedPopulation[countryIndex])
                extractedCases[countryIndex]['casesMA']             = extractedCases[countryIndex]['casesMA']               * (10000000.0 / extractedPopulation[countryIndex])
                extractedCases[countryIndex]['deathsMA']            = extractedCases[countryIndex]['deathsMA']              * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['cases']              = extractedDeaths[countryIndex]['cases']                * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['deaths']             = extractedDeaths[countryIndex]['deaths']               * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['casesCumulative']    = extractedDeaths[countryIndex]['casesCumulative']      * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['deathsCumulative']   = extractedDeaths[countryIndex]['deathsCumulative']     * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['casesCumulativeMA']  = extractedDeaths[countryIndex]['casesCumulativeMA']    * (10000000.0 / extractedPopulation[countryIndex])
                extractedDeaths[countryIndex]['deathsCumulativeMA'] = extractedDeaths[countryIndex]['deathsCumulativeMA']   * (10000000.0 / extractedPopulation[countryIndex])

            clen = np.maximum(clen, extractedCases[countryIndex].shape[0])
            dlen = np.maximum(dlen, extractedDeaths[countryIndex].shape[0])

            countryIndex = countryIndex+1

                                                    # Select daily or cumulative results
        if (cumulativeResultsFlag == True):
            casesType        = 'casesCumulative'
            deathsType       = 'deathsCumulative'
            casesMAType      = 'casesCumulativeMA'
            deathsMAType     = 'deathsCumulativeMA'
            percentageType   = 'fatalityPercentageCumulative'
            percentageMAType = 'fatalityPercentageCumulativeMA'
        else:
            casesType        = 'cases'
            deathsType       = 'deaths'
            casesMAType      = 'casesMA'
            deathsMAType     = 'deathsMA'
            percentageType   = 'fatalityPercentage'
            percentageMAType = 'fatalityPercentageMA'

                                                    # Get last date in combinedCases
        lastDate = str(covidData.first_valid_index())
        lastDate = lastDate.replace(' 00:00:00','')


        if len(countries) == 1:                     # Single country - Cases And Deaths
                                                    # Plot titles
            titleStr='Covid-19 '
            if (cumulativeResultsFlag == True):
                titleStr=titleStr + ' Cumulative Cases And Deaths: '
            else:
                titleStr=titleStr + ' Daily Cases And Deaths: '

            ax = plt.gca()                          # Create plot - get current axis
            ax.autoscale(enable=True, tight=True)

                                                    # Plot daily cases and deaths AND moving average
            extractedCases[0].plot(kind='line', y=casesType,  label='Cases', color='blue',ax=ax)
            extractedCases[0].plot(kind='line', y=casesMAType,  label='Cases - 7 Day Moving Average', color='dodgerblue',ax=ax)
            extractedDeaths[0].plot(kind='line',y=deathsType, label='Deaths',color='seagreen',ax=ax)
            extractedDeaths[0].plot(kind='line',y=deathsMAType, label='Deaths - 7 Day Moving Average',color='lime',ax=ax)

                                                    # Plot daily mortality rate
            ax2 = plt.gca().twinx()
            extractedDeaths[0].plot(kind='line',y=percentageType,label='Mortality Rate (%)',color='red',linewidth=.75,ax=ax2)
            extractedDeaths[0].plot(kind='line',y=percentageMAType,label='Mortality Rate Moving Average (%)',color='orange',linewidth=.75,ax=ax2)
            ax2.set_ylabel('Mortality Rate (%)', color='red')
            ax2.tick_params(axis='y', labelcolor='red')
            ax2.set_ylim(0, 40)
            ax2.get_legend().remove()

            plt.title(extractedCountry[0] + '\n' + titleStr + str(lastDate) + "\nSource: European Centre for Disease Prevention and Control")

            if (fileSavePlotFlag == True):
                titleStr = titleStr.split(':', 1)[0]
                plt.savefig(titleStr.replace(" ", "")+extractedCountry[0]+'.png')

            if noPlotFlag == False:                 # Plot the data
                plt.show()
            else:
                plt.close()

        else:                                       # Multiple countries
                                                    # Plot titles
            titleStr='Covid-19 '
            if (noAlignFlag == False):
                titleStr=titleStr + 'Aligned '
            if (popNormalizeFlag == True):
                titleStr=titleStr + 'Pop Noramlized '
            if (cumulativeResultsFlag == True):
                titleStr=titleStr + 'Cumulative Cases: '
            else:
                titleStr=titleStr + 'Daily Cases: '

            ax = plt.gca()                          # Create plot - get current axis
            countryIndex = 0
            for country in countries:
                extractedCases[countryIndex].plot(kind='line',y=casesType,title=titleStr + str(lastDate) + "\nSource: European Centre for Disease Prevention and Control",label=extractedCountry[countryIndex],color=colours[countryIndex],ax=ax)
                countryIndex = countryIndex+1

            if (fileSavePlotFlag == True):
                titleStr = titleStr.split(':', 1)[0]
                plt.savefig(titleStr.replace(" ", "")+'.png')

            if noPlotFlag == False:                 # Plot the data
                plt.show()
            else:
                plt.close()

                                                    # Plot titles
            titleStr='Covid-19 '
            if (noAlignFlag == False):
                titleStr=titleStr + 'Aligned '
            if (popNormalizeFlag == True):
                titleStr=titleStr + 'Pop Noramlized '
            if (cumulativeResultsFlag == True):
                titleStr=titleStr + 'Cumulative Deaths: '
            else:
                titleStr=titleStr + 'Daily Deaths: '

            ax = plt.gca()                          # Create plot - get current axis
            countryIndex = 0
            for country in countries:
                extractedDeaths[countryIndex].plot(kind='line',y=deathsType,title=titleStr + str(lastDate) + "\nSource: European Centre for Disease Prevention and Control",label=extractedCountry[countryIndex],color=colours[countryIndex],ax=ax)
                countryIndex = countryIndex+1

            if (fileSavePlotFlag == True):
                titleStr = titleStr.split(':', 1)[0]
                plt.savefig(titleStr.replace(" ", "")+'.png')

            if noPlotFlag == False:                 # Plot the data
                plt.show()
            else:
                plt.close()

                                                    # Plot titles
            titleStr='Covid-19 '
            if (noAlignFlag == False):
                titleStr=titleStr + 'Aligned '
            if (popNormalizeFlag == True):
                titleStr=titleStr + 'Pop Noramlized '
            if (cumulativeResultsFlag == True):
                titleStr=titleStr + 'Cumulative Fatality Percentage: '
            else:
                titleStr=titleStr + 'Daily Fatality Percentage: '

            ax = plt.gca()                          # Create plot - get current axis
            countryIndex = 0
            for country in countries:
                extractedDeaths[countryIndex].plot(kind='line',y=percentageType,title=titleStr + str(lastDate) + "\nSource: European Centre for Disease Prevention and Control",label=extractedCountry[countryIndex],color=colours[countryIndex],ax=ax)
                countryIndex = countryIndex+1

            if (fileSavePlotFlag == True):
                titleStr = titleStr.split(':', 1)[0]
                plt.savefig(titleStr.replace(" ", "")+'.png')

            if noPlotFlag == False:                 # Plot the data
                plt.show()
            else:
                plt.close()

    else:
        print("Cached spreadsheet file not found on computer")
        exit()


if __name__ == '__main__':
    useCachedFileFlag       = False
    cumulativeResultsFlag   = False
    noAlignFlag             = False
    noPlotFlag              = False
    fileSavePlotFlag        = False
    popNormalizeFlag        = False

    if ((len(countries) > 1) and (len(countries) != len(colours))):
        print("The number of colours must equal the number of countries")
        exit()


    parser = argparse.ArgumentParser(description='Covid-19 Visualizer')
    parser.add_argument("-c", "--cumulative",     action="store_true", help="Display cumulative results")
    parser.add_argument("-l", "--local",          action="store_true", help="Use local cached Covid-19.csv")
    parser.add_argument("-n", "--noalign",        action="store_true", help="Do not align first instance dates - all graphs start 2019-12-31")
    parser.add_argument("-f", "--file",           action="store_true", help="Save plot to file")
    parser.add_argument("-p", "--population",     action="store_true", help="Use population to normalize data to cases per 10 Million")
    parser.add_argument("-q", "--quiet",          action="store_true", help="Quiet - Do not plot graphs")
    parser.add_argument("-s", "--single",         nargs='?', const=1,  help="Process a single country - Specify the countriesAndTerritories string used in the spreadsheet")
    parser.add_argument("-m", "--ma",             action="store_true", help="Plot Moving Average (only availeble for single courntry plot)")
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

    if (args.file):
        fileSavePlotFlag = True
        print("Save plot graphs to file = True")

    if (args.population):
        popNormalizeFlag = True
        print("Normalize to population = True")

    if (args.quiet):
        noPlotFlag = True
        print("Do not plot graphs = True")

    if (args.single):                   # Process single country - if no country specified use default country at top of file
        if (args.single != 1):
            country_single[0] = args.single

        countries = country_single      # Overwrite the countries array
        noAlignFlag = True
        print("Process single country: " + str(countries))
        print("Do not align first instance date = True")

    main(useCachedFileFlag, cumulativeResultsFlag, noAlignFlag, noPlotFlag, fileSavePlotFlag, popNormalizeFlag)
