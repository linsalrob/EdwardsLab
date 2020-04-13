#!/usr/bin/python

# Download WHO geographic distribution of COVID-19 cases worldwide
# Source : European Centre for Disease Prevention and Control
# Plot bar graph of cases and deaths for top N countries (Default N=10 but can be changed on command line)
# The downloaded spreadsheet is stored locally in Covid-19.csv
# To use cached local spreadsheet, use "-l" option

# Dependencies : pandas, matplotlib, numpy, google-auth-httplib2, beautifulsoup4, xlrd ffmpeg


import argparse
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import urllib.request

from bs4 import BeautifulSoup
from IPython.display import HTML


topLevelPage  = "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
localFileName = "Covid-19.csv"


# Extract cases and deaths and align day 0 to first date of detection or death
def extractCountries(covidData, desiredColumn, dates, popNormalizeFlag):

    countryData  = pd.DataFrame(index = dates)      # Create dataframe for country data

    # Get list of all countries
    countriesList = covidData['countriesAndTerritories'].tolist()
    countriesList = list(dict.fromkeys(countriesList))  # Remove duplicates
    countriesList.sort()                            # Sort alphabetically - just in case
#    print(countriesList)


                    # Extract the data for every country
                    # We need to copy it to the new countryData array so that dates are pre-pended back to 2019-12-31
    for country in countriesList:
        countryData_tmp = covidData[covidData["countriesAndTerritories"].str.match(country)].copy()        # Copy desired country

        if (countryData_tmp.dropna().empty == True):    # Drop empty frames
            continue

        population=countryData_tmp['popData2018'].iloc[0]
        # print("Country               : " + country)
        # print("Population (2018)     : %.2f (Million)" % (population / 1000000.))

        countryData_tmp = countryData_tmp.iloc[::-1]    # Invert table - top to bottom
                                                    # Create cumulative cases column and cumulative deaths column - Rename column titles
                                                    # We need to invert table top to bottom and back again to do cumulative sum
        countryDataCS = countryData_tmp.cumsum(axis = 0)
        countryDataCS = countryDataCS.rename(columns={"cases": "casesCumulative", "deaths": "deathsCumulative"})

                                                    # Copy cumulative columns to countryData
        countryData_tmp['casesCumulative']  = countryDataCS['casesCumulative']
        countryData_tmp['deathsCumulative'] = countryDataCS['deathsCumulative']

        if ((popNormalizeFlag == True) and (population > 0)):
            countryData_tmp['cases']            = countryData_tmp['cases']            * (10000000.0 / population)
            countryData_tmp['deaths']           = countryData_tmp['deaths']           * (10000000.0 / population)
            countryData_tmp['casesCumulative']  = countryData_tmp['casesCumulative']  * (10000000.0 / population)
            countryData_tmp['deathsCumulative'] = countryData_tmp['deathsCumulative'] * (10000000.0 / population)

        countryData_tmp = countryData_tmp.loc[:, countryData_tmp.columns.intersection([desiredColumn])]   # Delete all columns except desired


                                                    # Change column name to Country
        countryData_tmp = countryData_tmp.rename(columns={desiredColumn: country})
                                                    # Remove duplicate indices
        countryData_tmp = countryData_tmp.loc[~countryData_tmp.index.duplicated(keep='first')]

        countryData[country] = countryData_tmp[country]
        countryData=countryData.fillna(0)           # Replace NaN with 0


                                    # Calculate fatality rates and clip to 100%
    # countryData['fatalityPercentage'] = countryData['deaths'] * 100./countryData['cases']
    # countryData['fatalityPercentage']=countryData['fatalityPercentage'].where(countryData['fatalityPercentage'] <= 100., 100.)
    # countryData.loc[countryData.cases == 0, "fatalityPercentage"] = 0                           # When cases 0= 0 set percentage to 0

    # countryData['fatalityPercentageCumulative'] = countryData['deathsCumulative'] * 100./countryData['casesCumulative']
    # countryData['fatalityPercentageCumulative']=countryData['fatalityPercentageCumulative'].where(countryData['fatalityPercentageCumulative'] <= 100., 100.)
    # countryData.loc[countryData.casesCumulative == 0, "fatalityPercentageCumulative"] = 0       # When cases 0= 0 set percentage to 0

    countryData = countryData.transpose()
    countryData = countryData.reset_index()

    outputFileName = "compare.csv"
    countryData.to_csv(outputFileName, index=True)

    return(countryData)


# fig, ax = plt.subplots(figsize=(15, 8))
fig, ax = plt.subplots(figsize=(10, 5.5))
def draw_bargraph(i, countryData, nDisplay, title):
    dateList = list(countryData)                    # Get dates (column headers)
    dateList.pop(-1)                                # Remove "country"
    date=dateList[i]                                # Date to process

    sorted = countryData.sort_values(by=[date], ascending=False)
    topN = sorted[['country', date]].copy()         # Extract column
    topN = topN.head(nDisplay)
    topN = topN.reset_index()

    # print (topN)                                    # Some debug
    # if (topN[date].iloc[0] > 0):                    # If the largest value is > zero, print it
        # print('date: ' + date + ' country: ' + topN['country'].iloc[0] + ' value: ' + str(topN[date].iloc[0]))

    topN = topN.sort_values(by=[date], ascending=True)    # Flip - largest on top - need to do this for rotated bar graph

    ax.clear()
    ax.barh(topN['country'], topN[date], color='#228fcf')

    countryList = topN['country'].tolist()

    ml = countryData.max(axis = 0).tolist()         # Get largest value in data frame and use to scale x axis
    del ml[-1]                                      # Delete last entry (country)
    maximum = max(ml)

    ax.autoscale(tight=True)
    ax.set_title(title, size=18)
    ax.axis(xmax=maximum * 1.05)

                                                    # Iterate over the values to plot labels and values
    for i, (value, name) in enumerate(zip(topN[date], topN['country'])):
        if (value > 0) :
            ax.text(value, i, name+' ',       ha='right')  # name
            ax.text(value, i, ' '+str(value), ha='left')   # value


    displayYAxisLabels = False                      # Set to True to displaying y axis labels
    for x in range(nDisplay):                       # If the count is zero OR we are not displaying the y-axis labels then clear the countryList entry
        if (topN[date].iloc[x] <= 0) or (displayYAxisLabels == False):
            countryList[x] = ""
    countryIndex = np.arange(nDisplay)
    ax.set_yticks(countryIndex, minor=False)
    ax.axes.set_yticklabels(countryList)

    ax.text(1, 0.1, date + '    ', transform=ax.transAxes, size=12, ha='right') # Add date

    return ax,


# main
def main(nDisplay, useCachedFileFlag, cumulativeResultsFlag, noPlotFlag, fileSavePlotFlag, popNormalizeFlag, deathsResultsFlag):

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

        for link in soup.find_all(class_="btn btn-primary", href=True):
        # for link in soup.find_all('a', href=True):
            # print(link['href'])
            if ("csv" in link['href']):
                csvfileurl = link['href']
                urllib.request.urlretrieve(csvfileurl, localFileName)
                cachedFilePresentFlag = True
                print("Cached file updated")
                break
            elif ("xlsx" in link['href']):          # If data in .xlsx format then retrieve and store as local .csv format
                xlsxfileurl = link['href']
                xlsx_tmp = pd.read_excel(xlsxfileurl, index_col=0)
                xlsx_tmp.to_csv(localFileName, index=True)
                cachedFilePresentFlag = True
                print("Cached file updated")
                break

        if (cachedFilePresentFlag == False):
            print("No spreadsheet found at the URL, use \"-l\" to use local cached file")
            exit(0)


    if (cachedFilePresentFlag == True):
        covidData = pd.read_csv(localFileName, index_col=0, encoding="iso8859_1")
                    # Spreadsheet columns :
                    # dateRep	day	month	year	cases	deaths	countriesAndTerritories	geoId	countryterritoryCode	popData2018

        covidData=covidData.fillna(0)               # Replace NaN with 0

        clen = 0                                    # For longest sequency
        dlen = 0

                                                    # Extract Chinese dates to create index - this allows for countries that do not have full data supplied
        dates_tmp = covidData[covidData["countriesAndTerritories"].str.match("China")]
        dates_tmp = dates_tmp.iloc[::-1]            # Invert table - top to bottom
        dates_tmp=dates_tmp.reset_index()
        dates=list(dates_tmp['dateRep'])

        if (cumulativeResultsFlag == True):
            if (deathsResultsFlag == True):
                countryData = extractCountries(covidData, "deathsCumulative", dates, popNormalizeFlag)
            else:
                countryData = extractCountries(covidData, "casesCumulative", dates, popNormalizeFlag)
        else:
            if (deathsResultsFlag == True):
                countryData = extractCountries(covidData, "deaths", dates, popNormalizeFlag)
            else:
                countryData = extractCountries(covidData, "cases", dates, popNormalizeFlag)

        # print(countryData)
        countryData = countryData.rename(columns={"index": "country"})

        countryColumn = countryData.pop('country')  # Remove country column
        countryData['country']=countryColumn        # Insert at the end of the dataframe

        countryData = countryData.loc[:, '24/01/2020':] # Remove columns up to the given date

        # print(countryData)

        dateList = list(countryData)                # Get dates (column headers)
        dateList.pop(-1)                            # Remove "country"

        numDates = len(dateList)
        print('numDates: ' + str(numDates))

                                                    # Get last date in combinedCases
        lastDate = str(covidData.first_valid_index())
        lastDate = lastDate.replace(' 00:00:00','')

        titleStr='Covid-19 '
        if (popNormalizeFlag == True):
            titleStr=titleStr + 'Pop Noramlized '
        if (cumulativeResultsFlag == True):
            titleStr = titleStr + "Cumulative "
        else:
            titleStr = titleStr + "Daily "
        if (deathsResultsFlag == True):
            titleStr = titleStr + "Deaths "
        else:
            titleStr = titleStr + "Cases "

        animator = animation.FuncAnimation(fig, draw_bargraph, frames=numDates, fargs=(countryData, nDisplay, titleStr + str(lastDate),), interval=250, blit=False)

        if (noPlotFlag == False):
            plt.show()

        if (fileSavePlotFlag == True):
            # HTML(animator.to_jshtml())
            # animator.to_html5_video()
            # animator.save(titleStr.strip()+'.mp4')

            animator.save(titleStr.replace(" ", "")+'BarGraph.gif', writer='imagemagick')#, dpi=80, fps=10)

    else:
        print("Cached spreadsheet file not found on computer")
        exit()


if __name__ == '__main__':
    nDisplay                = 10
    useCachedFileFlag       = False
    deathsResultsFlag       = False
    cumulativeResultsFlag   = False
    noPlotFlag              = False
    fileSavePlotFlag        = False
    popNormalizeFlag        = False


    parser = argparse.ArgumentParser(description='Covid-19 Comparison')
    parser.add_argument("-n", "--number",         nargs='?', const=10, help="Number of countries to display")
    parser.add_argument("-c", "--cumulative",     action="store_true", help="Display cumulative results")
    parser.add_argument("-l", "--local",          action="store_true", help="Use local cached Covid-19.csv")
    parser.add_argument("-d", "--deaths",         action="store_true", help="Deaths")
    parser.add_argument("-m", "--mortality_rate", action="store_true", help="Mortality rate")
    parser.add_argument("-f", "--file",           action="store_true", help="Save plot to file")
    parser.add_argument("-p", "--population",     action="store_true", help="Use population to normalize data to cases per 10 Million")
    parser.add_argument("-q", "--quiet",          action="store_true", help="Quiet - Do not plot graphs")
    args = parser.parse_args()

    if (args.local):
        useCachedFileFlag = True
        print("Use cached file = True")

    if (args.number):
        if (args.number != 1):
            nDisplay = int(args.number)
        print("Number to display =" + str(nDisplay))

    if (args.cumulative):
        cumulativeResultsFlag = True
        print("Cumulative Results = True")

    if (args.deaths):
        deathsResultsFlag = True
        print("Deaths Results = True")

    if (args.file):
        fileSavePlotFlag = True
        print("Save plot graphs to file = True")

    if (args.population):
        popNormalizeFlag = True
        print("Normalize to population = True")

    if (args.quiet):
        noPlotFlag = True
        print("Do not plot graphs = True")

    main(nDisplay, useCachedFileFlag, cumulativeResultsFlag, noPlotFlag, fileSavePlotFlag, popNormalizeFlag, deathsResultsFlag)
