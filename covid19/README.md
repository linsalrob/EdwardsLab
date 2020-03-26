# Graph the COVID-19 outbreak data from the European Union

## Authorship

This code was written by John Edwards. Rob might make some comments at some point

# COVID-19 Data

There are a lot of graphs of the COVID-19 outbreaks, but this code uses the data from the [European Centre for Disease Prevention and Control](https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide).

# Installation

You will need:

- pandas
- matplotlib
- numpy
- beautifulsoup4

# Running the code

Simply call it with your flavor of python and you will see a plot of the latest data:

```bash
python3 nCoV-Viz.py
```

The downloaded spreadsheet is stored locally in Covid-19.xlsx
To use cached local spreadsheet, use "-l" option

Intermediate data for Cases/Deaths and also for each country are stored in relevant .csv files

All plots can be aligned to :
  First date of detection or death, in that country (default)
  First date of detection in China, 2019-12-31 (-n)

Data can be plotted as daily values (default) cumulative values (-c)

Countries to plot and line colours are specified in the appropriate tables at the top of this file

