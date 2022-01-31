"""
Get and parse one url
"""

import os
import sys
import requests
from bs4 import BeautifulSoup
import re

__author__ = 'Rob Edwards'

if True:
    url = 'https://www.altmetric.com/details/1386636'
    headers = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.110 Safari/537.36'}
    r = requests.get(url, headers=headers)
    if r.status_code != 200:
        sys.stderr.write(f"FATAL: There was an error retrieving {url} and I can't continue\n")
    soup = BeautifulSoup(r.text, 'lxml')
else:
    print("Reading file")
    text = ""
    with open("1386636.html", 'r') as f:
        for l in f:
            text += l
    soup = BeautifulSoup(text, 'lxml')
    print("Done")

dh = soup.find('div', class_="document-header").find("a")
for c in dh:
    if c.text:
        title = c.text
print(f"Title: {title}")

summ = soup.find_all('div', class_='summary')
summaries = []
for c in summ:
    summaries.append(c.text)
print(summaries)

ex = re.compile('(\d+)\s+(\S.*)$')
# extract mentions
mentions = {}
mc = soup.find_all('dl', class_='mention-counts')
for c in mc:
    for lnk in c.find_all('a'):
        m = ex.match(lnk.text)
        if m:
            mentions[m.groups()[1]] = m.groups()[0]

print(mentions)
# extract citations
citations = {}
mc = soup.find_all('dl', class_='scholarly-citation-counts')
for c in mc:
    for lnk in c.find_all('a'):
        m = ex.match(lnk.text)
        if m:
            citations[m.groups()[1]] = m.groups()[0]
print(citations)