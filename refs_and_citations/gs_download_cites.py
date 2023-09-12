"""
Given a list of GS URLs, download the page from google scholar and get the title and number  of citations

We use requests to get the page, and then sleep a random time to try and avoid bot detection!
"""
import argparse
import os
from roblib import message
import sys
import requests
import time
from bs4 import BeautifulSoup
import re
from random import randint

__author__ = 'Rob Edwards'


# set the list of URLs here, one per line. I use the cluster URL
urls = """https://scholar.google.com.au/scholar?cluster=13011484634294850749&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=15948142415166130473&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=4483447249227580912&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=11513119306974836518&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=4513723272904109625&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=12765392127648450959&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=16106135173513489063&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=8144540207175036240&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=15312967259686395118&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=3453650122218310673&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=10194930186593133536&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=2653269953959691334&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=7067622931716524932&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=5911286658349579584&hl=en&as_sdt=0,5
https://scholar.google.com.au/scholar?cluster=4555724394641856429&hl=en&as_sdt=0,5
"""




def url_to_cites(url, verbose=False):
    """
    convert a url to a tuple of title and citation
    """
    headers = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.110 Safari/537.36'}
    r = requests.get(url, headers=headers)
    if r.status_code != 200:
        print(f"FATAL: There was an error retrieving {url}. Return code was {r.status_code} and I can't continue", file=sys.stderr)
        return
    if verbose:
        message(f"Retrieved {url}\nParsing", "BLUE")
    soup = BeautifulSoup(r.text, 'lxml')
    tit = soup.find("div", {"class": "gs_ri"}).find("h3", {"class": "gs_rt"}).text
    if verbose:
        message(f"Title: {tit}", "BLUE")
    cites = [0]
    cb = re.compile('Cited by (\d+)')
    for c in cb.findall(str(soup)):
        cites.append(int(c))
    if verbose:
        message(f"All citations: {cites}. We used the most!", "BLUE")
    best_cite = sorted(cites, reverse=True)[0]

    return (tit, best_cite)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='file of one gs link per line', required=True)
    parser.add_argument('-o', help='output tsv file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.v:
        message(f"Parsing {args.f}", "GREEN")

    with open(args.f, 'r') as f, open(args.o, 'w') as out:
        for url in f:
            url = url.strip()
            tit, cit = url_to_cites(url, args.v)
            print(f"{tit}\t{cit}", file=out)
            time.sleep(randint(0,60))

