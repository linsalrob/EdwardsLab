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
from random import randint, choice

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


headers = {
    'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36',
    'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
    'ACCEPT-ENCODING': 'gzip, deflate, br',
    'ACCEPT-LANGUAGE': 'en-AU,en;q=0.9,en-US;q=0.8',
}

# we use a requests session to store cookies etc
session = requests.Session()

def url_to_cites(url, outfile=None, verbose=False):
    """
    convert a url to a tuple of title and citation
    """
    r = session.get(url, headers=headers)
    if r.status_code != 200:
        print(f"FATAL: There was an error retrieving {url}. Return code was {r.status_code} and I can't continue", file=sys.stderr)
        return
    if verbose:
        message(f"Retrieved {url}\nParsing", "BLUE")
    if outfile:
        with open(outfile, 'w') as out:
            out.write(r.text)
    soup = BeautifulSoup(r.text, 'lxml')
    try:
        tit = soup.find("div", {"class": "gs_ri"}).find("h3", {"class": "gs_rt"}).text
    except AttributeError as e:
        return ("UNABLE TO RETRIEVE", "0")
    
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
    parser.add_argument('-p', help='partial output, these will be skipped')
    parser.add_argument('-s', help='random sleep time (seconds). Default=100', type=int, default=100)
    parser.add_argument('-w', help='write the html to specified directory')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if os.path.exists(args.o):
        message(f"{args.o} ALREADY EXISTS, NOT OVERWRITING", "RED");
        sys.exit(10)



    skip = set()
    if args.p:
        with open(args.p, 'r') as f:
            for l in f:
                skip.add(l.strip().split("\t")[0])

    if args.v:
        message(f"Parsing {args.f}", "GREEN")
    if args.w:
        os.makedirs(args.w, exist_ok=True)

    words = set()
    s = re.compile(r'cluster=(\d+)')
    gsid = re.compile(r'^\d+$')
    with open(args.f, 'r') as f, open(args.o, 'w') as out:
        for url in f:
            url = url.strip()

            if gsid.match(url):
                # this is just a gs cluster id
                cid = url
                url = f"https://scholar.google.com/scholar?cluster={cid}&hl=en&as_sdt=0,5"
            else:
                cid = s.findall(url)[0]
            if cid in skip:
                continue
            outfile = None
            if args.w:
                outfile = os.path.join(args.w, f"{cid}.html")
            tit, cit = url_to_cites(url, outfile, args.v)
            print(f"{cid}\t{tit}\t{cit}", file=out)
            time.sleep(randint(0,args.s))
            # here we search for a random word from all our titles so far
            words.update(filter(lambda x: len(x) > 5, tit.split(" ")))
            srch = choice(list(words))
            srurl = f"https://scholar.google.com.au/scholar?hl=en&as_sdt=0%2C5&q={srch}&btnG="
            r = session.get(url, headers=headers)
            if r.status_code == 429:
                message("FATAL: Google thinks we are a robot")
                sys.exit(1)
            if r.status_code != 200:
                message(
                    f"FATAL: There was an error retrieving {url}. Return code was {r.status_code} and I can't continue",
                    "RED")
                sys.exit(1)

