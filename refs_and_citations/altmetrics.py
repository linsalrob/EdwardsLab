"""
a list of altmetrics links, we'll go ahead and parse some information about them

We expect a link, but will check for it
"""

import os
import sys
import argparse
from roblib import message
from roblib_tk import choose_a_file, write_a_file
import requests
import re
from bs4 import BeautifulSoup
from time import sleep
from random import randint


__author__ = 'Rob Edwards'


class AltmetricsReference():
    def __init__(self, number, alt_id, url):
        self.count = number
        self.altmetric_id = alt_id
        self.url = url
        self.title = None
        self.summaries = []
        self.mentions = {}
        self.citations = {}
        self.age_and_source = ""
        self.age = ""
        self.in_the_top = ""


def parse_altmetrics_files(filename, verbose=False):
    s = re.compile('\w+$')
    urlsearch = re.compile('www.altmetric.com/details/(\d+)')

    headers = {
        'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.110 Safari/537.36'
    }

    counter  = 0
    all_references = []
    max_summaries = 0
    all_mention_sources = set()
    all_citation_sources = set()
    os.makedirs('altmetrics_html', exist_ok=True)
    with open(filename, 'r', encoding="utf-8") as f:
        for l in f:
            url = ""
            altid = ""
            l = l.strip()
            if not l:
                continue
            if s.match(l):
                altid = l
                url = f"https://www.altmetric.com/details/{l}"
            else:
                url = l
                m = urlsearch.search(l)
                if m:
                    altid = m.groups()[0]
            counter += 1
            altref = AltmetricsReference(counter, altid, url)

            if verbose:
                sys.stderr.write(f"Getting {counter}: ID: {altid} URL: {url}\n")

            if os.path.exists(os.path.join("altmetrics_html", f"{altid}.html")):
                with open(os.path.join("altmetrics_html", f"{altid}.html"), 'r', encoding="utf-8") as f:
                    html = ""
                    for l in f:
                        html += l
                    soup = BeautifulSoup(html, 'lxml')
            else:
                sleep(randint(1,5))
                r = requests.get(url, headers=headers)
                if r.status_code != 200:
                    sys.stderr.write(f"ERROR: There was an error retrieving {url} and I can't continue\n")
                soup = BeautifulSoup(r.text, 'lxml')
                with open(os.path.join("altmetrics_html", f"{altid}.html"), "w", encoding="utf-8") as htmlout:
                    htmlout.write(r.text)
            if not soup:
                sys.stderr.write(f"Skipping {altid}\n")
                continue

            # get the title
            dh = soup.find('div', class_="document-header").find("a")
            if dh:
                for c in dh:
                    if c.text:
                        altref.title = c.text
                print(f"Title: {altref.title}")

            # get the summaries
            summ = soup.find_all('div', class_='summary')
            if summ:
                for c in summ:
                    if 'compared to outputs of the same age and source' in c.text:
                        altref.age_and_source = c.text
                    elif 'compared to outputs of the same age' in c.text:
                        altref.age = c.text
                    elif 'In the top' in c.text:
                        altref.in_the_top = c.text
                    else:
                        altref.summaries.append(c.text)
                if len(altref.summaries) > max_summaries:
                    max_summaries = len(altref.summaries)

            ex = re.compile('(\d+)\s+(\S.*)$')
            # extract mentions
            mc = soup.find_all('dl', class_='mention-counts')
            if mc:
                for c in mc:
                    for lnk in c.find_all('a'):
                        m = ex.match(lnk.text)
                        if m:
                            mentionby = m.groups()[1]
                            if mentionby.endswith('s'):
                                mentionby = mentionby[:-1]
                            altref.mentions[mentionby] = m.groups()[0]
                            all_mention_sources.add(mentionby)

            # extract citations
            mc = soup.find_all('dl', class_='scholarly-citation-counts')
            if mc:
                for c in mc:
                    for lnk in c.find_all('a'):
                        m = ex.match(lnk.text)
                        if m:
                            altref.citations[m.groups()[1]] = m.groups()[0]
                            all_citation_sources.add(m.groups()[1])
            all_references.append(altref)

    return all_references, max_summaries, all_mention_sources, all_citation_sources

def print_references(all_references, max_summaries, all_mention_sources, all_citation_sources, outputfile):
    out = open(outputfile, 'w', encoding="utf-8")
    print("Number\tAltmetric ID\tAltmetric URL\tTitle\tTop\tAge and Source\tAge", end="", file=out)
    for i in range(max_summaries):
        print(f"\tSummary {i+1}", end="", file=out)
    acs = sorted(list(all_citation_sources))
    for s in acs:
        print(f"\tCitations ({s})", end="", file=out)
    ams = sorted(list(all_mention_sources))
    for s in ams:
        print(f"\tMentions ({s})", end="", file=out)
    print(file=out)


    for altref in all_references:
        data = [
            altref.count,
            altref.altmetric_id,
            altref.url,
            altref.title,
            altref.in_the_top,
            altref.age_and_source,
            altref.age
        ]
        for s in sorted(altref.summaries):
            data.append(s)
        for i in range(max_summaries - len(altref.summaries)):
            data.append("")
        for s in acs:
            if s in altref.citations:
                data.append(altref.citations[s])
            else:
                data.append("")
        for s in ams:
            if s in altref.mentions:
                data.append(altref.mentions[s])
            else:
                data.append("")

        print("\t".join(map(str, data)), file=out)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    verbose = True
    if args.f:
        filename = args.f
        verbose = args.v
    else:
        filename = choose_a_file("File of altmetric links")
        verbose = True

    if args.o:
        outputfile = args.o
    else:
        outputfile = write_a_file("tsv file to write")

    all_references, max_summaries, all_mention_sources, all_citation_sources = parse_altmetrics_files(filename, verbose)
    print_references(all_references, max_summaries, all_mention_sources, all_citation_sources, outputfile)
