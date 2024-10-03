"""
Given a list of pubmed IDs in a file, create a CSV from the JSON data
"""
import csv
import os
import sys
import argparse
import requests
import xmltodict
import json

__author__ = 'Rob Edwards'

def fetch_pubmed_citations_with_abstracts(pmids):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "rettype": "abstract"
    }
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # Check for request errors
    return xmltodict.parse(response.content)



def fetch_pubmed_citations(pmids):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json"
    }
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # Check for request errors
    return response.json()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pubmed to csv')
    parser.add_argument('-f', help='input file')
    parser.add_argument('-p', help='pubmed ID', action='append')
    parser.add_argument('-n', help='number to get at once (default=100)', type=int, default=100)
    parser.add_argument('-o', help='output file (default = sys.stdout)')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.p:
        if not args.f:
            print("No arguments provided. Using our test case, pubmed ID: 31285584", file=sys.stderr)
            args.p = ['31285584', '25058116']
        else:
            args.p = []

    if args.f:
        with open(args.f, 'r') as f:
            for l in f:
                args.p.append(l.strip())

    rows = []
    for pmid_list in [args.p[i:i+args.n] for i in range(0, len(args.p), args.n)]:
        # Fetch citation information
        citations = fetch_pubmed_citations_with_abstracts(pmid_list)

        # Extract relevant fields
        articles = citations['PubmedArticleSet']['PubmedArticle']
        if isinstance(articles, dict):
            articles = [articles]

        for article in articles:
            medline_citation = article['MedlineCitation']
            pmid = medline_citation['PMID']['#text']
            article_info = medline_citation['Article']
            journal = article_info['Journal']['Title']
            article_title = article_info['ArticleTitle']
            abstract_text = article_info.get('Abstract', {}).get('AbstractText', "")

            if args.v:
                print(abstract_text, file=sys.stderr)

            if isinstance(abstract_text, dict) and '#text' in abstract_text:
                abstract_text = abstract_text['#text']

            # Handle multiple abstract sections if present
            if isinstance(abstract_text, list):
                if isinstance(abstract_text[0], dict):
                    abstract_text = [i['#text'] for i in abstract_text]
                try:
                    abstract_text = " ".join(abstract_text)
                except:
                    print(abstract_text, file=sys.stderr)
                    sys.exit(1)


            row = [pmid, journal, article_title, abstract_text]
            rows.append(row)

    # Define the CSV file header
    header = ['PMID', 'Journal', 'Article Title', 'Abstract']

    # Write to CSV file
    if args.o:
        csvfile = open(args.o, 'w', newline='', encoding='utf-8')
    else:
        csvfile = sys.stdout

    csvwriter = csv.writer(csvfile, delimiter='\t')
    csvwriter.writerow(header)
    csvwriter.writerows(rows)

    if args.o:
        csvfile.close()