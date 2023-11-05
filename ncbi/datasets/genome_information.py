"""
Get genome information in bulk. See one_genome_information for debugging purposes!
"""
import json
import os
import sys
import argparse
import requests
from roblib import colours

__author__ = 'Rob Edwards'


def get_accessions(allkeys, data, accessions):
    """
    Get the data for some accesssions
    """
    NCBI_API_KEY = os.environ['NCBI_API_KEY']
    rest_url = "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml"
    url = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha'

    summary = f"/genome/accession/{accessions}/dataset_report"

    r = requests.get(url + summary, headers={"Content-Type": "text", "api-key": NCBI_API_KEY})
    with open("temp.txt", 'w') as out:
        out.write(r.text)
    d = json.loads(r.text)

    wanted_keys = ['organism', 'assembly_info',
                   'assembly_stats', 'annotation_info', 'wgs_info', 'checkm_info', 'average_nucleotide_identity']

    for report in d['reports']:
        accession = report['current_accession']
        data[accession] = {}

        for wk in wanted_keys:
            if wk not in report:
                continue
            for item_key in report[wk].keys():
                if wk == 'assembly_info' and item_key == 'biosample':
                    continue
                keyname = f"{wk} : {item_key}"
                data[accession][keyname] = report[wk][item_key]
                allkeys.add(keyname)

        wk = 'assembly_info'
        for item_key in ['biosample']:
            for sub_key in report[wk][item_key].keys():
                if sub_key == 'attributes':
                    for kvpair in report[wk][item_key][sub_key]:
                        keyname = f"{wk} : {item_key} : {sub_key} : {kvpair['name']}"
                        data[accession][keyname] = kvpair['value']
                        allkeys.add(keyname)
                else:
                    keyname = f"{wk} : {item_key} : {sub_key}"
                    data[accession][keyname] = report[wk][item_key][sub_key]
                    allkeys.add(keyname)

    return allkeys, data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-a', '--accessions', help='file of GenBank assembly accessions', required=True)
    parser.add_argument('-n', '--num_requests', help='number of simultaneous requests default: 1000', default=1000, type=int)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if 'NCBI_API_KEY' not in os.environ:
        print(f"{colours.RED}FATAL: Please set the environment variable NCBI_API_KEY by ", file=sys.stderr)
        print(f"'export NCBI_API_KEY=xxxxx'. You can get the key from your NCBI profile.{colours.ENDC}")
        sys.exit(2)

    allkeys = set()
    data = {}

    with open(args.accessions, 'r') as f:
        acc_count = 0
        accessions = None
        for l in f:
            l = l.strip()
            if acc_count > 0:
                accessions = f"{accessions},{l}"
            else:
                accessions = l
            acc_count += 1
            if acc_count >= args.num_requests:
                allkeys, data = get_accessions(allkeys, data, accessions)
                accessions = None
                acc_count = 0

    thekeys = sorted(allkeys)

    with open('output.tsv', 'w') as out:
        print("\t".join(["accession"] + thekeys), file=out)
        for acc in data:
            print(acc, end="", file=out)
            for key in thekeys:
                print("\t" + str(data[acc].get(key, "")).replace('\n', ' '), end="", file=out)
            print(file=out)