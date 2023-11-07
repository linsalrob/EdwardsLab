"""
Use the openAPI to get some genome information
"""
import json
import os
import sys
import argparse
import requests


__author__ = 'Rob Edwards'

if True:
    NCBI_API_KEY = os.environ['NCBI_API_KEY']
    rest_url = "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml"
    url = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha'
    accession = 'GCA_008244535.1,GCA_019119935.1'

    summary = f"/genome/accession/{accession}/dataset_report"

    r = requests.get(url + summary, headers={"Content-Type": "text", "api-key": NCBI_API_KEY})
    with open("temp.txt", 'w') as out:
        out.write(r.text)
    d = json.loads(r.text)
else:
    with open("temp.txt", "r") as f:
        d = json.load(f)



data = {}

report_keys = ['accession', 'current_accession', 'paired_accession', 'source_database', 'organism', 'assembly_info',
               'assembly_stats', 'annotation_info', 'wgs_info', 'checkm_info', 'average_nucleotide_identity']

wanted_keys = ['organism', 'assembly_info',
               'assembly_stats', 'annotation_info', 'wgs_info', 'checkm_info', 'average_nucleotide_identity']




allkeys = set()

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



thekeys = sorted(allkeys)

with open('output.tsv', 'w') as out:
    print("\t".join(["accession"] + thekeys), file=out)
    for acc in data:
        print(acc, end="", file=out)
        for key in thekeys:
            print("\t" + str(data[acc].get(key, "")).replace('\n', ' '), end="", file=out)
        print(file=out)

"""

for acc in data:
    for key in thekeys:
        val = str(data[acc].get(key, ''))
        print(f"{key}\t{val}")
    print()


"""