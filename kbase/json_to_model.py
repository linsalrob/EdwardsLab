"""
Link functions in a JSON file from KBase to model data using PyFBA
"""

import os
import sys
sys.path.insert(0, '/data/PyFBA')

import re
import argparse
import PyFBA
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse the JSON file downloaded from KBase")
    parser.add_argument('-f', help='JSON file', required=True)
    args = parser.parse_args()

    data = json.load(open(args.f))


    # print keys associated with features
    feats = data['features']

    # parse the plant data from the model seed
    compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes('plant')

    # find out which enzymes have EC numbers
    ecnos = {}
    allroles = {}
    for e in enzymes:
        for ec in enzymes[e].ec_number:
            if ec in ecnos:
                ecnos[ec].add(e)
            else:
                ecnos[ec] = {e}
        for r in enzymes[e].roles:
            if r.lower() in allroles:
               allroles[r.lower()].add(e)
            else:
                allroles[r.lower()] = {e}



    cds=0
    allrids = set()
    for f in feats:
        if 'CDS' in f['id']:
            cds+=1
        for p in f['function'].split(' ; '):
            p = p.strip()
            m = re.match('\s*[\d\-\.]+$', p)
            if m and m.end() == len(p):
                # p is an ec number
                if p in ecnos:
                    for e in ecnos[p]:
                        for rid in enzymes[e].reactions:
                            # print("\t".join([f['id'], 'EC ' + p, rid]))
                            allrids.add(rid)
                # else:
                #     print("NOT FOUND |" + p + "|")
            if p.lower() in allroles:
                for e in allroles[p.lower()]:
                    for rid in enzymes[e].reactions:
                        allrids.add(rid)



    print("There are " + str(cds) + " coding features that map to " + str(len(allrids)) + " reaction ids")