"""
Does every sequence in the crAss001 genome, is there a gp number
"""

import os
import sys
from roblib import genbank_seqio, feature_id

__author__ = 'Rob Edwards'

gbkf = "C:\\Users\\edwa0468\\Downloads\\crass001.gbk"

out = open("C:\\Users\\edwa0468\\Downloads\\crass001.faa", "w")
for seq in genbank_seqio(gbkf):
    for feat in seq.features:
        if feat.type != 'CDS':
            continue
        nt = feat.qualifiers["note"]
        if 'Gp' not in nt[0]:
            print(f"{nt}")
        (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)
        prtmtd = {
            'EC_number': "",
            'locus_tag': "",
            'old_locus_tag': "",
            'note': "",
            'product': "",
            'protein_id': "",
            'ribosomal_slippage': "",
            'transl_table': 11,
            'translation': ""
        }

        cid = feature_id(seq, feat)


        loc = f"{start}_{stop}"
        if strand < 0:
            loc = f"{stop}_{start}"

        cid += f' [{seq.id}] '
        if 'organism' in seq.annotations:
            cid += f' [{seq.annotations["organism"]}]'
        cid += f' [{seq.id}_{loc}]'
        if 'product' in feat.qualifiers:
            cid += f' [{feat.qualifiers["product"][0]}]'
        else:
            cid += f' [hypothetical protein]'
        lt = None
        if 'old_locus_tag' in feat.qualifiers:
            lt = feat.qualifiers['old_locus_tag'][0]
            lt = lt.replace('crAss001_', '')
            if lt.startswith('t'):
                lt = f"g{lt}"
            else:
                lt = f"gp{lt}"

        if 'note' in feat.qualifiers:
            cid += " [" + " ".join(feat.qualifiers['note']) + "]"
        else:
            print(f"No note in {cid}", file=sys.stderr)

        transp = ""
        if 'translation' in feat.qualifiers:
            transp = feat.qualifiers['translation'][0]
        else:
            transp = str(feat.extract(seq).translate().seq)

        if lt:
            print(f">{lt}\n{transp}", file=out)
        else:
            print(f">{cid}\n{transp}", file=out)


out.close()